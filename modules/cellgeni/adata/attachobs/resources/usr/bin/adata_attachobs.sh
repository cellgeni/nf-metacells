#!/usr/bin/env bash
# Attach cell metadata to an AnnData store's obs and/or make barcodes unique per sample, using adata-cli and duckdb.
#
#   adata_attachobs.sh --adata IN --output OUT --sample-id ID [--metadata CSV] [--barcode-column obs_names|COL]
#                      [--join inner|left] [--delimiter STR] [--sample-column NAME] [--allow-dtype-change]
#
# Metadata CSV: a `barcode` column plus any annotations. With a `sample_id` column, only the rows of --sample-id are
# used; without one, every row is. Rows are joined on --barcode-column (default: obs_names). An inner join drops
# cells without metadata; a left join keeps them.
# Metadata columns that already exist in obs are attached with a `_metadata` suffix.
# --delimiter renames obs_names to <barcode><delimiter><sample-id> (the join always uses the original names) and
# keeps the original names in a `barcode` column (`barcode_original` if obs already has `barcode`).
# --sample-column adds (or overwrites) a column holding --sample-id.
#
# obs is round-tripped through CSV because adata-cli has no join. adata-cli 0.17.6 turns categorical columns with
# digit-only labels into integers (cellgeni/adata-cli#124), so the script fails if any categorical column would
# change type, unless --allow-dtype-change is given. Column names containing commas are not supported.
set -euo pipefail

die() { echo "Error: $*" >&2; exit 1; }
warn() { echo "Warning: $*" >&2; }
sql_str() { printf "'%s'" "${1//\'/\'\'}"; }   # SQL string literal
sql_id() { printf '"%s"' "${1//\"/\"\"}"; }     # SQL identifier

ADATA="" OUTPUT="" METADATA="" SAMPLE_ID="" BARCODE_COLUMN="obs_names" JOIN="inner" DELIMITER="" SAMPLE_COLUMN=""
ALLOW_DTYPE_CHANGE=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --adata) ADATA="$2"; shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --metadata) METADATA="$2"; shift 2 ;;
        --sample-id) SAMPLE_ID="$2"; shift 2 ;;
        --barcode-column) BARCODE_COLUMN="$2"; shift 2 ;;
        --join) JOIN="$2"; shift 2 ;;
        --delimiter) DELIMITER="$2"; shift 2 ;;
        --sample-column) SAMPLE_COLUMN="$2"; shift 2 ;;
        --allow-dtype-change) ALLOW_DTYPE_CHANGE=true; shift ;;
        *) die "unknown argument: $1" ;;
    esac
done

[[ -n "$ADATA" && -n "$OUTPUT" && -n "$SAMPLE_ID" ]] || die "--adata, --output and --sample-id are required"
[[ -e "$ADATA" ]] || die "input store not found: $ADATA"
[[ -z "$METADATA" || -f "$METADATA" ]] || die "metadata file not found: $METADATA"
[[ "$JOIN" == inner || "$JOIN" == left ]] || die "--join must be inner or left, got '$JOIN'"
[[ -n "$METADATA" || -n "$DELIMITER" || -n "$SAMPLE_COLUMN" ]] \
    || die "nothing to do: give --metadata, --delimiter or --sample-column"

TMP=$(mktemp -d attachobs.XXXXXX)   # scratch files, kept apart from staged inputs
DB="$TMP/obs.duckdb"
sql() { duckdb -noheader -list "$DB" -c "$1"; }
categorical_columns() {  # categorical obs columns of a store, sorted, one per line
    adata view "$1" obs | sed -n 's/^ *\([^ ]*\) <categorical>.*/\1/p' | sort
}

# 1. Export obs and number its rows, so the original order survives the join
adata export dataframe "$ADATA" obs -o $TMP/obs.csv
HEADER=$(head -n 1 $TMP/obs.csv | tr -d '"')
INDEX_COLUMN=${HEADER%%,*}
mapfile -t OBS_COLUMNS < <(tr ',' '\n' <<< "$HEADER" | tail -n +2)
awk 'NR == 1 { print "__row," $0; next } { print NR - 1 "," $0 }' $TMP/obs.csv > $TMP/obs_rows.csv
sql "CREATE TABLE o AS SELECT * FROM read_csv('$TMP/obs_rows.csv', header = true, all_varchar = true);"
N_OBS=$(sql "SELECT count(*) FROM o;")
categorical_columns "$ADATA" > $TMP/categorical_before.txt

has_obs_column() { printf '%s\n' "${OBS_COLUMNS[@]}" | grep -qxF -- "$1"; }
IDX=$(sql_id "$INDEX_COLUMN")
KEY=$IDX
if [[ "$BARCODE_COLUMN" != obs_names ]]; then
    has_obs_column "$BARCODE_COLUMN" || die "barcode column '$BARCODE_COLUMN' not found in obs"
    KEY=$(sql_id "$BARCODE_COLUMN")
fi

# 2. Select list: new index, optional barcode/sample columns, existing obs columns, metadata columns
SELECT="o.$IDX AS $IDX"
BARCODE_TARGET=""
if [[ -n "$DELIMITER" ]]; then
    SELECT="o.$IDX || $(sql_str "$DELIMITER") || $(sql_str "$SAMPLE_ID") AS $IDX"
    BARCODE_TARGET=barcode
    if has_obs_column barcode; then
        BARCODE_TARGET=barcode_original
        warn "obs already has a 'barcode' column; original names kept in 'barcode_original'"
    fi
    SELECT+=", o.$IDX AS $(sql_id "$BARCODE_TARGET")"
fi
if [[ -n "$SAMPLE_COLUMN" ]]; then
    SELECT+=", $(sql_str "$SAMPLE_ID") AS $(sql_id "$SAMPLE_COLUMN")"
    has_obs_column "$SAMPLE_COLUMN" && warn "obs column '$SAMPLE_COLUMN' is overwritten with '$SAMPLE_ID'"
fi
for column in "${OBS_COLUMNS[@]}"; do
    [[ -n "$SAMPLE_COLUMN" && "$column" == "$SAMPLE_COLUMN" ]] && continue
    SELECT+=", o.$(sql_id "$column")"
done

FROM="FROM o"
if [[ -n "$METADATA" ]]; then
    sql "CREATE TABLE meta_all AS SELECT * FROM read_csv($(sql_str "$METADATA"), header = true, all_varchar = true);"
    mapfile -t META_COLUMNS < <(sql "SELECT column_name FROM (DESCRIBE meta_all);")
    printf '%s\n' "${META_COLUMNS[@]}" | grep -qx barcode || die "metadata must have a 'barcode' column"
    if printf '%s\n' "${META_COLUMNS[@]}" | grep -qx sample_id; then
        sql "CREATE TABLE m AS SELECT * FROM meta_all WHERE sample_id = $(sql_str "$SAMPLE_ID");"
        ROWS="rows for sample '$SAMPLE_ID'"
    else
        sql "CREATE TABLE m AS SELECT * FROM meta_all;"
        ROWS="rows (no sample_id column, all used)"
    fi

    N_META=$(sql "SELECT count(*) FROM m;")
    [[ "$N_META" -gt 0 ]] || die "no metadata rows with sample_id '$SAMPLE_ID'"
    DUPLICATES=$(sql "SELECT barcode FROM m GROUP BY barcode HAVING count(*) > 1 LIMIT 5;" | paste -sd' ' -)
    [[ -z "$DUPLICATES" ]] || die "duplicated barcodes in the metadata of sample '$SAMPLE_ID': $DUPLICATES"

    for column in "${META_COLUMNS[@]}"; do
        [[ "$column" == barcode || "$column" == sample_id ]] && continue
        target="$column"
        if has_obs_column "$column" || [[ "$column" == "$BARCODE_TARGET" || "$column" == "$SAMPLE_COLUMN" ]]; then
            target="${column}_metadata"
            warn "metadata column '$column' already exists in obs; attached as '$target'"
        fi
        SELECT+=", m.$(sql_id "$column") AS $(sql_id "$target")"
    done

    N_MATCHED=$(sql "SELECT count(*) FROM o JOIN m ON o.$KEY = m.barcode;")
    echo "Metadata: $N_META $ROWS; $N_MATCHED of $N_OBS cells matched; join=$JOIN" >&2
    [[ "$N_MATCHED" -gt 0 ]] || die "no cell matched the metadata (obs '${BARCODE_COLUMN}' vs metadata 'barcode')"
    if (( N_META > N_MATCHED )); then warn "$(( N_META - N_MATCHED )) metadata rows match no cell"; fi
    FROM="FROM o $([[ "$JOIN" == inner ]] && echo JOIN || echo 'LEFT JOIN') m ON o.$KEY = m.barcode"
fi

sql "COPY (SELECT $SELECT $FROM ORDER BY CAST(o.__row AS INTEGER)) TO '$TMP/obs_new.csv' (HEADER, DELIMITER ',');"

# 3. Inner join that drops cells: subset first (adata subset keeps the source order)
BASE="$ADATA"
if [[ -n "$METADATA" && "$JOIN" == inner && "$N_MATCHED" -lt "$N_OBS" ]]; then
    sql "SELECT o.$IDX FROM o JOIN m ON o.$KEY = m.barcode ORDER BY CAST(o.__row AS INTEGER);" > $TMP/keep.txt
    BASE="$TMP/subset.${OUTPUT##*.}"
    adata subset "$ADATA" -o "$BASE" --obs $TMP/keep.txt
fi

# 4. Write the new obs, asking adata-cli to keep the originally categorical columns categorical
CATEGORICAL=$(paste -sd, $TMP/categorical_before.txt)
adata import dataframe "$BASE" obs $TMP/obs_new.csv -o "$OUTPUT" ${CATEGORICAL:+--categorical "$CATEGORICAL"}

# 5. Guard against dtype changes from the CSV round trip (cellgeni/adata-cli#124)
categorical_columns "$OUTPUT" > $TMP/categorical_after.txt
CHANGED=$(comm -23 $TMP/categorical_before.txt $TMP/categorical_after.txt | paste -sd, -)
if [[ -n "$CHANGED" ]]; then
    MESSAGE="categorical obs columns lost their categorical type in the CSV round trip: $CHANGED (cellgeni/adata-cli#124)"
    if [[ "$ALLOW_DTYPE_CHANGE" == true ]]; then
        warn "$MESSAGE"
    else
        rm -rf "$OUTPUT"
        die "$MESSAGE. Pass --allow-dtype-change to accept the change."
    fi
fi
rm -rf "$TMP"
echo "Wrote $OUTPUT" >&2
