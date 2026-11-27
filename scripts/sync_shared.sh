#!/usr/bin/env bash
# Copy shared/python/scio.py into every module that carries it (modules/cellgeni/<tool>/<step>/resources/usr/bin/).
# Modules are published to the registry on their own, so each one ships its own copy; edit only the shared file.
#
#   scripts/sync_shared.sh                    update every copy
#   scripts/sync_shared.sh <module dir> ...   also add a copy to these modules
#   scripts/sync_shared.sh --check            fail when a copy differs from the shared file (CI)
set -euo pipefail
cd "$(dirname "$0")/.."

shared=shared/python/scio.py
check=false
copies=()
for arg in "$@"; do
    case "$arg" in
        --check) check=true ;;
        *) [[ -f "$arg/main.nf" ]] || { echo "not a module directory: $arg" >&2; exit 2; }
           copies+=("$arg/resources/usr/bin/scio.py") ;;
    esac
done
while IFS= read -r copy; do copies+=("$copy"); done < <(find modules -path '*/resources/usr/bin/scio.py' | sort)

drift=0
for copy in $(printf '%s\n' "${copies[@]}" | sort -u); do
    if cmp -s "$shared" "$copy"; then
        continue
    elif $check; then
        echo "out of date: $copy (run scripts/sync_shared.sh)" >&2
        drift=1
    else
        mkdir -p "$(dirname "$copy")"
        cp "$shared" "$copy"
        echo "updated $copy"
    fi
done
exit $drift
