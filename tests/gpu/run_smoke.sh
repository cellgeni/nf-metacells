#!/usr/bin/env bash
# GPU checks on the Sanger farm: one LSF job per GPU model runs, with Singularity and --nv,
#   smoke_rapids.py and smoke_seacells.py   the GPU images work on the farm driver
#   parity_rna.py and parity_atac.py        the step scripts' GPU paths match their CPU paths
# Run after changing a GPU image or a step script's GPU path. Each check prints one SMOKE_REPORT JSON line.
#
#   tests/gpu/run_smoke.sh [gmodel ...]     default: NVIDIAA100_SXM4_80GB NVIDIAH10080GBHBM3
#   logs: ${SMOKE_LOGS:-.dev/gpu_smoke}/<jobid>.<gmodel>.out
set -euo pipefail
cd "$(dirname "$0")/../.."

models=("${@:-NVIDIAA100_SXM4_80GB NVIDIAH10080GBHBM3}")
cache=/nfs/cellgeni/singularity/images
logs=${SMOKE_LOGS:-.dev/gpu_smoke}
mkdir -p "$logs"

# The image file in the Nextflow cache, pulled under the name Nextflow gives it if missing
sif() {
    local file
    file=$cache/$(tr '/:' '--' <<< "$1").img
    [[ -f "$file" ]] || singularity pull "$file" "docker://$1" >&2
    echo "$file"
}
rapids=$(sif quay.io/cellgeni/rapids-singlecell:0.18.1)
seacells=$(sif quay.io/cellgeni/seacells:a1173c8-gpu)
snapatac2=$(sif community.wave.seqera.io/library/python_pip_gcc_conda-gcc-specs_pruned:b8199b0e3c90fd81)
# Wave images set no environment: keep ~/.local out and the numba and matplotlib caches writable
run="singularity exec --nv -B /lustre,/nfs --env PYTHONPATH=$PWD/shared/python,PYTHONNOUSERSITE=1"
run+=" --env NUMBA_CACHE_DIR=/tmp/numba,MPLCONFIGDIR=/tmp/matplotlib"
atac="$PWD/$logs/parity_atac.\$LSB_JOBID"
checks=(
    "$run $rapids python tests/gpu/smoke_rapids.py"
    "$run $seacells python tests/gpu/smoke_seacells.py"
    "$run $rapids python tests/gpu/parity_rna.py"
    "$run $snapatac2 python tests/gpu/parity_atac.py prepare $atac && $run $rapids python tests/gpu/parity_atac.py compare $atac"
)

for model in ${models[@]}; do
    bsub -q gpu -n 4 -M 64000 -R "select[mem>64000] rusage[mem=64000] span[hosts=1]" \
        -gpu "mode=shared:j_exclusive=yes:num=1:gmodel=$model" \
        -J "nfmc_gpu_checks_$model" -o "$logs/%J.$model.out" -e "$logs/%J.$model.err" \
        "nvidia-smi --query-gpu=name,driver_version --format=csv,noheader; $(IFS=';'; echo "${checks[*]}")"
done
