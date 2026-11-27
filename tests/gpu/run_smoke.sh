#!/usr/bin/env bash
# Smoke-test the GPU images on the Sanger farm: one LSF job per GPU model runs tests/gpu/smoke_*.py in each image
# (Singularity, --nv). Run after changing a GPU image; the report lines are SMOKE_REPORT JSON.
#
#   tests/gpu/run_smoke.sh [gmodel ...]     default: NVIDIAA100_SXM4_80GB NVIDIAH10080GBHBM3
#   logs: ${SMOKE_LOGS:-.dev/gpu_smoke}/<jobid>.<gmodel>.out
set -euo pipefail
cd "$(dirname "$0")/../.."

images=(quay.io/cellgeni/rapids-singlecell:0.18.1 quay.io/cellgeni/seacells:a1173c8-gpu)
scripts=(tests/gpu/smoke_rapids.py tests/gpu/smoke_seacells.py)
models=("${@:-NVIDIAA100_SXM4_80GB NVIDIAH10080GBHBM3}")
cache=/nfs/cellgeni/singularity/images
logs=${SMOKE_LOGS:-.dev/gpu_smoke}
mkdir -p "$logs"

# Pull into the Nextflow cache, under the name Nextflow gives the image
commands=()
for i in "${!images[@]}"; do
    sif=$cache/$(tr '/:' '--' <<< "${images[$i]}").img
    [[ -f "$sif" ]] || singularity pull "$sif" "docker://${images[$i]}"
    commands+=("singularity exec --nv -B /lustre,/nfs --env PYTHONPATH=$PWD/shared/python $sif python ${scripts[$i]}")
done

for model in ${models[@]}; do
    bsub -q gpu -n 4 -M 32000 -R "select[mem>32000] rusage[mem=32000] span[hosts=1]" \
        -gpu "mode=shared:j_exclusive=yes:num=1:gmodel=$model" \
        -J "nfmc_smoke_$model" -o "$logs/%J.$model.out" -e "$logs/%J.$model.err" \
        "nvidia-smi --query-gpu=name,driver_version --format=csv,noheader; $(IFS=';'; echo "${commands[*]}")"
done
