#!/bin/bash
# submit_mivia.sh - submit one algorithm for the MIVIA dataset
# Usage: bash hpc/submit_mivia.sh <algo> [--with-ref]
# --with-ref: also run reference binary and compare

ALGO=$1
WITH_REF=0
if [ "$2" = "--with-ref" ]; then
    WITH_REF=1
fi

if [ -z "$ALGO" ]; then
    echo "Usage: bash hpc/submit_mivia.sh <algo>"
    echo "Algorithms: mcsplit rl ll, ll_lsm, ll_lum, dal, dal_dal, dal_rl, dsb, dsb_always, dsb_never, rrsplit, rrsplit_noveq, rrsplit_nomax, rrsplit_nobound, symsplit, symsplit_valonly, symsplit_varonly"
    exit 1
fi

mkdir -p hpc/logs hpc/scripts hpc/results

LIST="hpc/instances_${ALGO}.txt"

# Generate instance list if it doesn't exist
if [ ! -f "$LIST" ]; then
    echo "Generating instance list for $ALGO..."
    > "$LIST"
    find instances/MIVIA -name "*.A0[0-9]" | sort | while read APATH; do
        BPATH="${APATH/.A/.B}"
        [ -f "$BPATH" ] && echo "$APATH $BPATH" >> "$LIST"
    done
    echo "Generated $(wc -l < $LIST) instance pairs"
fi

echo "MIVIA instances: $(wc -l < $LIST)"

NJOBS=$(($(wc -l < "$LIST") - 1))

if [ "$WITH_REF" = "1" ]; then
    # Validation mode: use reference comparison script
    SCRIPT="hpc/scripts/job_${ALGO}_ref.slurm"
    sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g" \
        hpc/scripts/job_array_ref.slurm > "$SCRIPT"
    sed -i 's/\r$//' "$SCRIPT"
    echo "Submitting $ALGO MIVIA (with reference comparison): $((NJOBS+1)) jobs"
else
    # Experiment mode: unified solver only
    SCRIPT="hpc/scripts/job_${ALGO}.slurm"
    sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g" \
        hpc/scripts/job_array.slurm > "$SCRIPT"
    sed -i 's/\r$//' "$SCRIPT"
    echo "Submitting $ALGO MIVIA: $((NJOBS+1)) jobs"
fi

if grep -q "{ALGO}" "$SCRIPT"; then
    echo "ERROR: {ALGO} was not replaced in $SCRIPT"
    exit 1
fi

echo "Submitting $ALGO MIVIA: $((NJOBS+1)) jobs"
sbatch "$SCRIPT"