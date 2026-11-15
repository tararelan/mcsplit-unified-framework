#!/bin/bash
# resubmit_missing.sh - rerun only the pairs that have no usable result for one algorithm on one dataset
# Usage: bash hpc/resubmit_missing.sh <algo> <mivia|lv|bi> [--with-ref] [--swap] [--mem=40G] [--time=02:00:00] [--dry-run]
#   First refresh the pair lists with:  python hpc/find_missing.py
#   Use the SAME --with-ref / --swap options as the original run, so the results land in the same folder.
#   The pair list is hpc/missing/<algo>_<dataset>.txt (pairs with no unified result in hpc/publication).
#   --dry-run prints what would happen and changes nothing (no files moved, nothing submitted).
#
# run_instance.sh skips any pair whose result file already exists, and a killed run leaves an empty result
# file behind. Those files are moved to <results folder>/_empty_backup/ first so the pair is really rerun;
# nothing is deleted.

ALGO=$1
DS=$2
WITH_REF=0; SWAP=0; MEM=""; TIME=""; DRY=0
for opt in "${@:3}"; do
    case "$opt" in
        --with-ref) WITH_REF=1 ;;
        --swap)     SWAP=1 ;;
        --mem=*)    MEM="${opt#--mem=}" ;;
        --time=*)   TIME="${opt#--time=}" ;;
        --dry-run)  DRY=1 ;;
        *) echo "Unknown option: $opt"; exit 1 ;;
    esac
done

LIST="hpc/missing/${ALGO}_${DS}.txt"
case "$DS" in mivia|lv|bi) ;; *) DS_BAD=1 ;; esac
if [ -z "$ALGO" ] || [ -n "$DS_BAD" ] || [ ! -f "$LIST" ]; then
    echo "Usage: bash hpc/resubmit_missing.sh <algo> <mivia|lv|bi> [--with-ref] [--swap] [--mem=40G] [--time=02:00:00] [--dry-run]"
    echo "Pair lists available (run 'python hpc/find_missing.py' to refresh): $(ls hpc/missing/*.txt 2>/dev/null | xargs -n1 basename 2>/dev/null | sed 's/\.txt$//' | tr '\n' ' ')"
    exit 1
fi

# same result folder as run_instance.sh
if [ "$WITH_REF" = "1" ]; then OUTDIR="hpc/results/${ALGO}_${DS}_with_ref"; else OUTDIR="hpc/results/${ALGO}_${DS}"; fi
[ "$SWAP" = "1" ] && OUTDIR="${OUTDIR}_swap"

moved=0
[ "$DRY" = "1" ] || mkdir -p hpc/logs "$OUTDIR/_empty_backup"
while read -r A B; do
    f="${OUTDIR}/$(basename "$A").$(basename "$B").csv"
    if [ -f "$f" ]; then
        [ "$DRY" = "1" ] || mv -n "$f" "$OUTDIR/_empty_backup/"
        moved=$((moved+1))
    fi
done < "$LIST"
echo "$ALGO $DS: $(wc -l < "$LIST") pairs to rerun, $moved existing result files $( [ "$DRY" = "1" ] && echo 'would be moved' || echo 'moved') to $OUTDIR/_empty_backup/"

NJOBS=$(($(wc -l < "$LIST") - 1))
SCRIPT="hpc/job_${ALGO}_${DS}_missing.slurm"
sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g; s/{DATASET}/$DS/g; s/{WITH_REF}/$WITH_REF/g; s#{LIST}#$LIST#g" \
    hpc/job_array.slurm | sed "s/{SWAP}/$SWAP/g" > "$SCRIPT"
sed -i 's/\r$//' "$SCRIPT"
[ -n "$MEM" ]  && sed -i "s/^#SBATCH --mem=.*/#SBATCH --mem=$MEM/" "$SCRIPT"
[ -n "$TIME" ] && sed -i "s/^#SBATCH --time=.*/#SBATCH --time=$TIME/" "$SCRIPT"

if grep -qE "\{ALGO\}|\{NJOBS\}|\{DATASET\}|\{LIST\}|\{WITH_REF\}|\{SWAP\}" "$SCRIPT"; then
    echo "ERROR: placeholder not replaced in $SCRIPT"; exit 1
fi
if [ "$DRY" = "1" ]; then echo "dry run: would sbatch $SCRIPT"; grep '^#SBATCH' "$SCRIPT"; rm -f "$SCRIPT" 2>/dev/null; exit 0; fi
sbatch "$SCRIPT"
