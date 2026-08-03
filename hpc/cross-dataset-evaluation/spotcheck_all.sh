#!/bin/bash
# Spot-check: unified vs reference, all 7 algorithms, on a handful of BI pairs,
# both undirected and directed. Prints a clean PASS/MISMATCH/CRASH summary table.
#
# Usage: bash spotcheck_all.sh   (run from repo root, same place you run submit_bi.sh)

PAIRS=(
  "176.txt 192.txt"
  "176.txt 196.txt"
  "177.txt 192.txt"
  "177.txt 196.txt"
  "192.txt 196.txt"
)
ALGOS=("mcsplit" "rl" "dal" "ll" "dsb" "rrsplit" "symsplit")
INSTDIR="instances/SIP/biochemicalReactions"
LOGDIR="spotcheck_logs"
mkdir -p "$LOGDIR"

run_reference() {
  local algo=$1 dirflag=$2 A=$3 B=$4
  case "$algo" in
    mcsplit)
      REF_OUT=$(timeout 60 ./reference/mcsplit/mcsp $dirflag -l -q "$A" "$B" 2>&1)
      REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
      ;;
    rl)
      REF_OUT=$(timeout 60 ./reference/mcsplit-rl/mcsp+RL min_max $dirflag -l -q "$A" "$B" 2>&1)
      REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
      ;;
    dal)
      REF_OUT=$(timeout 60 ./reference/mcsplit-dal/mcspDAL min_max $dirflag -l -l -q "$A" "$B" 2>&1)
      REF_SIZE=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $2}')
      ;;
    ll)
      REF_OUT=$(timeout 60 ./reference/mcsplit-ll/mcsp+ll min_max $dirflag -l -q "$A" "$B" 2>&1)
      REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
      ;;
    dsb)
      REF_OUT=$(timeout 60 ./reference/mcsplit-dsb/bin/run.o min_max $dirflag -l -q "$A" "$B" 2>&1)
      REF_SIZE=$(echo "$REF_OUT" | awk -F',' '{print $1}' | tr -d ' ')
      ;;
    rrsplit)
      REF_OUT=$(timeout 60 ./reference/rrsplit/mcsp min_max $dirflag -l -q "$A" "$B" 2>&1)
      REF_SIZE=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $2}')
      ;;
    symsplit)
      REF_OUT=$(timeout 60 ./reference/symsplit/bin/run.o min_max $dirflag -l -q "$A" "$B" 2>&1)
      REF_SIZE=$(echo "$REF_OUT" | awk -F',' '{print $1}' | tr -d ' ')
      ;;
  esac
}

CSVFILE="spotcheck_results.csv"
echo "algo,pair,directed,unified_size,ref_size,status" > "$CSVFILE"

printf "%-10s %-18s %-10s %-10s %-10s %-10s\n" "ALGO" "PAIR" "DIRECTED" "UNIFIED" "REF" "STATUS"
printf '%.0s-' {1..75}; echo

for algo in "${ALGOS[@]}"; do
  case "$algo" in
    mcsplit) REF_DIRFLAG="-r" ;;
    *)       REF_DIRFLAG="-i" ;;
  esac

  for pair in "${PAIRS[@]}"; do
    A_NAME=$(echo $pair | cut -d' ' -f1)
    B_NAME=$(echo $pair | cut -d' ' -f2)
    A="$INSTDIR/$A_NAME"
    B="$INSTDIR/$B_NAME"
    PAIRTAG="${A_NAME%.txt}-${B_NAME%.txt}"

    for directed in "no" "yes"; do
      TAG="${algo}_${PAIRTAG}_${directed}"
      LOGFILE="$LOGDIR/${TAG}.log"

      if [ "$directed" = "yes" ]; then
        UNIFIED_FLAG="-i"
        DIRFLAG="$REF_DIRFLAG"
      else
        UNIFIED_FLAG=""
        DIRFLAG=""
      fi

      UNIFIED_OUT=$(timeout 60 ./solvers/bin/mcsp -A "$algo" -t 60 -l $UNIFIED_FLAG -q "$A" "$B" 2>&1)
      UNIFIED_SIZE=$(echo "$UNIFIED_OUT" | tail -1 | awk '{print $1}')

      run_reference "$algo" "$DIRFLAG" "$A" "$B"

      {
        echo "=== $TAG ==="
        echo "--- unified ---"
        echo "$UNIFIED_OUT"
        echo "--- reference ---"
        echo "$REF_OUT"
      } > "$LOGFILE"

      STATUS="MATCH"
      if echo "$REF_OUT" | grep -qi "mistaching\|Invalid solution"; then
        STATUS="REF_CRASH"
      elif [ -z "$UNIFIED_SIZE" ]; then
        STATUS="UNIFIED_EMPTY"
      elif [ -z "$REF_SIZE" ]; then
        STATUS="REF_EMPTY"
      elif [ "$UNIFIED_SIZE" != "$REF_SIZE" ]; then
        STATUS="MISMATCH"
      fi

      printf "%-10s %-18s %-10s %-10s %-10s %-10s\n" "$algo" "$PAIRTAG" "$directed" "$UNIFIED_SIZE" "$REF_SIZE" "$STATUS"
      echo "$algo,$PAIRTAG,$directed,$UNIFIED_SIZE,$REF_SIZE,$STATUS" >> "$CSVFILE"
    done
  done
done

echo
echo "Full logs (unified + reference raw output) in: $LOGDIR/"
echo "CSV summary: $CSVFILE"