#!/bin/bash
# Superseded by hpc/resubmit_missing.sh (takes the dataset as an argument). Kept so old commands still work.
# Usage: bash hpc/resubmit_missing_bi.sh <algo> [options]   ==   bash hpc/resubmit_missing.sh <algo> bi [options]
exec bash hpc/resubmit_missing.sh "$1" bi "${@:2}"
