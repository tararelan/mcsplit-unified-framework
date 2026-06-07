#!/bin/bash
ALGOS="mcsplit rl rllum dal ll dsb rrsplit symsplit"

for ALGO in $ALGOS; do
    OUTFILE="results/hpc/${ALGO}_all.csv"
    echo "instance_a,instance_b,solution_size,solution_edges,nodes,time_elapsed,aborted,root_ub,nodes_to_best,time_to_best,cut_branches,bound_pruned,sym_pruned,conflicts" > $OUTFILE
    for f in results/hpc/${ALGO}/*.csv; do
        cat $f >> $OUTFILE
    done
    echo "Collected $ALGO: $(wc -l < $OUTFILE) instances"
done