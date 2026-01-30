#!/usr/bin/sh

PREFIX_GSE138734="/nfs/data3/CIRCEST/runs/test_new_benchmarking/tidy_benchmarking/data/GSE138734/"
PREFIX_DEEP="/nfs/data3/CIRCEST/runs/test_new_benchmarking/tidy_benchmarking/data/DEEP/"
PLOT_OUT="/nfs/data3/CIRCEST/runs/test_new_benchmarking/tidy_benchmarking/results/overview"


# 1. Filter length 
# python3 ./preprocessing/filter_length.py $PREFIX_GSE138734
# python3 ./preprocessing/filter_min_threshold.py $PREFIX_GSE138734
# python3 ./preprocessing/filter_length.py $PREFIX_DEEP
# python3 ./preprocessing/filter_min_threshold.py $PREFIX_DEEP


# 2. Sort, Merge, Intersect
#python3 ./preprocessing/bedtools.py $PREFIX_GSE138734
#python3 ./preprocessing/bedtools.py $PREFIX_DEEP


python3 ./overview.py $PLOT_OUT

