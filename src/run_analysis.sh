#!/usr/bin/sh

PREFIX_GSE138734="/nfs/data3/CIRCEST/runs/test_new_benchmarking/tidy_benchmarking/data/GSE138734/"
PREFIX_DEEP="/nfs/data3/CIRCEST/runs/test_new_benchmarking/tidy_benchmarking/data/DEEP/"
PLOT_OUT="/nfs/data3/CIRCEST/runs/test_new_benchmarking/tidy_benchmarking/results/overview"
RRNA_GTF="/nfs/data3/CIRCEST/runs/test_new_benchmarking/tidy_benchmarking/data/rrna_gtf/ensembl_rRNA.gtf"


# 1. Filter length 
echo "Filtering BSJ based on length and evidence"
python3 ./preprocessing/filter_length.py $PREFIX_GSE138734
python3 ./preprocessing/filter_min_threshold.py $PREFIX_GSE138734
python3 ./preprocessing/filter_length.py $PREFIX_DEEP
python3 ./preprocessing/filter_min_threshold.py $PREFIX_DEEP


# 2. Sort, Merge, Intersect
echo "Merging and intersecting BEDs"
python3 ./preprocessing/bedtools.py $PREFIX_GSE138734
python3 ./preprocessing/bedtools.py $PREFIX_DEEP

# 3. BED Analysis
echo "Running BED analysis"
python3 ./overview.py $PLOT_OUT


# 4. RRNA Corr Analysis
PREFIX_DEEP_="/nfs/data3/CIRCEST/runs/test_new_benchmarking/benchmarking_analysis/data_with_blacklist/"
PREFIX_GSE138734_="/nfs/data3/CIRCEST/runs/test_new_benchmarking/benchmarking_analysis/data_NIH/"

count_reads () {
    local bam_dir="$1"
    local meta_dir="$2"
    local outfile="${meta_dir}/counts.txt"

    mkdir -p "$meta_dir"

    if [ ! -s "$outfile" ]; then
        echo "Counting mapped primary reads in ${bam_dir}"
        bash ./util/count_entries_in_bam.sh "$bam_dir" > "$outfile"
    else
        echo "Skipping ${bam_dir} (meta $outfile already exists and is non-empty)"
    fi
}

# get number of primary mapped reads only
echo "Counting primary mapped reads..."
count_reads "${PREFIX_DEEP_%/}/total/bam"     "${PREFIX_DEEP%/}/total/bam_meta"
count_reads "${PREFIX_DEEP_%/}/polya/bam"     "${PREFIX_DEEP%/}/polya/bam_meta"
count_reads "${PREFIX_GSE138734_%/}/total/bam" "${PREFIX_GSE138734%/}/total/bam_meta"
count_reads "${PREFIX_GSE138734_%/}/polya/bam" "${PREFIX_GSE138734%/}/polya/bam_meta"

echo "Counting rrna spanning reads..."
python3 ./util/count_rrna_reads.py $PREFIX_DEEP_ $RRNA_GTF $PREFIX_DEEP
python3 ./util/count_rrna_reads.py $PREFIX_GSE138734_ $RRNA_GTF $PREFIX_GSE138734

echo "Running rRNA analysis"
python3 ./rrna_analysis/rrna_corr.py $PREFIX_GSE138734  ../results/rrna_analysis/GSE138734
python3 ./rrna_analysis/rrna_corr.py $PREFIX_DEEP ../results/rrna_analysis/DEEP


echo "Plotting correlation dumbbell plot"
python3 ./rrna_analysis/plot_rrna_corr.py ../results/rrna_analysis/DEEP
python3 ./rrna_analysis/plot_rrna_corr.py ../results/rrna_analysis/GSE138734


echo "Plotting BSJ sum distribution"
python3 ./rrna_analysis/misc/plot_bsj_sum.py ../results/rrna_analysis/GSE138734
python3 ./rrna_analysis/misc/plot_bsj_sum.py ../results/rrna_analysis/DEEP
 
 
echo "Plotting rRNA %"
python3 ./rrna_analysis/misc/plot_norm_rrna_reads.py $PREFIX_DEEP ../results/rrna_analysis/DEEP 
python3 ./rrna_analysis/misc/plot_norm_rrna_reads.py $PREFIX_GSE138734 ../results/rrna_analysis/GSE138734


echo "DONE"