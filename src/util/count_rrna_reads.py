import os
import subprocess
import sys

MAIN_DATA_DIR = sys.argv[1]
RRNA_GTF = sys.argv[2] 
MAIN_DATA_DIR_ = sys.argv[3]

def run_feature_counts(d_type: str):
    outdir_prefix = os.path.join(MAIN_DATA_DIR_, d_type, "rrna_reads")
    os.makedirs(outdir_prefix, exist_ok=True)

    for dirpath, e, files in os.walk(os.path.join(MAIN_DATA_DIR, d_type, "bam")):
        for file in files:
            sample_name = os.path.basename(file)
            sample_out = os.path.join(
                outdir_prefix, f"{sample_name}.{d_type}.featurecounts.txt"
            )
            path_to_sample = os.path.join(dirpath, file)

            run_cmd = [
                "featureCounts",
                "-p",
                "--countReadPairs",
                "-M",
                "--fraction",
                "-a",
                RRNA_GTF,
                "-o",
                sample_out,
                "-t",
                "exon",
                "-g",
                "gene_id",
                "-T",
                "4",
                path_to_sample,
            ]
            subprocess.run(run_cmd)


if __name__ == "__main__":
    run_feature_counts("total")
    run_feature_counts("polya")