import os
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor

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

            run_cmd = [ "featureCounts", "-p", "--countReadPairs", "-M", "--fraction", "-a", RRNA_GTF, "-o", sample_out, "-t", "exon", "-g", "gene_id", "-T", "4", path_to_sample, ]
            subprocess.run(run_cmd)

def count_rrna_reads(d_type: str):
    outdir_prefix = os.path.join(MAIN_DATA_DIR_, d_type, "rrna_reads")
    os.makedirs(outdir_prefix, exist_ok=True)
    
    tasks = []
    for dirpath, _, files in os.walk(os.path.join(MAIN_DATA_DIR, d_type, "bam")):
        for file in files:
            if not file.endswith(".bam"):
                continue
            path_to_sample = os.path.join(dirpath, file)
            tasks.append((file, path_to_sample))

    def count_rrna_reads(args):
        sample_name, path_to_sample = args
        run_cmd = [
            "bash", "-c",
            f"samtools view -b -F 260 {path_to_sample}"
            f" | bedtools intersect -u -abam - -b {RRNA_GTF}"
            f" | samtools view -c"
        ]
        result = subprocess.run(run_cmd, capture_output=True, text=True)
        count = result.stdout.strip()
        return sample_name, count

    tsv_out = os.path.join(outdir_prefix, "rrna_spanning_reads.tsv")
    
    with open(tsv_out, "a+") as f:
        f.seek(0)
        processed = set()
        for line in f:
            if line.startswith("File") or not line.strip():
                continue
            processed.add(line.split("\t")[0])
        
        for task in processed:
            print(f"Skipping sample {task} because already processed.")
         
        if len(processed) > 0:
            print(f"Delete {tsv_out} if you want to force recompute.")

        remaining_tasks = [t for t in tasks if t[0] not in processed]
        
        if not processed:
            f.write("File\tCount\n")
        
        if remaining_tasks:
            print(f"Counting rRNA reads in {len(remaining_tasks)} files")
            with ThreadPoolExecutor(max_workers=10) as executor:
                for sample_name, count in executor.map(count_rrna_reads, remaining_tasks):
                    f.write(f"{sample_name}\t{count}\n")


if __name__ == "__main__":
    #run_feature_counts("total")
    #run_feature_counts("polya")
    count_rrna_reads("total")
    count_rrna_reads("polya")