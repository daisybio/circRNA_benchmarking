import os
import subprocess
import sys

MAIN_DATA_DIR = sys.argv[1]
# RRNA_GTF = sys.argv[2]
BEDTYPE = "filtered_bed_min_blacklist"


def run_feature_counts(d_type: str):
    outdir_prefix = os.path.join("featurecounts", d_type)

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
            print(" ".join(run_cmd))
            subprocess.run(run_cmd)


def run_merge(d_type: str):
    outdir_prefix = os.path.join(MAIN_DATA_DIR, "merge_concatenated_beds", d_type)
    os.makedirs(outdir_prefix, exist_ok=True)

    for dirpath, _, files in os.walk(
        os.path.join(MAIN_DATA_DIR, "./concatenated_bed_per_tool/", d_type)
    ):
        for file in files:
            in_path = os.path.join(os.path.abspath(dirpath), file)
            tool = file.split(".")[0]
            out_path = f"{outdir_prefix}/{tool}.{d_type}.merged.bed"
            cmd = f"bedtools merge -i {in_path} -c 4,5,6 -o distinct -s > {out_path}"
            # print(cmd)
            subprocess.run(cmd, shell=True, check=True)


def run_sort_concat(d_type: str):
    outdir_prefix = os.path.join(MAIN_DATA_DIR, "concatenated_bed_per_tool", d_type)
    os.makedirs(outdir_prefix, exist_ok=True)

    for dirpath, _, files in os.walk(os.path.join(MAIN_DATA_DIR, d_type, BEDTYPE)):
        if dirpath.endswith("filtered_blacklist"):
            in_path = os.path.abspath(dirpath)
            tool = dirpath.split("/")[-1]
            out_path = f"{outdir_prefix}/{tool}.{d_type}.all.bed"
            cmd = f"sort -k1,1 -k2,2n {in_path}/* > {out_path}"
            subprocess.run(cmd, shell=True, check=True)


def run_sort_concat_all(d_type: str):
    outdir_prefix = os.path.join(MAIN_DATA_DIR, "concatenated_all", d_type)
    os.makedirs(outdir_prefix, exist_ok=True)

    in_path = os.path.join(MAIN_DATA_DIR, "merge_concatenated_beds", d_type)
    out_path = f"{outdir_prefix}/{d_type}_concat_all.bed"
    cmd = f"sort -k1,1 -k2,2n {in_path}/* > {out_path}"
    subprocess.run(cmd, shell=True, check=True)

def run_merge_all(d_type: str):
    outdir_prefix = os.path.join(MAIN_DATA_DIR, "merged_all", d_type)
    os.makedirs(outdir_prefix, exist_ok=True)
    in_path = os.path.join(MAIN_DATA_DIR, "concatenated_all/", d_type, f"{d_type}_concat_all.bed")
    out_path = f"{outdir_prefix}/{d_type}_merged_all.bed"
    cmd = f"bedtools merge -i {in_path} -c 4,5,6 -o distinct -s > {out_path}"
    subprocess.run(cmd, shell=True, check=True)


def run_intersect(tools):
    os.makedirs(os.path.join(MAIN_DATA_DIR, "intersect_bed_pairs"), exist_ok=True)
    for tool in tools:
        out_path = os.path.join(MAIN_DATA_DIR, "intersect_bed_pairs", f"{tool}.intersect.bed")
        polya_path = os.path.join(MAIN_DATA_DIR, f"merge_concatenated_beds/polya/{tool}.polya.merged.bed")
        total_path = os.path.join(MAIN_DATA_DIR, f"merge_concatenated_beds/total/{tool}.total.merged.bed")

        cmd = f"bedtools intersect -a {polya_path} -b {total_path} > {out_path}"
        subprocess.run(cmd, shell=True, check=True)


if __name__ == "__main__":
    #run_feature_counts("total")
    #run_feature_counts("polya")
    print("Counted rRNA spanning reads")
    
    #run_sort_concat("total")
    #run_sort_concat("polya")
    #print("Sorted and concatenated beds")


    #run_merge("total")
    #run_merge("polya")
    #print("Merged concatenated beds")

    ### extract tools
    #tools = set()
    #for dirpath, _, files in os.walk(os.path.join(MAIN_DATA_DIR, "polya", BEDTYPE)):
    #    if dirpath.endswith("filtered_blacklist"):
    #        in_path = os.path.abspath(dirpath)
    #        tool = dirpath.split("/")[-1]
    #        tools.add(tool)

    #run_intersect(tools)
    #print("Intersected merged beds")

    run_sort_concat_all("total")
    run_sort_concat_all("polya")
    
    run_merge_all("total")
    run_merge_all("polya")
