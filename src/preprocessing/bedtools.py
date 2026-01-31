import os
import subprocess
import sys

MAIN_DATA_DIR = sys.argv[1]
BEDTYPE = "filtered_bed_min_blacklist"


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
    run_sort_concat("total")
    run_sort_concat("polya")
    print("Sorted and concatenated beds")


    run_merge("total")
    run_merge("polya")
    print("Merged concatenated beds")

    ## extract tools
    tools = set()
    for dirpath, _, files in os.walk(os.path.join(MAIN_DATA_DIR, "polya", BEDTYPE)):
        if dirpath.endswith("filtered_blacklist"):
            in_path = os.path.abspath(dirpath)
            tool = dirpath.split("/")[-1]
            tools.add(tool)

    run_intersect(tools)
    print("Intersected merged beds")

    run_sort_concat_all("total")
    run_sort_concat_all("polya")
    print("Merged tool-wise merged beds")
    
    run_merge_all("total")
    run_merge_all("polya")
    print("Intersected merged tool-wise merged beds")
