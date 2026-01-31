from adjustText import adjust_text
from collections import defaultdict
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

from scipy.stats import pearsonr
import sys



MAIN_DATA_DIR = sys.argv[1]
OUT_DIR = sys.argv[2]

# parse read amount per sample bam
READS_AMOUNT = dict()

def parse_bed(bed_file_path):
    """
    Reads bed file using pyranges and returns a two keys and the sum
    of bsj reads from the bed.
    :returns: sample, tool, sum
    """
    # get file name and extract tool and sample
    filename = os.path.basename(bed_file_path)
    # filename needs to be formatted like this: sample.tool.bed
    comp = filename.split(".")
    sample = "_".join(comp[0].split("_")[0:2])

    tool = bed_file_path.split("/")[-2].replace("_filtered_blacklist", "")

    # read bed
    df = pd.read_csv(
        bed_file_path,
        sep="\t",
        header=None,
        usecols=range(6),
        names=[
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
        ],
        index_col=False,
    )

    # sum # of bsj_reads
    num_bsj_reads = df["score"].sum()

    return sample, tool, int(num_bsj_reads)


def parse_featurecount(featurecounts_file_path):
    """
    Reads featurecounts file and returns a key and the sum of rRNA spanning reads per tool.
    :returns: tool, sum
    """
    # get file name and extract tool and sample
    filename = os.path.basename(featurecounts_file_path)
    # filename needs to be formatted like this: sample.origin.featurecounts.txt.summary
    # (Important is that the first element after split(".") is the sample name)
    comp = filename.split(".")
    # sample = comp[0]
    sample = "_".join(comp[0].split("_")[0:2])


    df = pd.read_csv(featurecounts_file_path, sep="\t", header=None, index_col=0)
    df = df.transpose()
    df.iloc[:, 1:] = df.iloc[:, 1:].astype(int)
    # num_rRNA_sp_reads = df["Assigned"] + df["Unassigned_MultiMapping"]
    num_rRNA_sp_reads = df["Assigned"] # now I am calling feature counts correctly (with -M and --fraction)
    return sample, num_rRNA_sp_reads.iloc[0]


def compute_corr(bed_paths, featurecounts_paths, origin):
    bsj_dict = defaultdict(dict)
    rRNA_dict = dict()

    bsj_samples = set()
    featurecounts_samples = set()
    tools = set()

    # load num_bsj_reads into dict:
    # tools -> samples -> nums
    for bed_file_path in bed_paths:
        sample, tool, num_bsj_reads = parse_bed(bed_file_path.strip("\n"))
        bsj_dict[tool][sample] = num_bsj_reads
        bsj_samples.add(sample)
        tools.add(tool)

    # extract number of rRNA spanning reads per sample
    rrna_reads = []
    paths = []
    for featurecounts_file in featurecounts_paths:
        sample, num_rRNA_sp_reads = parse_featurecount(featurecounts_file.strip("\n"))
        rRNA_dict[sample] = num_rRNA_sp_reads
        rrna_reads.append(str((num_rRNA_sp_reads / READS_AMOUNT[sample]) * 1000000))
        paths.append(featurecounts_file)
        featurecounts_samples.add(sample)

    if bsj_samples != featurecounts_samples:
        print("Samples are mismatching!")
        print(f"Featruecount samples: {featurecounts_samples}")
        print(f"BSJ_samples: {bsj_samples}")
        exit(1)

    rrna_cmp_out = os.path.join(OUT_DIR,f"{origin}.rrnareads_cpm.txt")
    with open(rrna_cmp_out, "w") as f:
        for i in range(len(paths)):
            f.write(f"{os.path.basename(paths[i])},{rrna_reads[i]}\n")
    
    # save bsj_dict
    bsj_out = os.path.join(OUT_DIR,f"{origin}.bsj_amount.json")
    with open(bsj_out, "w") as f:
        json.dump(bsj_dict, f, indent=4)


    # compute corr per tool
    correlations = dict()
    for tool in tools:
        X = []
        Y = []
        for sample in bsj_samples:
            X.append((bsj_dict[tool][sample], sample))
            rpm = (rRNA_dict[sample] / READS_AMOUNT[sample]) * 1000000
            Y.append((rpm, sample))

        # here we sort tupes of X (num_bsj_reads, sample) by num_bsj_reads and
        # apply same (sample based) order to Y
        # then we can check for linear relation ship
        X_sorted = sorted(X, key=lambda x: x[0])
        sorted_samples = [sample for _, sample in X_sorted]
        Y_sorted = [next(y for y in Y if y[1] == sample) for sample in sorted_samples]

        # last sanity check
        X_corr = []
        Y_corr = []
        for i in range(len(X_sorted)):
            if X_sorted[i][1] != Y_sorted[i][1]:
                print("Order of samples incorrect. Can't calculate correlation")
                print(f"X: {X_sorted}")
                print(f"Y: {Y_sorted}")
                exit(1)
            X_corr.append(X_sorted[i][0])
            Y_corr.append(Y_sorted[i][0])


        # compute correlation and significance
        r, p_value = pearsonr(X_corr, Y_corr)
        correlations[tool] = {"r": r, "p": p_value}

    return correlations


def rrna_analysis(total_bed_paths, total_featurecounts_paths, data_origin):
    json_path = f"{data_origin}.rna_rrna_corr.json"
    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, json_path)
    corr_dict = compute_corr(total_bed_paths, total_featurecounts_paths, data_origin)
    with open(out, "w") as f:
        json.dump(corr_dict, f, indent=4)


def read_bam_meta(path_to_meta):
    with open(path_to_meta) as f:
        lines = f.readlines()
        for line in lines:
            if "File" in line:
                continue
            sample = "_".join(line.split("\t")[0].split("/")[-1].split(".")[0].split("_")[0:2])
            amount = int(line.split("\t")[-1].strip("\n"))
            READS_AMOUNT[sample] = amount


if __name__ == "__main__":
    polya_bam_amount = os.path.join(MAIN_DATA_DIR, "polya", "bam_meta", "counts.txt")
    total_bam_amount = os.path.join(MAIN_DATA_DIR, "total", "bam_meta", "counts.txt")
    read_bam_meta(polya_bam_amount)
    read_bam_meta(total_bam_amount)


    for data_origin in ["total", "polya"]:
        bed_files = []
        featurecounts_files = []
        featurecounts_dir = os.path.join(MAIN_DATA_DIR, data_origin, "rrna_reads")
        for di, _, files in os.walk(featurecounts_dir):
            for file in files:
                if "summary" in file:
                    path_to_file = os.path.join(os.path.abspath(di), file)
                    featurecounts_files.append(path_to_file)

        bed_dir = os.path.join(MAIN_DATA_DIR, data_origin, "filtered_bed_min_blacklist")
        for di, _, files in os.walk(bed_dir):
            for file in files:
                path_to_file = os.path.join(os.path.abspath(di), file)
                bed_files.append(path_to_file)

        rrna_analysis(bed_files, featurecounts_files, data_origin)