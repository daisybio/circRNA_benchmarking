from adjustText import adjust_text
from collections import defaultdict
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

from scipy.stats import pearsonr



# GLOBAL DICT TO STORE ADDITIONAL DATA FOR CORR PLOTS
MAIN_DICT = dict()
MAIN_DICT["LiHe"] = defaultdict(dict)
MAIN_DICT["BlCM"] = defaultdict(dict)
MAIN_DICT["BlEM"] = defaultdict(dict)
MAIN_DICT["BlTN"] = defaultdict(dict)
MAIN_DICT["mixed"] = defaultdict(dict)

MAIN_DICT["Liver"] = defaultdict(dict)
MAIN_DICT["TCell"] = defaultdict(dict)

# TODO currently have to manually create this dict -> for later need automation of course
reads_amount = {
"10_polyA":	65410853,
"11_polyA":	59622163,
"12_polyA":	60185671,
"13_polyA":	59765749,
"1_polyA":	47182096,
"2_polyA":	48097240,
"3_polyA":	46053189,
"4_polyA":	45655797,
"5_polyA":	58277442,
"6_polyA":	58796153,
"7_polyA":	58261358,
"8_polyA":	65080068,
"9_polyA":	65555480,
"10_total":	67067974,
"11_total":	66489369,
"12_total":	63801782,
"13_total":	63822771,
"1_total":	50795124,
"2_total":	51950222,
"3_total":	49704465,
"4_total":	49542984,
"5_total":	61553073,
"6_total":	59115762,
"7_total":	59093054,
"8_total":	70143244,
"9_total":	67273436
}

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
    tool = comp[1]

    # read bed
    df = pd.read_csv(
        bed_file_path,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
        ],
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


def compute_corr(bed_paths, featurecounts_paths, origin, tissue_type):
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
        print(f"{sample} {num_rRNA_sp_reads}")
        rRNA_dict[sample] = num_rRNA_sp_reads
        rrna_reads.append(str((num_rRNA_sp_reads / reads_amount[sample]) * 1000000))
        paths.append(featurecounts_file)
        featurecounts_samples.add(sample)

    if bsj_samples != featurecounts_samples:
        print("Samples are mismatching!")
        print(f"Featruecount samples: {featurecounts_samples}")
        print(f"BSJ_samples: {bsj_samples}")
        exit(1)

    with open(f"{origin}_{tissue_type}.rrnareads.txt", "w") as f:
        for i in range(len(paths)):
            f.write(f"{os.path.basename(paths[i])},{rrna_reads[i]}\n")

    # save bsj_dict
    with open(f"{origin}_{tissue_type}.bsj_amount.json", "w") as f:
        json.dump(bsj_dict, f, indent=4)
    # compute corr per tool
    correlations = dict()
    for tool in tools:
        X = []
        Y = []
        for sample in bsj_samples:
            X.append((bsj_dict[tool][sample], sample))
            rpm = (rRNA_dict[sample] / reads_amount[sample]) * 1000000
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


        MAIN_DICT[tissue_type][origin][tool] = {
            "bsj_scores": X_corr,
            "rrna_reads": Y_corr,
            "order": [o[1] for o in X_sorted],
        }
        print(tool)
        print(origin)
        print(tissue_type)
        print(X_corr)
        print(Y_corr)
        print(sorted_samples)

        print(np.corrcoef(X_corr, Y_corr)[0, 1])
        r, p = pearsonr(X_corr, Y_corr)
        print("r =", r, "p-value =", p)
        correlations[tool] = {'r': r, 'p' : p}
    return correlations


def rrna_analysis(total_bed_paths, total_featurecounts_paths, data_origin, tissue_type):
    json_path = f"{data_origin}_{tissue_type}.rna_rrna_corr.json"
    corr_dict = compute_corr(total_bed_paths, total_featurecounts_paths, data_origin, tissue_type)
    with open(json_path, "w") as f:
        json.dump(corr_dict, f, indent=4)


if __name__ == "__main__":
    tissue_types  = ["LiHe", "BlCM", "BlEM", "BlTN"]

    for tissue_type in tissue_types:
        for data_origin in ["total", "polya"]:
            bed_files = []
            featurecounts_files = []
            for di, _, files in os.walk(
                f"../../data_transformations_with_blacklist/featurecounts/{data_origin}"
            ):
                for file in files:
                    path_to_file = os.path.join(os.path.abspath(di), file)
                    if file.endswith("summary") and tissue_type in file:
                        featurecounts_files.append(path_to_file)

            for di, _, files in os.walk(f"../../data_with_blacklist/{data_origin}/filtered_bed_min"):
                for file in files:
                    path_to_file = os.path.join(os.path.abspath(di), file)
                    if path_to_file.endswith("bed") and tissue_type in path_to_file:
                        bed_files.append(path_to_file)

            rrna_analysis(bed_files, featurecounts_files, data_origin, tissue_type)
    
    tissue_groups  = [{"LiHe"}, {"BlCM", "BlEM", "BlTN"}]
    for tissue_group in tissue_groups:
        for data_origin in ["total", "polya"]:
            bed_files = []
            featurecounts_files = []
            for di, _, files in os.walk(
                f"../../data_transformations_with_blacklist/featurecounts/{data_origin}"
            ):
                for file in files:
                    path_to_file = os.path.join(os.path.abspath(di), file)
                    if file.endswith("summary") and  any(tissue in file for tissue in tissue_group):
                        featurecounts_files.append(path_to_file)

            for di, _, files in os.walk(f"../../data_with_blacklist/{data_origin}/filtered_bed_min"):
                for file in files:
                    path_to_file = os.path.join(os.path.abspath(di), file)
                    if path_to_file.endswith("bed") and any(tissue in file for tissue in tissue_group):
                        bed_files.append(path_to_file)

            if len(tissue_group) == 1:
                rrna_analysis(bed_files, featurecounts_files, data_origin, "Liver")
            else:
                rrna_analysis(bed_files, featurecounts_files, data_origin, "TCell")


    for data_origin in ["total", "polya"]:
        bed_files = []
        featurecounts_files = []
        for di, _, files in os.walk(
            f"../../data_transformations_with_blacklist/featurecounts/{data_origin}"
        ):
            for file in files:
                path_to_file = os.path.join(os.path.abspath(di), file)
                if file.endswith("summary"):
                    featurecounts_files.append(path_to_file)

        for di, _, files in os.walk(f"../../data_with_blacklist/{data_origin}/filtered_bed_min"):
            for file in files:
                path_to_file = os.path.join(os.path.abspath(di), file)
                if path_to_file.endswith("bed"):
                    bed_files.append(path_to_file)

        rrna_analysis(bed_files, featurecounts_files, data_origin, "mixed")


# Prepare the data
tissue_types.append('mixed')
for tissue_type in tissue_types:
    records = []
    for category, tools in MAIN_DICT[tissue_type].items():
        for tool, values in tools.items():
            for bsj, rrna, order in zip(
                values["bsj_scores"], values["rrna_reads"], values["order"]
            ):
                records.append(
                    {
                        "Tool": tool,
                        "BSJ Score": bsj,
                        "rRNA Reads": rrna,
                        "Category": category,
                        "Order": order,
                    }
                )


    df = pd.DataFrame(records)

    # Define consistent colors
    palette = {
                "total": "#FC8D62",  # orange
                "polya": "#66C2A5",  # green
    }


    # Unique tools
    tools = df["Tool"].unique()

    # Plot for each tool separately with regression line per category
    for tool in tools:
        tool_df = df[df["Tool"] == tool]
        plt.figure(figsize=(8, 6))

        sns.scatterplot(
            data=tool_df,
            x="BSJ Score",
            y="rRNA Reads",
            hue="Category",
            palette=palette,
            s=100,
            alpha=0.8,
        )

        for category in tool_df["Category"].unique():
            cat_df = tool_df[tool_df["Category"] == category]
            sns.regplot(
                data=cat_df,
                x="BSJ Score",
                y="rRNA Reads",
                scatter=False,
                ci=None,
                color=palette[category],
            )

        texts = []
        for _, row in tool_df.iterrows():
            texts.append(
                plt.text(
                    row["BSJ Score"],
                    row["rRNA Reads"],
                    row["Order"],
                    fontsize=9,
                )
            )

        # adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

        plt.title(f"{tool}\nBSJ Score vs rRNA Reads with Regression Line per Category ({tissue_type})")
        plt.xlabel(f"BSJ score sum")
        plt.grid(True)
        plt.legend()
        plt.ylim(
            0,
        )
        plt.xlim(
            0,
        )
        plt.tight_layout()
        plt.savefig(f"tool_plots/{tissue_type}_{tool}.png", dpi=300)
