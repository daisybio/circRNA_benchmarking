import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from upsetplot import from_memberships, UpSet
import plotly.graph_objects as go
import numpy as np

import os
import subprocess
import shlex
import sys

PLOT_OUT = sys.argv[1]
os.makedirs(PLOT_OUT, exist_ok=True)


 
sns.set_theme(style="ticks", context="paper", palette="colorblind")

tool_order_a = ["total", "polyA"]

sns.set_context(
    "paper",
    rc={
        "font.size": 15,  # base font size
        "axes.titlesize": 17,  # title
        "axes.labelsize": 15,  # axis labels
        "xtick.labelsize": 14,  # x tick labels
        "ytick.labelsize": 14,  # y tick labels
        "legend.fontsize": 12,  # legend
        "lines.linewidth": 1.5,  # line width
        "lines.markersize": 6,  # marker size
    },
)

plt.rcParams["figure.figsize"] = (8, 4.5)  # width x height in inches

cb_palette = sns.color_palette("colorblind")

# update your palette
palette = {
    "total": cb_palette[0],   # 1st color from seaborn's colorblind palette
    "polya": cb_palette[1],       # keep your green
}

def run_jaccard_matrix(tools):
    # Initialize empty DataFrames (rows=tool1, cols=tool2)
    nih_df = pd.DataFrame(index=tools, columns=tools, dtype=float)
    deep_df = pd.DataFrame(index=tools, columns=tools, dtype=float)
    
    nih_df_poly = pd.DataFrame(index=tools, columns=tools, dtype=float)
    deep_df_poly = pd.DataFrame(index=tools, columns=tools, dtype=float)
    
    nih_df_total = pd.DataFrame(index=tools, columns=tools, dtype=float)
    deep_df_total = pd.DataFrame(index=tools, columns=tools, dtype=float)

    # polya vs total
    for tool1 in tools:
        for tool2 in tools:
            # NIH run
            cmd = f'bedtools jaccard -a ../data/GSE138734/merge_concatenated_beds/polya/{tool1}_filtered_blacklist.polya.merged.bed -b ../data/GSE138734/merge_concatenated_beds/total/{tool2}_filtered_blacklist.total.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            nih_res = proc.stdout.strip().split("\n")[-1].split("\t")
            jaccard_val = float(nih_res[2])  # column 3 = jaccard
            nih_df.loc[tool1, tool2] = jaccard_val

            # DEEP run
            cmd = f'bedtools jaccard -a ../data/DEEP/merge_concatenated_beds/polya/{tool1}_filtered_blacklist.polya.merged.bed -b ../data/DEEP/merge_concatenated_beds/total/{tool2}_filtered_blacklist.total.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            deep_res = proc.stdout.strip().split("\n")[-1].split("\t")
            jaccard_val = float(deep_res[2])
            deep_df.loc[tool1, tool2] = jaccard_val
    
    # polya vs polya
    for tool1 in tools:
        for tool2 in tools:
            # NIH run
            cmd = f'bedtools jaccard -a ../data/GSE138734/merge_concatenated_beds/polya/{tool1}_filtered_blacklist.polya.merged.bed -b ../data/GSE138734/merge_concatenated_beds/polya/{tool2}_filtered_blacklist.polya.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            nih_res = proc.stdout.strip().split("\n")[-1].split("\t")
            jaccard_val = float(nih_res[2])  # column 3 = jaccard
            nih_df_poly.loc[tool1, tool2] = jaccard_val

            # DEEP run
            cmd = f'bedtools jaccard -a ../data/DEEP/merge_concatenated_beds/polya/{tool1}_filtered_blacklist.polya.merged.bed -b ../data/DEEP/merge_concatenated_beds/polya/{tool2}_filtered_blacklist.polya.merged.bed'

            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            deep_res = proc.stdout.strip().split("\n")[-1].split("\t")
            jaccard_val = float(deep_res[2])
            deep_df_poly.loc[tool1, tool2] = jaccard_val
    
    # total vs total
    for tool1 in tools:
        for tool2 in tools:
            # NIH run
            cmd = f'bedtools jaccard -a ../data/GSE138734/merge_concatenated_beds/total/{tool1}_filtered_blacklist.total.merged.bed -b ../data/GSE138734/merge_concatenated_beds/total/{tool2}_filtered_blacklist.total.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            nih_res = proc.stdout.strip().split("\n")[-1].split("\t")
            jaccard_val = float(nih_res[2])  # column 3 = jaccard
            nih_df_total.loc[tool1, tool2] = jaccard_val

            # DEEP run
            cmd = f'bedtools jaccard -a ../data/DEEP/merge_concatenated_beds/total/{tool1}_filtered_blacklist.total.merged.bed -b ../data/DEEP/merge_concatenated_beds/total/{tool2}_filtered_blacklist.total.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            deep_res = proc.stdout.strip().split("\n")[-1].split("\t")
            jaccard_val = float(deep_res[2])
            deep_df_total.loc[tool1, tool2] = jaccard_val

    
    name_map = { 'segemehl': 'Segemehl',
    'dcc': 'DCC',
    'ciriquant': 'CIRIquant',
    'circexplorer2': 'CIRCexplorer2',
    'find_circ':'find_circ'
    }
    # Apply the mapping to both index and columns
    nih_df = nih_df.rename(index=name_map, columns=name_map)
    deep_df = deep_df.rename(index=name_map, columns=name_map)
    nih_df_total = nih_df_total.rename(index=name_map, columns=name_map)
    deep_df_total = deep_df_total.rename(index=name_map, columns=name_map)
    nih_df_poly = nih_df_poly.rename(index=name_map, columns=name_map)
    deep_df_poly = deep_df_poly.rename(index=name_map, columns=name_map)

    return nih_df, deep_df, nih_df_poly, deep_df_poly, nih_df_total, deep_df_total


def run_compare_matrix(tools):
    # Initialize empty DataFrames (rows=tool1, cols=tool2)
    nih_df = pd.DataFrame(index=tools, columns=tools, dtype=float)
    deep_df = pd.DataFrame(index=tools, columns=tools, dtype=float)
    
    nih_df_poly = pd.DataFrame(index=tools, columns=tools, dtype=float)
    deep_df_poly = pd.DataFrame(index=tools, columns=tools, dtype=float)
    
    nih_df_total = pd.DataFrame(index=tools, columns=tools, dtype=float)
    deep_df_total = pd.DataFrame(index=tools, columns=tools, dtype=float)

    # polya vs total
    for tool1 in tools:
        for tool2 in tools:
            # NIH run
            cmd = f'bash  ./util/compare.sh ../data/GSE138734/merge_concatenated_beds/polya/{tool1}_filtered_blacklist.polya.merged.bed ../data/GSE138734/merge_concatenated_beds/total/{tool2}_filtered_blacklist.total.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            nih_res = proc.stdout.strip()
            print(nih_res)
            jaccard_val = float(nih_res)  # column 3 = jaccard
            if jaccard_val > 1: # hard cap at 1 since some times merging A and B results in smaller file then intersection A B
                jaccard_val = 1
            nih_df.loc[tool1, tool2] = jaccard_val

            # DEEP run
            cmd = f'bash  ./util/compare.sh ../data/DEEP/merge_concatenated_beds/polya/{tool1}_filtered_blacklist.polya.merged.bed ../data/DEEP/merge_concatenated_beds/total/{tool2}_filtered_blacklist.total.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            deep_res = proc.stdout.strip()
            jaccard_val = float(deep_res)
            deep_df.loc[tool1, tool2] = jaccard_val
    
    # polya vs polya
    for tool1 in tools:
        for tool2 in tools:
            # NIH run
            cmd = f'bash  ./util/compare.sh ../data/GSE138734/merge_concatenated_beds/polya/{tool1}_filtered_blacklist.polya.merged.bed ../data/GSE138734/merge_concatenated_beds/polya/{tool2}_filtered_blacklist.polya.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            nih_res = proc.stdout.strip()
            jaccard_val = float(nih_res)  # column 3 = jaccard
            if jaccard_val > 1: # hard cap at 1 since some times merging A and B results in smaller file then intersection A B
                jaccard_val = 1
            nih_df_poly.loc[tool1, tool2] = jaccard_val

            # DEEP run
            cmd = f'bash  ./util/compare.sh ../data/DEEP/merge_concatenated_beds/polya/{tool1}_filtered_blacklist.polya.merged.bed ../data/DEEP/merge_concatenated_beds/polya/{tool2}_filtered_blacklist.polya.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            deep_res = proc.stdout.strip()
            jaccard_val = float(deep_res)
            if jaccard_val > 1: # hard cap at 1 since some times merging A and B results in smaller file then intersection A B
                jaccard_val = 1
            deep_df_poly.loc[tool1, tool2] = jaccard_val
    
    # total vs total
    for tool1 in tools:
        for tool2 in tools:
            # NIH run
            cmd = f'bash  ./util/compare.sh ../data/GSE138734/merge_concatenated_beds/total/{tool1}_filtered_blacklist.total.merged.bed ../data/GSE138734/merge_concatenated_beds/total/{tool2}_filtered_blacklist.total.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            nih_res = proc.stdout.strip()
            jaccard_val = float(nih_res)  # column 3 = jaccard
            if jaccard_val > 1: # hard cap at 1 since some times merging A and B results in smaller file then intersection A B
                jaccard_val = 1
            nih_df_total.loc[tool1, tool2] = jaccard_val

            # DEEP run
            cmd = f'bash  ./util/compare.sh ../data/DEEP/merge_concatenated_beds/total/{tool1}_filtered_blacklist.total.merged.bed  ../data/DEEP/merge_concatenated_beds/total/{tool2}_filtered_blacklist.total.merged.bed'
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            deep_res = proc.stdout.strip()
            jaccard_val = float(deep_res)
            if jaccard_val > 1: # hard cap at 1 since some times merging A and B results in smaller file then intersection A B
                jaccard_val = 1
            deep_df_total.loc[tool1, tool2] = jaccard_val

    name_map = { 'segemehl': 'Segemehl',
    'dcc': 'DCC',
    'ciriquant': 'CIRIquant',
    'circexplorer2': 'CIRCexplorer2',
    'find_circ':'find_circ'
    }
    
    # Apply the mapping to both index and columns
    nih_df = nih_df.rename(index=name_map, columns=name_map)
    deep_df = deep_df.rename(index=name_map, columns=name_map)
    nih_df_total = nih_df_total.rename(index=name_map, columns=name_map)
    deep_df_total = deep_df_total.rename(index=name_map, columns=name_map)
    nih_df_poly = nih_df_poly.rename(index=name_map, columns=name_map)
    deep_df_poly = deep_df_poly.rename(index=name_map, columns=name_map)

    return nih_df, deep_df, nih_df_poly, deep_df_poly, nih_df_total, deep_df_total

def plot_jaccard_heatmaps(nih_df, deep_df, name):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # NIH heatmap
    sns.heatmap(
        nih_df.astype(float), 
        ax=axes[0], 
        cmap="viridis", 
        annot=True, 
        fmt=".3f", 
        cbar=True,
        vmin=0,
        vmax=1

    )
    axes[0].set_title(f"Jaccard Index (GSE138734, {name})")

    # DEEP heatmap
    sns.heatmap(
        deep_df.astype(float), 
        ax=axes[1], 
        cmap="viridis", 
        annot=True, 
        fmt=".3f", 
        cbar=True,
        vmin=0,
        vmax=1

    )
    axes[1].set_title(f"Jaccard Index (DEEP, {name})")

    name = name.replace(" ", "_")
    if name == "polyA_vs_Total":
        axes[0].set_ylabel("Poly(A) BED")
        axes[0].set_xlabel("Total BED")
        axes[1].set_ylabel("Poly(A) BED")
        axes[1].set_xlabel("Total BED")
    elif name == "polyA":
        axes[0].set_ylabel("Poly(A) BED")
        axes[0].set_xlabel("Poly(A) BED")
        axes[1].set_ylabel("Poly(A) BED")
        axes[1].set_xlabel("Poly(A) BED")
    else:
        axes[0].set_ylabel("Total BED")
        axes[0].set_xlabel("Total BED")
        axes[1].set_ylabel("Total BED")
        axes[1].set_xlabel("Total BED")

    axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
    axes[1].set_yticklabels(axes[1].get_yticklabels(), rotation=0)

    plt.tight_layout()
    out = os.path.join(PLOT_OUT, f"{name}_jaccard_heat.png")
    plt.savefig(out, dpi = 300)

def plot_inter_heatmaps(nih_df, deep_df, name):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Row-wise normalization for colors
    def normalize_rows(df):
        return df.div(df.max(axis=1), axis=0)

    # Normalized for coloring
    nih_colors = normalize_rows(nih_df.astype(float))
    deep_colors = normalize_rows(deep_df.astype(float))

    # NIH heatmap
    sns.heatmap(
        nih_colors, 
        ax=axes[0], 
        cmap="viridis", 
        annot=nih_df.astype(float),  # Keep absolute values
        fmt=".0f", 
        cbar=True,
        vmin=0,
        vmax=1
    )
    axes[0].set_title(f"Intersecting BSJs (GSE138734, {name})")

    # DEEP heatmap
    sns.heatmap(
        deep_colors, 
        ax=axes[1], 
        cmap="viridis", 
        annot=deep_df.astype(float),  # Keep absolute values
        fmt=".0f", 
        cbar=True,
        vmin=0,
        vmax=1
    )
    axes[1].set_title(f"Intersecting BSJs (DEEP, {name})")

    name_clean = name.replace(" ", "_")
    if name_clean == "polyA_vs_Total":
        ylabel, xlabel = "polyA BED", "Total BED"
    elif name_clean == "polyA":
        ylabel, xlabel = "polyA BED", "polyA BED"
    else:
        ylabel, xlabel = "Total BED", "Total BED"

    for ax in axes:
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    plt.tight_layout()
    plt.savefig(f"{name_clean}_int_heat.png", dpi=300)

def _plot_inter_heatmaps(nih_df, deep_df, name):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # NIH heatmap
    sns.heatmap(
        nih_df.astype(float), 
        ax=axes[0], 
        cmap="viridis", 
        annot=True, 
        fmt=".3f", 
        cbar=True
    )
    axes[0].set_title(f"Intersecting BSJs (GSE138734, {name})")

    # DEEP heatmap
    sns.heatmap(
        deep_df.astype(float), 
        ax=axes[1], 
        cmap="viridis", 
        annot=True, 
        fmt=".3f", 
        cbar=True
    )
    axes[1].set_title(f"Intersecting BSJs  (DEEP, {name})")

    name = name.replace(" ", "_")
    if name == "polyA_vs_Total":
        axes[0].set_ylabel("polyA BED")
        axes[0].set_xlabel("Total BED")
        axes[1].set_ylabel("polyA BED")
        axes[1].set_xlabel("Total BED")
    elif name == "polyA":
        axes[0].set_ylabel("polyA BED")
        axes[0].set_xlabel("polyA BED")
        axes[1].set_ylabel("polyA BED")
        axes[1].set_xlabel("polyA BED")
    else:
        axes[0].set_ylabel("Total BED")
        axes[0].set_xlabel("Total BED")
        axes[1].set_ylabel("Total BED")
        axes[1].set_xlabel("Total BED")

    axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
    axes[1].set_yticklabels(axes[1].get_yticklabels(), rotation=0)

    plt.tight_layout()
    plt.savefig(f"{name}_int_heat.png", dpi = 300)

def plot_comp_heatmaps(nih_df, deep_df, name):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # NIH heatmap
    sns.heatmap(
        nih_df.astype(float), 
        ax=axes[0], 
        cmap="viridis", 
        annot=True, 
        fmt=".3f", 
        cbar=True,
        vmin=0,
        vmax=1

    )
    axes[0].set_title(f"Percentage of Overlapping BSJ \n Between BED Files (GSE138734, {name})")

    # DEEP heatmap
    sns.heatmap(
        deep_df.astype(float), 
        ax=axes[1], 
        cmap="viridis", 
        annot=True, 
        fmt=".3f", 
        cbar=True,
        vmin=0,
        vmax=1

    )
    axes[1].set_title(f"Percentage of Overlapping BSJ \n Between BED Files (DEEP, {name})")

    name = name.replace(" ", "_")
    if name == "polyA_vs_Total":
        axes[0].set_ylabel("Poly(A) BED")
        axes[0].set_xlabel("Total BED")
        axes[1].set_ylabel("Poly(A) BED")
        axes[1].set_xlabel("Total BED")
    elif name == "polyA":
        axes[0].set_ylabel("Poly(A) BED")
        axes[0].set_xlabel("Poly(A) BED")
        axes[1].set_ylabel("Poly(A) BED")
        axes[1].set_xlabel("Poly(A) BED")
    else:
        axes[0].set_ylabel("Total BED")
        axes[0].set_xlabel("Total BED")
        axes[1].set_ylabel("Total BED")
        axes[1].set_xlabel("Total BED")

    axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0)
    axes[1].set_yticklabels(axes[1].get_yticklabels(), rotation=0)

    plt.tight_layout()
    out = os.path.join(PLOT_OUT, f"{name}_comp_heat.png")
    plt.savefig(out, dpi = 300)

# I. Plot amount of bsjs in merged beds per dataset
def plot_merged_all():
    NIH_polya = pd.read_csv("../data/GSE138734/merged_all/polya/polya_merged_all.bed"
    , sep = "\t", 
    header=None,
    names=["chrom", "start", "end", "random", "scores", "stand"]
    )
    NIH_total = pd.read_csv("../data/GSE138734/merged_all/total/total_merged_all.bed"
    , sep = "\t", 
    header=None,
    names=["chrom", "start", "end", "random", "scores", "stand"]
    )
    
    DEEP_polya = pd.read_csv("../data/DEEP/merged_all/polya/polya_merged_all.bed"
    , sep = "\t", 
    header=None,
    names=["chrom", "start", "end", "random", "scores", "stand"]
    )
    DEEP_total = pd.read_csv("../data/DEEP/merged_all/total/total_merged_all.bed"
    , sep = "\t", 
    header=None,
    names=["chrom", "start", "end", "random", "scores", "stand"]
    )

    counts = {
        "DEEP": {
            "total": len(DEEP_total),
            "polya": len(DEEP_polya)
        },
        "GSE138734": {
            "total": len(NIH_total),
            "polya": len(NIH_polya)
        }
    }

    df_counts = pd.DataFrame(counts).T.reset_index().rename(columns={"index": "group"})

    fig, ax = plt.subplots(figsize=(8, 4.5))
    plt.subplots_adjust(top=0.85)  # move title higher up

    fig, ax = plt.subplots(figsize=(4.5, 8))
    plt.subplots_adjust(top=0.85)  # move title higher up

    for i, row in df_counts.iterrows():
        ax.plot([i, i], [row["total"], row["polya"]], color="gray", lw=3, zorder=1)

        ax.scatter(i, row["total"], color=palette["total"], s=200, label="Total" if i == 0 else "", edgecolors="black", linewidths=2)
        ax.scatter(i, row["polya"], color=palette["polya"], s=200, label="Poly(A)" if i == 0 else "", edgecolors="black", linewidths=2)

        ax.text(i+0.05, row["total"], f"{row['total']}", va="center", ha="left", fontsize=11)
        ax.text(i+0.05, row["polya"], f"{row['polya']}", va="center", ha="left", fontsize=11)

    ax.set_xticks(range(len(df_counts)))
    ax.set_xticklabels(df_counts["group"], ha="center")

    ax.yaxis.grid(True, linestyle="--", linewidth=1, color="lightgray")
    ax.set_axisbelow(True)

    ax.set_ylabel("Amount")
    ax.set_title("Unique BSJs per Dataset\n(Merged Across All Samples)")
    ax.set_ylim(0, max(df_counts[["total", "polya"]].max()) * 1.1)


    sns.despine()
    ax.legend()
    plt.tight_layout()

    out = os.path.join(PLOT_OUT, "overview_1.png")
    plt.savefig(out, dpi=300)


def plot_dataset(base_dir, output_file, dataset):
    base_dir = Path(base_dir)
    output_file = os.path.join(PLOT_OUT, output_file)

    counts = []
    for condition in ["total", "polya"]:
        files = (base_dir / condition).glob("*.bed")
        for f in files:
            tool = f.name.split("_")[0]  # extract tool name (before first "_")
            if tool == "find":
                tool = "find_circ"
            df = pd.read_csv(
                f,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", "random", "scores", "stand"]
            )
            counts.append({"tool": tool, "condition": condition, "count": len(df)})

    df_counts = pd.DataFrame(counts)
    print(df_counts.head())
    name_map = { 'segemehl': 'Segemehl',
    'dcc': 'DCC',
    'ciriquant': 'CIRIquant',
    'circexplorer2': 'CIRCexplorer2',
    'find_circ':'find_circ'
    }
    
    df_counts['tool'] =  df_counts["tool"].replace(name_map)

    print(df_counts.head())

    # Pivot so each tool has total and polya side by side
    df_wide = df_counts.pivot(index="tool", columns="condition", values="count").reset_index()

    fig, ax = plt.subplots(figsize=(9, 5))
    plt.subplots_adjust(top=0.85)  # move title higher up

    for i, row in df_wide.iterrows():
        # gray connector line
        ax.plot([row["total"], row["polya"]], [i, i], color="gray", lw=3, zorder=1)
        # scatter points with palette colors
        ax.scatter(row["total"], i, color=palette["total"], s=200, label="Total" if i == 0 else "", edgecolors="black", linewidths=2)
        ax.scatter(row["polya"], i, color=palette["polya"], s=200, label="Poly(A)" if i == 0 else "", edgecolors="black", linewidths=2)

        #ax.text(row["total"]-300, i+0.1, f"{row['total']}", va="center", ha="left", fontsize=12)
        #ax.text(row["polya"]+150, i+0.1, f"{row['polya']}", va="center", ha="right", fontsize=12)

    # Formatting
    dataset = "GSE138734" if dataset == "NIH" else dataset
    ax.set_yticks(range(len(df_wide)))
    ax.set_yticklabels(df_wide["tool"])

    ax.xaxis.grid(True, linestyle="--", linewidth=1, color="lightgray")
    ax.set_axisbelow(True)

    ax.set_xlabel("Amount")
    ax.set_title(f"Unique BSJs per Tool in {dataset} \n (Merged Across all Samples)")
    ax.legend()
    sns.despine()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def build_bsj_union(base_dir, condition="polya"):
    base_dir = Path(base_dir) / condition

    dfs = []
    for f in base_dir.glob("*.bed"):
        tool = f.name.split("_")[0]
        if tool == "find":
            tool = "find_circ"
        df = pd.read_csv(
            f,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "random", "scores", "strand"]
        )
        df = df[["chrom", "start", "end", "strand"]].copy()
        df["tool"] = tool
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def plot_tool_counts(base_dir, condition="polya", dataset="NIH", output_file="tool_counts.png"):
    df_union = build_bsj_union(base_dir, condition=condition)

    grouped = (
        df_union.groupby(["chrom", "start", "end", "strand"])["tool"]
        .nunique()
        .reset_index(name="tool_count")
    )

    counts = grouped["tool_count"].value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(6, 4))
    sns.barplot(x=counts.index, y=counts.values, ax=ax, color="#66C2A5")

    ax.set_xlabel("Number of Tools Detecting BSJ")
    ax.set_ylabel("Number of BSJs")
    ax.set_title(f"BSJ Detection by Tools in {dataset} ({condition})")

    for i, val in enumerate(counts.values):
        ax.text(i, val, str(val), ha="center", va="bottom")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def plot_upset(base_dir, condition, dataset, output_file):
    """
    Create an upset plot showing overlap of BSJs across tools.
    """

    output_file = os.path.join(PLOT_OUT, output_file)
    df_union = build_bsj_union(base_dir, condition=condition)
    name_map = { 'segemehl': 'Segemehl',
    'dcc': 'DCC',
    'ciriquant': 'CIRIquant',
    'circexplorer2': 'CIRCexplorer2',
    'find_circ':'find_circ'
    }
    
    df_union['tool'] =  df_union["tool"].replace(name_map)

    # Group by BSJ coordinates and collect tools
    grouped = (
        df_union.groupby(["chrom", "start", "end", "strand"])["tool"]
        .apply(set)
        .reset_index()
    )

    name_map = { 'segemehl': 'Segemehl',
    'dcc': 'DCC',
    'ciriquant': 'CIRIquant',
    'circexplorer2': 'CIRCexplorer2',
    'find_circ':'find_circ'
    }
    
    grouped['tool'] =  grouped["tool"].replace(name_map)

    # Convert to list for upsetplot
    memberships = grouped["tool"].apply(lambda s: tuple(sorted(s)))

    # Collapse identical tool sets by counting occurrences
    membership_counts = memberships.value_counts()

    # Build upset data
    upset_data = from_memberships(membership_counts.index, data=membership_counts.values)

    # Plot
    fig = plt.figure(figsize=(10, 6))
    upset = UpSet(upset_data, show_counts=True, sort_by="cardinality")
    upset.plot(fig=fig)
    plt.subplots_adjust(top=0.75)  # move title higher up
    for ax in fig.axes:
        if ax.get_ylabel() and "Intersection" in ax.get_ylabel():
            for text in ax.texts:
                text.set_rotation(90)
                text.set_ha('center')  # horizontal alignment
                x, y = text.get_position()
                text.set_position((x, y*1.1)) 
    
    cond_title = "Poly(A)" if condition == "polya" else "Total"
    dataset = "GSE138734" if dataset == "NIH" else dataset
    plt.suptitle(f"UpSet Plot of BSJs Across Tools in {dataset} ({cond_title})")
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_tool_counts_mirrored(base_dir, dataset, output_file):
    """
    Create a mirrored horizontal bar plot:
    left = total, right = polya
    categories = number of tools detecting a BSJ (1–5).
    """
    output_file = os.path.join(PLOT_OUT, output_file)
    df_total = build_bsj_union(base_dir, condition="total")
    df_polya = build_bsj_union(base_dir, condition="polya")

    def count_tools(df):
        grouped = (
            df.groupby(["chrom", "start", "end", "strand"])["tool"]
            .nunique()
            .reset_index(name="tool_count")
        )
        return grouped["tool_count"].value_counts().sort_index()

    counts_total = count_tools(df_total)
    counts_polya = count_tools(df_polya)


    # Align indices so both have 1–5
    all_idx = range(1, 6)
    counts_total = counts_total.reindex(all_idx, fill_value=0)
    counts_polya = counts_polya.reindex(all_idx, fill_value=0)

    fig, ax = plt.subplots(figsize=(9, 5))

    # Draw bars (total goes negative just for plotting)
    ax.barh(all_idx, -counts_total.values, color=palette["total"], label="Total", edgecolor="black", linewidth=1.5)
    ax.barh(all_idx, counts_polya.values, color=palette["polya"], label="Poly(A)", edgecolor="black", linewidth=1.5)

    # Formatting
    ax.set_yticks(all_idx)
    ax.set_yticklabels([f"# {i}" for i in all_idx])
    ax.invert_yaxis()  # so 1 tool is at top

    # Symmetric x-axis with 0 in the middle
    max_count = max(counts_total.max(),counts_polya.max() )
    ax.set_xlim(-max_count * 1.1,  counts_polya.max()* 2.1)

    # Custom x-ticks: show positive values both sides
    # xticks = np.linspace(-max_count, counts_polya.max(), 5, dtype=int)
    xticks = np.unique(np.concatenate([np.linspace(-max_count, counts_polya.max(), 5, dtype=int), [0]]))
    ax.set_xticks(xticks)
    ax.set_xticklabels([abs(x) for x in xticks])

    # Add vertical line at x=0
    ax.axvline(0, color="black", linewidth=2)

    # Titles and labels
    dataset = "GSE138734" if dataset == "NIH" else dataset
    plt.subplots_adjust(top=0.85)  # move title higher up
    ax.set_title(f"Distribution of Unique BSJs by Number \n of Detecting Tools in {dataset}")
    ax.set_xlabel("Amount")
    ax.set_ylabel("Confidence Level (# Tools)")

    ax.xaxis.grid(True, linestyle="--", linewidth=1, color="lightgray")
    ax.set_axisbelow(True)

    # Add counts as labels (always positive)
    for y, val in zip(all_idx, counts_total.values):
        if val > 0:
            ax.text(-val, y, str(val), va="center", ha="right", color="black", rotation=90)
    for y, val in zip(all_idx, counts_polya.values):
        if val > 0:
            ax.text(val, y, str(val), va="center", ha="left", color="black", rotation=-90)

    ax.legend(loc="lower left")
    sns.despine()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def compare_to_majority(base_dir, dataset, output_file):
    output_file = os.path.join(PLOT_OUT, output_file)
    df_total = build_bsj_union(base_dir, condition="total")
    df_polya = build_bsj_union(base_dir, condition="polya")
    majority_total = (
        df_total.groupby(["chrom", "start", "end", "strand"])["tool"]
        .nunique()
        .reset_index(name="tool_count")
    )

    #name_map = { 'segemehl': 'Segemehl',
    #'dcc': 'DCC',
    #'ciriquant': 'CIRIquant',
    #'circexplorer2': 'CIRCexplorer2',
    #'find_circ':'find_circ'
    #}
    
    #majority_total['tool'] =  majority_total["tool"].replace(name_map)

    majority_total = majority_total[majority_total["tool_count"] >= 4]
    print(majority_total.head())

    df_polya_majority = df_polya.merge(
        majority_total[["chrom", "start", "end", "strand"]],
        on=["chrom", "start", "end", "strand"],
        how="inner"
    )
    print(df_polya_majority)

    grouped = (
        df_polya_majority.groupby(["chrom", "start", "end", "strand"])["tool"]
        .apply(set)
        .reset_index()
    )

    memberships = grouped["tool"].apply(lambda s: tuple(sorted(s)))

    membership_counts = memberships.value_counts()

    upset_data = from_memberships(membership_counts.index, data=membership_counts.values)

    fig = plt.figure(figsize=(14, 8))
    upset = UpSet(upset_data, show_counts=False, sort_by="cardinality")
    upset.plot(fig=fig)
    plt.subplots_adjust(top=0.75)  # move title higher up

    dataset = "GSE138734" if dataset == "NIH" else dataset
    plt.suptitle(f"PolyA BSJs Matching \n TotalRNA Consensus \n ({dataset})")

    plt.savefig(output_file, dpi=300)
    plt.close()


if __name__ == '__main__':
    plot_merged_all()

    plot_dataset("../data/GSE138734/merge_concatenated_beds", "NIH_tools.png", "NIH")
    plot_dataset("../data/DEEP/merge_concatenated_beds", "DEEP_tools.png", "DEEP")

    plot_upset(
        base_dir="../data/GSE138734/merge_concatenated_beds",
        condition="polya",
        dataset="NIH",
        output_file="NIH_polya_upset.png"
    )
    plot_upset(
        base_dir="../data/GSE138734/merge_concatenated_beds",
        condition="total",
        dataset="NIH",
        output_file="NIH_total_upset.png"
    )
    plot_upset(
        base_dir="../data/DEEP/merge_concatenated_beds",
        condition="polya",
        dataset="DEEP",
        output_file="DEEP_polya_upset.png"
    )
    plot_upset(
        base_dir="../data/DEEP/merge_concatenated_beds",
        condition="total",
        dataset="DEEP",
        output_file="DEEP_total_upset.png"
    )
    
    plot_tool_counts_mirrored(
        base_dir="../data/GSE138734/merge_concatenated_beds",
        dataset="NIH",
        output_file="NIH_tool_intersection.png"
    )
    
    plot_tool_counts_mirrored(
        base_dir="../data/DEEP/merge_concatenated_beds",
        dataset="DEEP",
        output_file="DEEP_tool_intersection.png"
    )

    compare_to_majority("../data/DEEP/merge_concatenated_beds", "DEEP", "DEEP_gt.png")
    compare_to_majority("../data/GSE138734/merge_concatenated_beds", "NIH", "NIH_gt.png")

    tools = ["circexplorer2", "dcc", "segemehl", "ciriquant", "find_circ"] 
    n, d, n_p, d_p, n_t, d_t  =  run_jaccard_matrix(tools)
    plot_jaccard_heatmaps(n,d, "polyA vs Total")
    plot_jaccard_heatmaps(n_p,d_p, "polyA")
    plot_jaccard_heatmaps(n_t,d_t, "total")
    
    n, d, n_p, d_p, n_t, d_t  =  run_compare_matrix(tools)
    plot_comp_heatmaps(n,d, "polyA vs Total")
    plot_comp_heatmaps(n_p,d_p, "polyA")
    plot_comp_heatmaps(n_t,d_t, "total")
