import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
from matplotlib.lines import Line2D

import os
import json

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


def plot_corr(total_json_path, polya_json_path):
    with open(total_json_path, "r") as f:
        total_data = json.load(f)

    with open(polya_json_path, "r") as f:
        polya_data = json.load(f)

    # bar_plot(total_data, "total")
    # bar_plot(polya_data, "polya")
    tissue = polya_json_path.split(".")[0].split("_")[1]
    path = f"{tissue}_dumbbell.corr_bar.png"
    dumbell_plot(total_data, polya_data, path, tissue)


def bar_plot(data, data_origin):
    sorted_items = sorted(data.items(), key=lambda x: x[1])
    tools = [tool for tool, _ in sorted_items]
    corr_values = [val for _, val in sorted_items]

    plt.figure(figsize=(6, 4))
    plt.bar(tools, corr_values, edgecolor="black")
    plt.ylabel("Correlation")
    plt.title(
        f"Correlation between number of BSJ reads\nand number of rRNA-spanning reads\n(source: {data_origin})"
    )
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.savefig(f"{data_origin}.corr_bar.png", dpi=300)
    plt.close()


def dumbell_plot(total, polya, output_file, tissue):
    cb_palette = sns.color_palette("colorblind")
    palette = {"total": cb_palette[0], "polya": cb_palette[1]}

    df = pd.DataFrame(
        {
            "Tool": list(total.keys()),
            "total_r": [total[k]["r"] for k in total.keys()],
            "total_p": [total[k]["p"] for k in total.keys()],
            "polya_r": [polya[k]["r"] for k in total.keys()],
            "polya_p": [polya[k]["p"] for k in total.keys()],
        }
    )
    name_map = { 'segemehl_filtered': 'Segemehl',
    'dcc_filtered': 'DCC',
    'ciriquant_filtered': 'CIRIquant',
    'circexplorer2_filtered': 'CIRCexplorer2',
    'find_circ_filtered':'find_circ'
    }
    df['Tool'] = df['Tool'].replace(name_map)
    df = df.sort_values("Tool", ascending=True)

    fig, ax = plt.subplots(figsize=(6, 4))

    def p_to_stars(p):
        if p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"

    # Plot dumbbell lines
    for i, row in df.iterrows():
        ax.plot(
            [row["total_r"], row["polya_r"]],
            [row["Tool"], row["Tool"]],
            color="gray",
            linestyle="-",
            linewidth=3,
            zorder=1,
        )

    # Scatter points
    total_scatter = ax.scatter(
        df["total_r"], df["Tool"], label="Total", s=50,
        edgecolor="black", linewidths=2, facecolor=palette["total"], zorder=2
    )
    polya_scatter = ax.scatter(
        df["polya_r"], df["Tool"], label="Poly(A)", s=50,
        edgecolor="black", linewidths=2, facecolor=palette["polya"], zorder=2
    )

        # Map tool names to their numeric y positions
    yticks = ax.get_yticks()       # numeric positions
    ylabels = [tick.get_text() for tick in ax.get_yticklabels()]
    tool_to_y = {tool: y for tool, y in zip(ylabels, yticks)}
    

    for i, row in df.iterrows():
        y_pos = tool_to_y[row["Tool"]]
        ax.text(row["total_r"], y_pos + 0.14, p_to_stars(row["total_p"]),
                va="center", ha="right", fontsize=10, fontweight="bold")
        ax.text(row["polya_r"], y_pos + 0.14, p_to_stars(row["polya_p"]),
                va="center", ha="left", fontsize=10, fontweight="bold")

    tissue = "DEEP" if tissue == "mixed" else tissue

    ax.set_xlabel("Pearson Correlation", fontsize=15)
    ax.set_title(f"{tissue}",
                 fontsize=13, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    #ax.legend(handles=legend_handles, loc="best")
    ax.set_xlim(-1.2, 1.2)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)

    fig, ax = plt.subplots(figsize=(6, 2))

    cb_palette = sns.color_palette("colorblind")
    total_color = cb_palette[0]
    polya_color = cb_palette[1]

    legend_handles = [
        Line2D([0], [0], marker='o', color='w', label='Total',
               markerfacecolor=total_color, markeredgecolor='black', markersize=12),
        Line2D([0], [0], marker='o', color='w', label='Poly(A)',
               markerfacecolor=polya_color, markeredgecolor='black', markersize=12),
        Line2D([0], [0], color='none', marker='*', markerfacecolor='black',
               markersize=12, label='Significance: * p<0.05 / ** p<0.01 / *** p<0.001')
    ]

    ax.legend(handles=legend_handles, loc='center', frameon=False, ncol=1)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig("dumbell_legend.png", dpi=300, bbox_inches='tight')
    plt.close()

def _dumbell_plot(total, polya, output_file, tissue):
    cb_palette = sns.color_palette("colorblind")

    palette = {
        "total": cb_palette[0],   # 1st color from seaborn's colorblind palette
        "polya": cb_palette[1],       # keep your green
    }
    print(total)
    print(polya)

    df = pd.DataFrame(
        {
            "Tool": list(total.keys()),
            "total": [total[k]["r"] for k in total.keys()],
            "polya": [polya[k]["r"] for k in total.keys()],
        }
    )

    #df["mean"] = (df["total"]["r"] + df["polya"]["r"]) / 2
    df = df.sort_values("Tool", ascending=True)


    fig, ax = plt.subplots(figsize=(10, 6))

    for i, row in df.iterrows():
        ax.plot(
            [row["total"], row["polya"]],
            [row["Tool"], row["Tool"]],
            color="gray",
            linestyle="-",
            linewidth=3,
            zorder=1,
        )

    total_scatter = ax.scatter(
        df["total"],
        df["Tool"],
        label="total",
        s=100,
        edgecolor="black",
        linewidths=2,
        zorder=2,
        facecolor=palette["total"],
    )

    polya_scatter = ax.scatter(
        df["polya"],
        df["Tool"],
        label="polya",
        s=100,
        edgecolor="black",
        zorder=2,
        linewidths=2,
        facecolor=palette["polya"],
    )

    ax.set_xlabel("Pearson Correlation", fontsize=12)
    ax.set_title(f"Dumbbell Plot Comparing rRNA Correlation in {tissue}", fontsize=14, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(handles=[total_scatter, polya_scatter], loc="best")
    ax.set_xlim(-1.2,1.2)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)

    return fig, ax


if __name__ == "__main__":
    total_json_path = sys.argv[1]
    polya_json_path = sys.argv[2]

    plot_corr(total_json_path, polya_json_path)
