import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import sys

MAIN_DATA_DIR = sys.argv[1]
RRNA_DIR = sys.argv[2]

READS_AMOUNT = dict()
RRNA_READS_AMOUNT = dict()


dataset_name = MAIN_DATA_DIR.split("/")[-1]
if len(dataset_name) == 0:
    dataset_name = MAIN_DATA_DIR.split("/")[-2]
    
def read_bam_meta(path_to_meta, read_dict):
    with open(path_to_meta) as f:
        lines = f.readlines()
        for line in lines:
            if "File" in line:
                continue
            sample = "_".join(line.split("\t")[0].split("/")[-1].split(".")[0].split("_")[0:2])
            amount = int(line.split("\t")[-1].strip("\n"))
            read_dict[sample] = amount

polya_bam_amount = os.path.join(MAIN_DATA_DIR, "polya", "bam_meta", "counts.txt")
total_bam_amount = os.path.join(MAIN_DATA_DIR, "total", "bam_meta", "counts.txt")
read_bam_meta(polya_bam_amount, READS_AMOUNT)
read_bam_meta(total_bam_amount, READS_AMOUNT)

polya_bam_rrna_amount = os.path.join(MAIN_DATA_DIR, "polya", "rrna_reads", "rrna_spanning_reads.tsv")
total_bam_rrna_amount = os.path.join(MAIN_DATA_DIR, "total", "rrna_reads", "rrna_spanning_reads.tsv")
read_bam_meta(total_bam_rrna_amount, RRNA_READS_AMOUNT)
read_bam_meta(polya_bam_rrna_amount, RRNA_READS_AMOUNT)

if READS_AMOUNT.keys() != RRNA_READS_AMOUNT.keys():
    print(f"Samples missmatching between bam_meta files and rrna_spanning_reads files: \n BAM_META: {READS_AMOUNT.keys()} \n vs. \n RRNA: {RRNA_READS_AMOUNT.keys()}")

# Build df_pct from the dicts
rows = []
for key, rrna_count in RRNA_READS_AMOUNT.items():
    total_count = READS_AMOUNT.get(key)
    if total_count is None:
        print(f"Warning: {key} not found in READS_AMOUNT")
        continue
    parts = key.rsplit("_", 1)
    sample_id = parts[0]
    rna_type = parts[1].lower()   
    pct = (rrna_count / total_count) * 100
    rows.append({"Sample": sample_id, "Type": rna_type, "Percent_rRNA": pct})

df_pct = pd.DataFrame(rows)
df_pct["Sample"] = df_pct["Sample"].astype(int)
df_pct = df_pct.sort_values("Sample")
df_pct["Sample"] = df_pct["Sample"].astype(str)

sns.set_theme(style="ticks", context="paper", palette="colorblind")
sns.set_context(
    "paper",
    rc={
        "font.size": 15,
        "axes.titlesize": 17,
        "axes.labelsize": 15,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 12,
        "lines.linewidth": 1.5,
        "lines.markersize": 6,
    },
)
plt.rcParams["figure.figsize"] = (8, 4.5)
cb_palette = sns.color_palette("colorblind")
palette = {
    "total": cb_palette[0],
    "polya": cb_palette[1],
}

fig, ax = plt.subplots()
sns.barplot(
    data=df_pct,
    x="Sample",
    y="Percent_rRNA",
    hue="Type",
    hue_order=["total", "polya"],
    palette=palette,
    ax=ax,
    edgecolor="black"
)
ax.set_ylabel("% rRNA Reads")
ax.set_xlabel("Sample")
handles, labels = ax.get_legend_handles_labels()
label_map = {"total": "Total", "polya": "Poly(A)"}
ax.legend(handles, [label_map.get(l, l) for l in labels])
ax.yaxis.grid(True, linestyle="--", linewidth=1, color="lightgray")
sns.despine()
ax.set_title(f"rRNA Contamination per Sample in {dataset_name}")
plt.tight_layout()
plt.savefig(os.path.join(RRNA_DIR, "rrna_percentage.png"), dpi=300, bbox_inches="tight")
plt.show()