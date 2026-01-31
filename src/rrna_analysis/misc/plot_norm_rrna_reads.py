import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import sys

MAIN_DIR = sys.argv[1]
RRNA_DIR = sys.argv[2]

# Read total counts
total_reads_path = os.path.join(MAIN_DIR, "total/bam_meta/counts.txt")
df_total_total = pd.read_csv(total_reads_path, sep="\t")

polya_reads_path = os.path.join(MAIN_DIR, "polya/bam_meta/counts.txt")
df_total_polya = pd.read_csv(polya_reads_path, sep="\t")

total_total_counts = {
    os.path.basename(f).split(".")[0].split("_")[0]: count
    for f, count in zip(df_total_total["File"], df_total_total["Read_Count"])
}

total_polya_counts = {
    os.path.basename(f).split(".")[0].split("_")[0]: count
    for f, count in zip(df_total_polya["File"], df_total_polya["Read_Count"])
}

# Read rRNA CPM counts
total_rrna_path = os.path.join(RRNA_DIR, "total.rrnareads_cpm.txt")
df_rrna_total = pd.read_csv(total_rrna_path, header=None)
df_rrna_total.columns = ["File", "rRNA_CPM"]

polya_rrna_path = os.path.join(RRNA_DIR, "polya.rrnareads_cpm.txt")
df_rrna_polya = pd.read_csv(polya_rrna_path, header=None)
df_rrna_polya.columns = ["File", "rRNA_CPM"]

rrna_total_cpm = {
    os.path.basename(f).split(".")[0].split("_")[0]: cpm
    for f, cpm in zip(df_rrna_total["File"], df_rrna_total["rRNA_CPM"])
}

rrna_polya_cpm = {
    os.path.basename(f).split(".")[0].split("_")[0]: cpm
    for f, cpm in zip(df_rrna_polya["File"], df_rrna_polya["rRNA_CPM"])
}

# Convert CPM back to raw counts and calculate percentages
pct_data = []

for sample_id in total_total_counts.keys():
    if sample_id in rrna_total_cpm:
        total_count = total_total_counts[sample_id]
        rrna_cpm = rrna_total_cpm[sample_id]
        # Convert CPM back to raw count
        rrna_raw = (rrna_cpm * total_count) / 1e6
        # Calculate percentage
        pct_rrna = (rrna_raw / total_count) * 100
        
        pct_data.append({
            'Sample': sample_id,
            'Type': 'total',
            'Percent_rRNA': pct_rrna
        })

for sample_id in total_polya_counts.keys():
    if sample_id in rrna_polya_cpm:
        total_count = total_polya_counts[sample_id]
        rrna_cpm = rrna_polya_cpm[sample_id]
        # Convert CPM back to raw count
        rrna_raw = (rrna_cpm * total_count) / 1e6
        # Calculate percentage
        pct_rrna = (rrna_raw / total_count) * 100
        
        pct_data.append({
            'Sample': sample_id,
            'Type': 'polya',
            'Percent_rRNA': pct_rrna
        })

df_pct = pd.DataFrame(pct_data)

# Sort samples numerically
def sort_key(sample):
    """Extract numeric part for sorting"""
    import re
    match = re.search(r'(\d+)', str(sample))
    return int(match.group(1)) if match else 0

df_pct['sort_key'] = df_pct['Sample'].apply(sort_key)
df_pct = df_pct.sort_values('sort_key')

# Set up plot style
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

# Plot
fig, ax = plt.subplots()
sns.barplot(
    data=df_pct, 
    x='Sample', 
    y='Percent_rRNA', 
    hue='Type',
    hue_order=['total', 'polya'],
    palette=palette,
    ax=ax,
    edgecolor= "black"

)
ax.set_ylabel('% rRNA Reads')
ax.set_xlabel('Sample')
handles, labels = ax.get_legend_handles_labels()
label_map = {'total': 'Total', 'polya': 'Poly(A)'}
new_labels = [label_map.get(label, label) for label in labels]
ax.legend(handles, new_labels)
ax.yaxis.grid(True, linestyle="--", linewidth=1, color="lightgray")
sns.despine()
ax.set_title(f'rRNA Contamination by Sample in {RRNA_DIR.split("/")[-1]}')
plt.tight_layout()
plt.savefig(os.path.join(RRNA_DIR, 'rrna_percentage.png'), dpi=300, bbox_inches='tight')
plt.show()