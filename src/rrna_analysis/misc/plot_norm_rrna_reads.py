import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os


df_total_total = pd.read_csv(
    "./total.txt", sep="\t")  # total amount of reads per sample
df_total_polya = pd.read_csv("./polya.txt", sep="\t")

total_total_counts = {
    os.path.basename(f).split(".")[0].split("_")[0]: count
    for f, count in zip(df_total_total["File"], df_total_total["Read_Count"])
}

total_polya_counts = {
    os.path.basename(f).split(".")[0].split("_")[0]: count
    for f, count in zip(df_total_polya["File"], df_total_polya["Read_Count"])
}

df_rrna_total = pd.read_csv("./total.rrnareads.txt", header=None)
df_rrna_polya = pd.read_csv("./polya.rrnareads.txt", header=None)
df_rrna_total.columns = ["File", "rRNA_Reads"]
df_rrna_polya.columns = ["File", "rRNA_Reads"]

rrna_total_counts = {
    os.path.basename(f).split(".")[0].split("_")[0]: count
    for f, count in zip(df_rrna_total["File"], df_rrna_total["rRNA_Reads"])
}

rrna_polya_counts = {
    os.path.basename(f).split(".")[0].split("_")[0]: count
    for f, count in zip(df_rrna_polya["File"], df_rrna_polya["rRNA_Reads"])
}

norm_rrna_total = dict()
norm_rrna_polya = dict()

for sample in rrna_total_counts.keys():
    norm = (rrna_total_counts[sample] / total_total_counts[sample]) * 1000000
    norm_rrna_total[sample] = norm

for sample in rrna_polya_counts.keys():
    norm = (rrna_polya_counts[sample] / total_polya_counts[sample]) * 1000000
    norm_rrna_polya[sample] = norm

df_total = pd.DataFrame(
    {
        "Sample": list(norm_rrna_total.keys()),
        "RPM": list(norm_rrna_total.values()),
        "Type": "total",
    }
)

df_polya = pd.DataFrame(
    {
        "Sample": list(norm_rrna_polya.keys()),
        "RPM": list(norm_rrna_polya.values()),
        "Type": "polyA",
    }
)

palette = {
    "total": "#FC8D62",  # orange
    "polyA": "#66C2A5",  # green
}
# Combine both
df_all = pd.concat([df_total, df_polya], ignore_index=True)

# Optional: sort samples alphabetically or by RPM
#df_all.sort_values(by="Sample", inplace=True)

# Optional: sort samples alphabetically or by RPM
sample_order = sorted(df_all["Sample"].unique(), key=lambda x: int(x))
df_all["Sample"] = pd.Categorical(df_all["Sample"], categories=sample_order, ordered=True)

# Plot
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")

# Bar plot with borders (edgecolor)
barplot = sns.barplot(
    data=df_all,
    x="Sample",
    y="RPM",
    hue="Type",
    edgecolor="black",
    linewidth=1.2,
    palette=palette,
)

plt.title("Normalized rRNA Reads (RPM)")
plt.ylabel("Reads Per Million (RPM)")
plt.xlabel("Sample")
plt.tight_layout()
plt.legend(title="Type")
plt.tight_layout()
plt.savefig("normalized_rrna_read_amount.png", dpi=300)
plt.close()


df_total = pd.DataFrame(
    {
        "Sample": list(total_total_counts.keys()),
        "Amount": list(total_total_counts.values()),
        "Type": "total",
    }
)

df_polya = pd.DataFrame(
    {
        "Sample": list(total_polya_counts.keys()),
        "Amount": list(total_polya_counts.values()),
        "Type": "polyA",
    }
)


palette = {
    "total": "#FC8D62",  # orange
    "polyA": "#66C2A5",  # green
}
# Combine both
df_all = pd.concat([df_total, df_polya], ignore_index=True)


# Optional: sort samples alphabetically or by RPM
sample_order = sorted(df_all["Sample"].unique(), key=lambda x: int(x))
df_all["Sample"] = pd.Categorical(df_all["Sample"], categories=sample_order, ordered=True)


# Plot
plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")

# Bar plot with borders (edgecolor)
barplot = sns.barplot(
    data=df_all,
    x="Sample",
    y="Amount",
    hue="Type",
    edgecolor="black",
    linewidth=1.2,
    palette=palette,)

plt.title("Total Amount of reads in BAM file")
plt.ylabel("Reads")
plt.xlabel("Sample")
plt.tight_layout()
plt.legend(title="Type")
plt.tight_layout()
plt.savefig("bam_read_amount.png", dpi=300)
plt.close()
