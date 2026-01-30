import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

types  = ["LiHe", "BlCM", "BlEM", "BlTN", "mixed"]

for t in types:
    # Read the data
    total = pd.read_csv(
        f"./total_{t}.rrnareads.txt", sep=",", header=None, names=["Sample", "Reads"]
    )
    total["Sample"] = total["Sample"].str.replace(
        ".Aligned.out.bam.total.featurecounts.txt.summary", "", regex=False
    )
    total["Type"] = "Total"

    poly = pd.read_csv(
        f"./polya_{t}.rrnareads.txt", sep=",", header=None, names=["Sample", "Reads"]
    )
    poly["Sample"] = poly["Sample"].str.replace(".Aligned.out.bam.polya.featurecounts.txt.summary", "", regex=True)
    poly["Type"] = "Poly(A)"

    cb_palette = sns.color_palette("colorblind")
    palette = {"Total": cb_palette[0], "Poly(A)": cb_palette[1]}

    # Combine dataframes
    df_long = pd.concat([total, poly], ignore_index=True)

    # Set Type as categorical with specific order

    # Sort Sample column and set as categorical for proper x-axis order
    df_long["Type"] = pd.Categorical(
        df_long["Type"], categories=["Total", "Poly(A)"], ordered=True
    )
    # sample_order = sorted(
    #     df_long["Sample"].unique()
    # )  # or use key=lambda x: int(x) if samples are numeric
    # df_long["Sample"] = pd.Categorical(
    #     df_long["Sample"], categories=sample_order, ordered=True
    # )


    df_long["Sample"] = df_long["Sample"].str.replace(
        "_.*", "", regex=True
    )

    sample_order = sorted(df_long["Sample"].unique(), key=lambda x: int(x))
    df_long["Sample"] = pd.Categorical(
                df_long["Sample"], categories=sample_order, ordered=True
                )

    # Plot
    plt.figure(figsize=(6, 4))
    ax = sns.barplot(
        data=df_long,
        x="Sample",
        y="Reads",
        hue="Type",
        palette=palette,
        hue_order=["Total", "Poly(A)"], 
        dodge=True,
        edgecolor="black",
        linewidth=1.2,
    )

    ax.yaxis.grid(True, linestyle='--', color='gray', alpha=0.7)

    plt.title(f"RRNA Spanning Reads in DEEP", fontsize=16)
    plt.xlabel("Sample")
    plt.ylabel("Read Count (CPM)")
    plt.legend(title="Type")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{t}_rrna_reads_per_sample_cpm.png", dpi=300)
    plt.close()




total = pd.read_csv(
            "./total.txt", sep="\t",
            )
total["Sample"] = total["File"].str.split("/").str[-1]
total["Type"] = "total"

polya = pd.read_csv(
            "./polya.txt", sep="\t",
            )
polya["Sample"] = polya["File"].str.split("/").str[-1]
polya["Type"] = "polyA"


df_long = pd.concat([total, polya], ignore_index=True)

df_long["Type"] = pd.Categorical(
            df_long["Type"], categories=["total", "polyA"], ordered=True
            )


df_long["Sample"] = df_long["Sample"].str.replace(
            "_.*", "", regex=True
            )

sample_order = sorted(df_long["Sample"].unique(), key=lambda x: int(x))
df_long["Sample"] = pd.Categorical(
                    df_long["Sample"], categories=sample_order, ordered=True
                                )


palette = {
            "total": "#FC8D62",  # orange
                "polyA": "#66C2A5",  # green
                }

# Plot
plt.figure(figsize=(12, 6))
sns.barplot(
            data=df_long,
                x="Sample",
                    y="Read_Count",
                        hue="Type",
                            hue_order=["polyA", "total"],
                                palette=palette,
                                    dodge=True,
                                        edgecolor="black",
                                            linewidth=1.2,
                                            )

plt.title(f"Raw Read Amount by Sample and Type", fontsize=16)
plt.xlabel("Sample")
plt.ylabel("Read Count")
plt.legend(title="Type")
plt.tight_layout()
plt.savefig("rrna_reads_per_sample_raw.png", dpi=300)
plt.close()

