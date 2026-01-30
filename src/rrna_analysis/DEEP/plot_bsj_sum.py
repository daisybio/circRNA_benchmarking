import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Custom palette
palette = {
    "total": "#FC8D62",  # orange
    "polyA": "#66C2A5",  # green
}

# Input files
polya_file = "./polya_mixed.bsj_amount.json"
total_file = "./total_mixed.bsj_amount.json"

# Load JSONs
with open(polya_file) as f:
    polya_data = json.load(f)

with open(total_file) as f:
    total_data = json.load(f)

tools = polya_data.keys()

# Plot per tool
for tool in tools:
    continue
    polya_vals = polya_data[tool]
    total_vals = total_data[tool]

    # Collect samples and sort numerically
    samples = sorted([k.split("_")[0] for k in polya_vals.keys()])
    data = []
    for s in samples:
        s_str = str(s)
        data.append({"Sample": s, "Type": "polyA", "Count": polya_vals.get(f"{s_str}_polyA", 0)})
        data.append({"Sample": s, "Type": "total", "Count": total_vals.get(f"{s_str}_total", 0)})

    df = pd.DataFrame(data)

    # Plot
    plt.figure(figsize=(10,6))
    ax = sns.barplot(
        data=df,
        x="Sample", y="Count", hue="Type",
        palette=palette,
        edgecolor="black"
    )

    ax.set_title(f"{tool} - PolyA vs Total")
    ax.set_ylabel("BSJ Sum")
    ax.set_xlabel("Sample")
    ax.set_yscale("log")
    plt.xticks()
    plt.tight_layout()
    
    # Save per tool
    plt.savefig(f"{tool}_bsj_amount_sum.png", dpi=300)
    plt.close()

tools = polya_data.keys()

#all_data = []
#for tool in tools:
#    polya_vals = polya_data[tool]
#    total_vals = total_data[tool]
#
#    samples = sorted([k.split("_")[0] for k in polya_vals.keys()])
#    for s in samples:
#        s_str = str(s)
#        all_data.append({
#            "Tool": tool,
#            "Type": "polyA",
#            "Count": polya_vals.get(f"{s_str}_polyA", 0)
#        })
#        all_data.append({
#            "Tool": tool,
#            "Type": "total",
#            "Count": total_vals.get(f"{s_str}_total", 0)
#        })
#
#df = pd.DataFrame(all_data)
#
#
#plt.figure(figsize=(10,6))
#ax = sns.boxplot(
#    data=df,
#    x="Tool", y="Count", hue="Type",
#    palette=palette
#)
#
#ax.set_title("PolyA vs Total across tools")
#ax.set_ylabel("BSJ Sum")
#ax.set_xlabel("Tool")
##ax.set_yscale("log")
#plt.xticks(rotation=45)
#plt.tight_layout()
#
#plt.savefig("all_tools_bsj_amount_boxplot.png", dpi=300)
#plt.close()


all_data = []
for tool in tools:
    polya_vals = polya_data[tool]
    total_vals = total_data[tool]

    # Extract per-sample values
    samples = sorted([k.split("_")[0] for k in polya_vals.keys()])
    for s in samples:
        s_str = str(s)
        all_data.append({
            "Tool": tool,
            "Type": "polyA",
            "Count": polya_vals.get(f"{s_str}_polyA", 0)
        })
        all_data.append({
            "Tool": tool,
            "Type": "total",
            "Count": total_vals.get(f"{s_str}_total", 0)
        })

df = pd.DataFrame(all_data)

name_map = { 'segemehl_filtered': 'Segemehl',
    'dcc_filtered': 'DCC',
    'ciriquant_filtered': 'CIRIquant',
    'circexplorer2_filtered': 'CIRCexplorer2',
    'find_circ_filtered':'find_circ'
    }

name_map_1 = { 'polyA': 'Poly(A)',
'total':'Total'
    }
df['Tool'] = df['Tool'].replace(name_map)
df['Type'] = df['Type'].replace(name_map_1)
tool_order = sorted(df["Tool"].unique())
plt.figure(figsize=(6,4))

cb_palette = sns.color_palette("colorblind")
palette = {"Total": cb_palette[0], "Poly(A)": cb_palette[1]}

ax = sns.boxplot(
    data=df,
    x="Tool", y="Count", hue="Type",
    palette=palette,
    order=tool_order,            # enforce alphabetical order
    showcaps=True,showfliers=False, 
    boxprops={'edgecolor': 'black'},       # ⬅ add black edge
    hue_order=["Total", "Poly(A)"],
    medianprops={'color': 'black'},        # optional: make median line black
    capprops={'color': 'black'},           # optional: black whisker caps
    whiskerprops={'color': 'black', 'linewidth': 1} 
)

sns.stripplot(
    data=df,
    x="Tool", y="Count", hue="Type",
    order=tool_order,
    hue_order=["Total", "Poly(A)"],
    dodge=True, jitter=True, alpha=0.9, size=4, color="black"
)

ax.set_title("BSJ Evidence per Tool in DEEP", fontsize=16)
ax.set_ylabel("Total BSJ Evidence")
ax.set_xlabel("Tool")
ax.yaxis.grid(True, linestyle='--', color='gray', alpha=0.7)
#ax.set_yscale("log")
# Fix double legend (boxplot + stripplot)
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], title="Type")
sns.despine()

plt.tight_layout()
plt.savefig("all_tools_bsj_amount_boxplot_jitter.png", dpi=300)
plt.close()
