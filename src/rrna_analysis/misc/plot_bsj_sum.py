import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import os

# Custom palette
palette = {
    "total": "#FC8D62",  # orange
    "polyA": "#66C2A5",  # green
}

MAIN_DATA_DIR = sys.argv[1]

# Input files
polya_file = os.path.join(MAIN_DATA_DIR, "polya.bsj_amount.json")
total_file = os.path.join(MAIN_DATA_DIR, "total.bsj_amount.json")

# Load JSONs
with open(polya_file) as f:
    polya_data = json.load(f)

with open(total_file) as f:
    total_data = json.load(f)

tools = polya_data.keys()

tools = polya_data.keys()



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
            "Count": polya_vals.get(f"{s_str}_polya", 0)
        })
        all_data.append({
            "Tool": tool,
            "Type": "total",
            "Count": total_vals.get(f"{s_str}_total", 0)
        })

df = pd.DataFrame(all_data)


name_map = { 'segemehl': 'Segemehl',
    'dcc': 'DCC',
    'ciriquant': 'CIRIquant',
    'circexplorer2': 'CIRCexplorer2',
    'find_circ':'find_circ'
    }

name_map_1 = { 'polyA': 'Poly(A)',
'total':'Total'
    }
df['Tool'] = df['Tool'].replace(name_map)
df['Type'] = df['Type'].replace(name_map_1)
tool_order = sorted(df["Tool"].unique())
plt.figure(figsize=(6, 4))
    
cb_palette = sns.color_palette("colorblind")
palette = {"Total": cb_palette[0], "Poly(A)": cb_palette[1]}

ax = sns.boxplot(
    data=df,
    x="Tool", y="Count", hue="Type",
    palette=palette,
    order=tool_order,            # enforce alphabetical order
    hue_order=["Total", "Poly(A)"],
    showcaps=True,showfliers=False, 
    boxprops={'edgecolor': 'black'},       # ⬅ add black edge
    medianprops={'color': 'black'},        # optional: make median line black
    capprops={'color': 'black'},           # optional: black whisker caps
    whiskerprops={'color': 'black', 'linewidth': 1} 
)
ax.yaxis.grid(True, linestyle='--', color='gray', alpha=0.7)
sns.stripplot(
    data=df,
    x="Tool", y="Count", hue="Type",
    order=tool_order,            # enforce alphabetical order
    hue_order=["Total", "Poly(A)"],
    dodge=True, jitter=True, alpha=0.9, size=4, color="black"
)

ax.set_title("BSJ Evidence per Tool in GSE138734", fontsize=16)
ax.set_ylabel("Total BSJ Evidence")
ax.set_xlabel("Tool")
#ax.set_yscale("log")
#plt.xticks(rotation=45)
# Fix double legend (boxplot + stripplot)
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], title="Type")

sns.despine()
plt.tight_layout()
out = os.path.join(MAIN_DATA_DIR, "all_tools_bsj_amount_boxplot_jitter.png")
plt.savefig(out, dpi=300)
plt.close()
