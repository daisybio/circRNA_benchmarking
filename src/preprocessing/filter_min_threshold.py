import os
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys


DATA_DIR=sys.argv[1]
if DATA_DIR[-1] != "/":
    DATA_DIR = f"{DATA_DIR}/"

for t in ["blacklist"]:
    root_dirs = [
        f"{DATA_DIR}polya/filtered_bed_{t}/circexplorer2_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/circexplorer2_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/segemehl_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/segemehl_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/dcc_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/dcc_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/ciriquant_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/ciriquant_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/find_circ_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/find_circ_filtered_{t}/",
    ]
    cutoff = 5

    print("filtered\ttool\torigin\ttotal\tpct\tremaining\ttype")
    columns = ["filtered", "tool", "origin", "total", "pct", "remaining", "type"]
    df = pd.DataFrame(columns=columns)

    for root_dir in root_dirs:
        i = 0
        total = 0
        origin = root_dir.split("/")[0]

        for filename in os.listdir(root_dir):
            if filename.endswith(".bed"):
                input_path = os.path.join(root_dir, filename)
                tool_name = root_dir.split("/")[-2]
                sample = filename.split("/")[-1].split(".")[0].split("_")[0]
                output_filename = os.path.join(root_dir.replace("bed", "bed_min"), filename)
                output_path = output_filename.replace(
                    f"{tool_name}_filtered", f"{tool_name}_filtered_min"
                )
                f_name = os.path.basename(output_path)
                out_dir = output_path.replace(f_name, "")
                os.makedirs(out_dir, exist_ok=True)

                if os.path.abspath(input_path) == os.path.abspath(output_path):
                    raise RuntimeError(
                        "Output path is same as input path — refusing to overwrite!"
                    )
                with open(input_path, "r") as infile, open(output_path, "w") as outfile:
                    for line in infile:
                        if line.strip() == "":
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) < 3:
                            continue
                        total += 1
                        try:
                            score = int(parts[4])
                            if score < cutoff:
                                i += 1
                            else:
                                outfile.write(line)
                        except ValueError:
                            continue

        tool_name = root_dir.split("/")[-2]
        print(f"{i}\t{tool_name}\t{origin}\t{total}\t{i/total}\t{total-i}\t{t}")
        df.loc[len(df)] = [i, tool_name, origin, total, i / total, total - i, t]