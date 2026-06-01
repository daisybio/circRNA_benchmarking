import os
import sys
import csv


DATA_DIR=sys.argv[1]
if DATA_DIR[-1] != "/":
    DATA_DIR = f"{DATA_DIR}/"

stats = []

for t in ["blacklist"]:
    root_dirs = [
        f"{DATA_DIR}polya/filtered_bed_{t}/circexplorer2_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/circexplorer2_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/segemehl_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/segemehl_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/circtools_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/circtools_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/ciri_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/ciri_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/find_circ_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/find_circ_filtered_{t}/",
        f"{DATA_DIR}polya/filtered_bed_{t}/circrna_finder_filtered_{t}/",
        f"{DATA_DIR}total/filtered_bed_{t}/circrna_finder_filtered_{t}/",
    ]
    cutoff = 5

    columns = ["filtered", "tool", "origin", "total", "pct", "remaining", "type"]

    for root_dir in root_dirs:
        origin = root_dir.split("/")[-4]
        for filename in os.listdir(root_dir):
            retained = 0
            skipped = 0
            total = 0
            if filename.endswith(".bed"):
                input_path = os.path.join(root_dir, filename)
                tool_name = root_dir.split("/")[-2]
                sample = filename.split("/")[-1].split(".")[0].split("_")[0]
                sample_name = filename.split("/")[-1].split(".")[0]
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
                                skipped += 1
                            else:
                                retained += 1
                                outfile.write(line)
                        except ValueError:
                            continue
            stats.append((origin, tool_name, sample_name, skipped, retained, total, retained/total))

        print(f"Stats for min evidence filtering for {sample_name}")
        print(f"Retained: {retained}, Removed: {skipped} Ratio: {retained/total}")


col_names = ["origin", "tool", "sample", "skipped", "retained", "total", "ratio"]

out_path = os.path.join(DATA_DIR, "evidence_filter_stats.tsv")
with open(out_path, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(col_names)
    writer.writerows(stats)