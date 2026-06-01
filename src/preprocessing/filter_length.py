import os
import sys
import csv


DATA_DIR = sys.argv[1]
if DATA_DIR[-1] != "/":
    DATA_DIR = f"{DATA_DIR}/"

CUTOFF = 100000

print(f"Filtering bsj by length cutoff {CUTOFF}bp")
for t in ["blacklist"]:
    root_dirs = [
        f"{DATA_DIR}polya/bed/circexplorer2/{t}",
        f"{DATA_DIR}total/bed/circexplorer2/{t}",
        f"{DATA_DIR}polya/bed/segemehl/{t}",
        f"{DATA_DIR}total/bed/segemehl/{t}",
        f"{DATA_DIR}polya/bed/circtools/{t}",
        f"{DATA_DIR}total/bed/circtools/{t}",
        f"{DATA_DIR}polya/bed/ciri/{t}",
        f"{DATA_DIR}total/bed/ciri/{t}",
        f"{DATA_DIR}polya/bed/find_circ/{t}",
        f"{DATA_DIR}total/bed/find_circ/{t}",
        f"{DATA_DIR}polya/bed/circrna_finder/{t}",
        f"{DATA_DIR}total/bed/circrna_finder/{t}"
    ]

    stats = []

    for root_dir in root_dirs:
        origin = root_dir.split("/")[-4]

        for filename in os.listdir(root_dir):
            skipped = 0
            retained = 0
            total = 0
            if filename.endswith(".bed"):
                input_path = os.path.join(root_dir, filename)
                tool_name = root_dir.split("/")[-2]
                sample = filename.split("/")[-1].split(".")[0].split("_")[0]
                sample_name = filename.split("/")[-1].split(".")[0]

                output_filename = os.path.join(root_dir.replace("blacklist", "").replace("bed", f"filtered_bed_{t}"), filename)
                output_path = output_filename.replace(tool_name, f"{tool_name}_filtered_{t}")
                f_name = os.path.basename(output_path)
                out_dir = output_path.replace(f_name, "")
                os.makedirs(out_dir, exist_ok=True)

                if os.path.abspath(input_path) == os.path.abspath(output_path):
                    raise RuntimeError("Output path is same as input path — refusing to overwrite!")
                with open(input_path, "r") as infile, open(output_path, "w") as outfile:
                    for line in infile:
                        if line.strip() == "":
                            continue
                        parts = line.strip().split("\t")
                        if len(parts) < 3:
                            continue
                        total+=1
                        try:
                            start = int(parts[1])
                            end = int(parts[2])
                            length = end - start
                            if length < CUTOFF:
                                outfile.write(line)
                                retained += 1
                            else:
                                skipped += 1
                        except ValueError:
                            continue

            ratio = retained/total
            print(f"Stats for length filtering on {filename}")
            print(f"Retained BSJ: {retained}, Removed BSJ: {skipped}, Ratio (%retained): {ratio}")
            stats.append((origin, tool_name, sample_name, skipped, retained, total, ratio))
    
    
col_names = ["origin", "tool", "sample", "skipped", "retained", "total", "ratio"]

out_path = os.path.join(DATA_DIR, "length_filter_stats.tsv")
with open(out_path, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(col_names)
    writer.writerows(stats)