import os
import json
import sys


DATA_DIR=sys.argv[1]
if DATA_DIR[-1] != "/":
    DATA_DIR = f"{DATA_DIR}/"


print(f"Filtered,Tool,Total,Pct,Blacklist")
for t in ["blacklist"]:
    root_dirs = [
        f"{DATA_DIR}polya/bed/circexplorer2/{t}",
        f"{DATA_DIR}total/bed/circexplorer2/{t}",
        f"{DATA_DIR}polya/bed/segemehl/{t}",
        f"{DATA_DIR}total/bed/segemehl/{t}",
        f"{DATA_DIR}polya/bed/dcc/{t}",
        f"{DATA_DIR}total/bed/dcc/{t}",
        f"{DATA_DIR}polya/bed/ciriquant/{t}",
        f"{DATA_DIR}total/bed/ciriquant/{t}",
        f"{DATA_DIR}polya/bed/find_circ/{t}",
        f"{DATA_DIR}total/bed/find_circ/{t}"
    ]
    cutoff = 100000

    length_dict = dict()

    for root_dir in root_dirs:
        i = 0
        total = 0
        origin = root_dir.split("/")[0]
        if origin not in length_dict.keys():
            length_dict[origin] = dict()

        for filename in os.listdir(root_dir):
            if filename.endswith(".bed"):
                input_path = os.path.join(root_dir, filename)
                tool_name = root_dir.split("/")[-2]
                sample = filename.split("/")[-1].split(".")[0].split("_")[0]
                if tool_name not in length_dict[origin].keys():
                    length_dict[origin][tool_name] = dict()
                
                if sample not in length_dict[origin][tool_name].keys():
                    length_dict[origin][tool_name][sample] = []

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
                            length_dict[origin][tool_name][sample].append(length)
                            if length < cutoff:
                                outfile.write(line)
                            else:
                                # print(f"Filtering out {start}\t{end}")
                                i += 1
                        except ValueError:
                            continue

                # if os.path.abspath(input_path) == os.path.abspath(output_path):
                #         raise RuntimeError("Output path is same as input path — refusing to overwrite!")
                # with open(input_path, "r") as infile:
                #     for line in infile:
                #         if line.strip() == "":
                #             continue
                #         parts = line.strip().split("\t")
                #         if len(parts) < 3:
                #             continue
                #         try:
                #             start = int(parts[1])
                #             end = int(parts[2])
                #             length = end - start
                #             if length < cutoff:
                #                 continue
                #             else:
                #                 #print(f"Filtering out {start}\t{end}\tlength:{end-start}")
                #                 i += 1
                #         except ValueError:
                #             continue

        print(f"{i},{tool_name},{total},{i/total},{t}")