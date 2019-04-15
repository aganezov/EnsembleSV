configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

tools_methods = config["tools_methods"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]
merged_dir = config["data_output"]["merge"]["dir"]

rule run_merge_all:
	input: os.path.join(merged_dir, config["data_sample_name"] + ".survivor.vcf")

rule run_merge_survivor:
	output: os.path.join(merged_dir, "{sample_name}.survivor.vcf")
	input: os.path.join(merged_dir, "{sample_name}.config.survivor")
	message: "running SURVIVOR with input {input}"
	conda: tools_methods["survivor"]["conda"]
	params:
		survivor=tools_methods["survivor"]["path"],
		max_distance=tools_methods["survivor"]["sv_max_distance"],
		min_caller_cnt=tools_methods["survivor"]["min_caller_cnt"],
		sv_type_consider=tools_methods["survivor"]["sv_type_consideration"],
		sv_strands_consider=tools_methods["survivor"]["sv_strands_consideration"],
		distance_estimate=tools_methods["survivor"]["distance_estimate"],
		min_sv_size=tools_methods["survivor"]["min_sv_size"]
	shell:
		"{params.survivor} merge {input} {params.max_distance} {params.min_caller_cnt} {params.sv_type_consider} {params.sv_strands_consider} {params.distance_estimate} {params.min_sv_size} {output}"

rule run_prepare_survivor_file_input:
	output: os.path.join(merged_dir, "{sample_name}.config.survivor")
	message: "putting all vcf files from the premerge aggregate folder for SURVIVOR merging into {output}"
	run:
		with open(output[0], "wt") as destination:
			for file_name in (fn for fn in os.listdir(pre_merge_aggregate_dir) if fn.endswith(".vcf")):
				full_path = os.path.join(pre_merge_aggregate_dir, file_name)
				print(full_path, file=destination)

