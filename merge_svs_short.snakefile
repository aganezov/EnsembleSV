configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

tools_methods = config["tools_methods"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]
merged_dir = config["data_output"]["merge"]["dir"]
rck_dir = os.path.join(merged_dir, "rck")
aggreagate_merged_dir = os.path.join(merged_dir, "merged")
raw_sv_calls_dir = os.path.join(config["data_output"]["raw_sv_calls"]["dir"], "raw")
illumina_methods = [method for method in config["tools_enabled_methods"] if method in config["tools_read_type_to_method"]["illumina"]]
linked_methods = [method for method in config["tools_enabled_methods"] if method in config["tools_read_type_to_method"]["linked"]]
short_methods = [method for method in config["tools_enabled_methods"] if method in config["tools_read_type_to_method"]["illumina"] + config["tools_read_type_to_method"]["linked"]]
short_methods_regex = "(" + "|".join(short_methods) + ")"

def expected_linked():
	if "linked" not in config["data_input"]["bams"]:
		return []
	linked_read_bams = config["data_input"]["bams"]["linked"]
	linked_read_bases = [os.path.basename(name).split(".")[0] for name in linked_read_bams]
	result = [os.path.join(aggreagate_merged_dir, base + ".rck.vcf") for base in linked_read_bases]	
	return result

def expected_illumina():
	if illumina not in config["data_input"]["bams"]:
		return []
	illumina_read_bams = config["data_input"]["bams"]["illumina"]
	illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
	result = [os.path.join(aggreagate_merged_dir, base + ".rck.vcf") for base in illumina_read_bases]	
	return result

def expected_short():
	return []
	if "illumina" in config["data_input"]["bams"] or "linked" in config["data_input"]["bams"]:
		result = os.path.join(aggreagate_merged_dir, config["data_sample_name"] + "_short.sens.rck.vcf")
	return result

def short_rck():
	linked_read_bams = []
	linked_read_bases = []
	if "linked" in config["data_input"]["bams"]:
		linked_read_bams = config["data_input"]["bams"]["linked"]
		linked_read_bases = [os.path.basename(name).split(".")[0] for name in linked_read_bams]
	illumina_read_bams = []
	illumina_read_bases = []
	if "illumina" in config["data_input"]["bams"]:
		illumina_read_bams = config["data_input"]["bams"]["illumina"]
		illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
	result = []
	for method in illumina_methods:
		for base in illumina_read_bases:
			result.append(os.path.join(rck_dir, base + "_" + method + ".sens.rck.adj.tsv"))
	for method in linked_methods:
		for base in linked_read_bases:
			result.append(os.path.join(rck_dir, base + "_" + method + ".sens.rck.adj.tsv"))
	return result

def short_rck_vcf():
	linked_read_bams = []
	linked_read_bases = []
	if "linked" in config["data_input"]["bams"]:
		linked_read_bams = config["data_input"]["bams"]["linked"]
		linked_read_bases = [os.path.basename(name).split(".")[0] for name in linked_read_bams]
	illumina_read_bams = []
	illumina_read_bases = []
	if "illumina" in config["data_input"]["bams"]:
		illumina_read_bams = config["data_input"]["bams"]["illumina"]
		illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
	result = []
	for method in illumina_methods:
		for base in illumina_read_bases:
			result.append(os.path.join(rck_dir, base + "_" + method + ".sens.rck.vcf"))
	for method in linked_methods:
		for base in linked_read_bases:
			result.append(os.path.join(rck_dir, base + "_" + method + ".sens.rck.vcf"))
	return result

def survivor_samples():
	result = []
	if "linked" in config["data_input"]["bams"]:
		linked_read_bams = config["data_input"]["bams"]["linked"]
		linked_read_bases = [os.path.basename(name).split(".")[0] for name in linked_read_bams]
		for method in linked_methods:
			for base in linked_read_bases:
				result.append(base + "_" + method)
	if "illumina" in config["data_input"]["bams"]:
		illumina_read_bams = config["data_input"]["bams"]["illumina"]
		illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
		for method in illumina_methods:
			for base in illumina_read_bases:
				result.append(base + "_" + method)
	return result	



rule get_short_merged_all:
	input: expected_short()

rule get_short_sens_vcf:
	input: os.path.join(rck_dir, config["data_sample_name"] + "_short.sens.rck.adj.tsv")
	output: os.path.join(aggreagate_merged_dir, config["data_sample_name"] + "_short.sens.rck.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: config["data_sample_name"] + "_short_sens",
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"


rule get_short_sens_rck:
	input: survivor=os.path.join(merged_dir, config["data_sample_name"] + "_short.sens.survivor.vcf"),
		   rck_files=short_rck(),
	output: os.path.join(rck_dir, config["data_sample_name"] + "_short.sens.rck.adj.tsv")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		samples=lambda wc: ",".join(survivor_samples()),
		samples_source=lambda wc: ",".join(short_rck()),
		suffix=config["data_sample_name"] + "_short_sens",
		chr_include=lambda wc: ("--chrs-include " + ",".join(config["data_merge_spec"]["chr_include"]["regions"])) if ("chr_include" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_include"]) else "",
		chr_include_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_include"]["file"]) if ("chr_include" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_include"]) else "",
		chr_exclude=lambda wc: ("--chrs-exclude " + ",".join(config["data_merge_spec"]["chr_exclude"]["regions"])) if ("chr_exclude" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_exclude"]) else "",
		chr_exclude_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_exclude"]["file"]) if ("chr_exclude" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_exclude"]) else "",
	shell:
		"{params.rck_adj_x2rck} survivor {input.survivor} --id-suffix {params.suffix} {params.chr_include} {params.chr_include_file} {params.chr_exclude} {params.chr_exclude_file} --samples-suffix-extra --samples {params.samples} --samples-source {params.samples_source} --survivor-prefix {params.suffix} -o {output}"


rule get_short_sens_survivor:
	input: os.path.join(merged_dir, config["data_sample_name"] + "_short.sens.survivor")
	output: os.path.join(merged_dir, config["data_sample_name"] + "_short.sens.survivor.vcf"),
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["survivor"]["conda"])
	params:
		survivor=tools_methods["survivor"]["path"],
		max_distance=config["data_merge_sens"]["survivor"]["max_distance"],
		min_caller_cnt=config["data_merge_sens"]["min_cnt"]["short"],
		sv_type_consider=0,
		sv_strands_consider=1,
		distance_estimate=0,
		min_sv_size=lambda wc: min([config["data_merge_sens"]["min_len"]["illumina"], config["data_merge_sens"]["min_len"]["linked"]]),
	shell:
		"{params.survivor} merge {input} {params.max_distance} {params.min_caller_cnt} {params.sv_type_consider} {params.sv_strands_consider} {params.distance_estimate} {params.min_sv_size} {output}" 

rule get_short_survivor_config:
	output: os.path.join(merged_dir, config["data_sample_name"] + "_short.sens.survivor")
	input: short_rck_vcf()
	run:
		with open(output[0], "wt") as dest:
			for file_name in input:
				print(file_name, file=dest)

rule get_short_sens_vcf_for_survivor:
	input: os.path.join(rck_dir, "{base}_{method}.sens.rck.adj.tsv")
	output: os.path.join(rck_dir, "{base}_{method," + short_methods_regex + "}.sens.rck.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_" + wc.method,
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"


rule get_short_sens_rck_for_survivor:
	input: os.path.join(rck_dir, "{base}_{method}.rck.adj.tsv")
	output: os.path.join(rck_dir, "{base}_{method," + short_methods_regex +"}.sens.rck.adj.tsv")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		min_size=config["data_merge_sens"]["min_len"]["short"],
	shell:
		"{params.rck_adj_process} filter {input} --size-extra-field svlen --min-size {params.min_size} -o {output}"

rule get_short_initial_naibr:
	output: os.path.join(rck_dir, "{base}_naibr.rck.adj.tsv")
	input: os.path.join(raw_sv_calls_dir, "{base}_naibr.bedpe")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		method="naibr",
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_naibr",
		chr_include=lambda wc: ("--chrs-include " + ",".join(config["data_merge_spec"]["chr_include"]["regions"])) if ("chr_include" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_include"]) else "",
		chr_include_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_include"]["file"]) if ("chr_include" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_include"]) else "",
		chr_exclude=lambda wc: ("--chrs-exclude " + ",".join(config["data_merge_spec"]["chr_exclude"]["regions"])) if ("chr_exclude" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_exclude"]) else "",
		chr_exclude_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_exclude"]["file"]) if ("chr_exclude" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_exclude"]) else "",
	shell:
		"{params.rck_adj_x2rck} {params.method} {input} --id-suffix {params.suffix} {params.chr_include} {params.chr_include_file} {params.chr_exclude} {params.chr_exclude_file} -o {output}"


rule get_short_initial_rck_longranger:
	input: dels=os.path.join(rck_dir, "{base}_longranger_dels.rck.adj.tsv"),
		   large=os.path.join(rck_dir, "{base}_longranger_large_svs.rck.adj.tsv")
	output: dels=os.path.join(rck_dir, "{base}_longranger.rck.adj.tsv")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
	shell:
		"{params.rck_adj_process} cat {input.large} {input.dels} -o {output}"

rule get_short_initial_rck_longranger_generic:
	output: os.path.join(rck_dir, "{base}_longranger_{sv_type,(large_svs|dels)}.rck.adj.tsv")
	input: os.path.join(raw_sv_calls_dir, "{base}_longranger_{sv_type}.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		method="longranger",
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_longranger_" + wc.sv_type,
		chr_include=lambda wc: ("--chrs-include " + ",".join(config["data_merge_spec"]["chr_include"]["regions"])) if ("chr_include" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_include"]) else "",
		chr_include_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_include"]["file"]) if ("chr_include" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_include"]) else "",
		chr_exclude=lambda wc: ("--chrs-exclude " + ",".join(config["data_merge_spec"]["chr_exclude"]["regions"])) if ("chr_exclude" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_exclude"]) else "",
		chr_exclude_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_exclude"]["file"]) if ("chr_exclude" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_exclude"]) else "",
	shell:
		"{params.rck_adj_x2rck} {params.method} {input} --id-suffix {params.suffix} {params.chr_include} {params.chr_include_file} {params.chr_exclude} {params.chr_exclude_file} -o {output}"


rule get_short_initial_rck_svaba:
	output: os.path.join(rck_dir, "{base}_svaba.rck.adj.tsv")
	input:  indel=os.path.join(rck_dir, "{base}_svaba_indel.rck.adj.tsv"),
			sv=os.path.join(rck_dir, "{base}_svaba_sv.rck.adj.tsv")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
	shell:
		"{params.rck_adj_process} cat {input.sv} {input.indel} -o {output}"

rule get_short_initial_svaba_generic:
	output: os.path.join(rck_dir, "{base}_svaba_{itype,(sv|indel)}.rck.adj.tsv")
	input:  os.path.join(raw_sv_calls_dir, "{base}_{itype}_svaba.vcf") 
	conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_svaba_" + wc.itype,
		itype=lambda wc: wc.itype,
		chr_include=lambda wc: ("--chrs-include " + ",".join(config["data_merge_spec"]["chr_include"]["regions"])) if ("chr_include" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_include"]) else "",
		chr_include_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_include"]["file"]) if ("chr_include" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_include"]) else "",
		chr_exclude=lambda wc: ("--chrs-exclude " + ",".join(config["data_merge_spec"]["chr_exclude"]["regions"])) if ("chr_exclude" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_exclude"]) else "",
		chr_exclude_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_exclude"]["file"]) if ("chr_exclude" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_exclude"]) else "",
	shell:
		"{params.rck_adj_x2rck} svaba {input} --i-type {params.itype} --id-suffix {params.suffix} {params.chr_include} {params.chr_include_file} {params.chr_exclude} {params.chr_exclude_file} -o {output}"


rule get_short_initial_rck_generic:
	output: os.path.join(rck_dir, "{base}_{method," + short_methods_regex +"}.rck.adj.tsv")
	input: os.path.join(raw_sv_calls_dir, "{base}_{method}.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		method=lambda wc: wc.method,
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_" + wc.method,
		chr_include=lambda wc: ("--chrs-include " + ",".join(config["data_merge_spec"]["chr_include"]["regions"])) if ("chr_include" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_include"]) else "",
		chr_include_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_include"]["file"]) if ("chr_include" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_include"]) else "",
		chr_exclude=lambda wc: ("--chrs-exclude " + ",".join(config["data_merge_spec"]["chr_exclude"]["regions"])) if ("chr_exclude" in config["data_merge_spec"] and "regions" in config["data_merge_spec"]["chr_exclude"]) else "",
		chr_exclude_file=lambda wc: ("--chrs-include-file " + config["data_merge_spec"]["chr_exclude"]["file"]) if ("chr_exclude" in config["data_merge_spec"] and "file" in config["data_merge_spec"]["chr_exclude"]) else "",
		sample_string=lambda wc: "--samples Tumor" if wc.method == "grocsvs" else ""
	shell:
		"{params.rck_adj_x2rck} {params.method} {input} --id-suffix {params.suffix} {params.chr_include} {params.chr_include_file} {params.chr_exclude} {params.chr_exclude_file} -o {output}"

ruleorder: get_short_initial_rck_longranger_generic > get_short_initial_svaba_generic > get_short_initial_rck_generic

include: "call_svs.snakefile"