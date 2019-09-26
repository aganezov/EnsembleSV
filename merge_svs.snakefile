configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

tools_methods = config["tools_methods"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]
merged_dir = config["data_output"]["merge"]["dir"]
rck_dir = os.path.join(merged_dir, "rck")
aggreagate_merged_dir = os.path.join(merged_dir, "merged")
raw_sv_calls_dir = os.path.join(config["data_output"]["raw_sv_calls"]["dir"], "raw")
long_methods = [method for method in config["tools_enabled_methods"] if method in config["tools_read_type_to_method"]["long"]]
short_methods = [method for method in config["tools_enabled_methods"] if method in config["tools_read_type_to_method"]["illumina"] + config["tools_read_type_to_method"]["linked"]]

long_read_bams = config["data_input"]["bams"]["long"] if "long" in config["data_input"]["bams"] else []
long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams] if len(long_methods) > 0 and "long" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["long"]) > 0 else []

def merge_input_rck_vcf_files():
	result = []
	if len(short_methods) > 0 and (("illumina" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["illumina"]) > 0) or ("linked" in config["data_input"]["bams"] and  len(config["data_input"]["bams"]["linked"]) > 0)):
		result.append(os.path.join(aggreagate_merged_dir, config["data_sample_name"] + "_short.sens.rck.vcf"))
	if len(long_methods) > 0 and len(config["data_input"]["bams"]["long"]) > 0:
		long_read_bams = config["data_input"]["bams"]["long"]
		long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
		for base in long_read_bases:
			result.append(os.path.join(aggreagate_merged_dir, base + ".sens.rck.vcf"))
	return result

def merge_input_rck_files():
	result = []
	if len(short_methods) > 0 and (("illumina" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["illumina"]) > 0) or ("linked" in config["data_input"]["bams"] and  len(config["data_input"]["bams"]["linked"]) > 0)):
		result.append(os.path.join(rck_dir, config["data_sample_name"] + "_short.sens.rck.adj.tsv"))
	if len(long_methods) > 0 and "long" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["long"]) > 0:
		long_read_bams = config["data_input"]["bams"]["long"]
		long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
		for base in long_read_bases:
			result.append(os.path.join(rck_dir, base + ".sens.rck.adj.tsv"))
	return result

def stats_files():
	result = []
	if len(merge_input_rck_files()) > 0:
		result.append(os.path.join(aggreagate_merged_dir, config["data_sample_name"] + ".spes.main_stats.txt"))
	if len(short_methods) > 0 and (("illumina" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["illumina"]) > 0) or ("linked" in config["data_input"]["bams"] and  len(config["data_input"]["bams"]["linked"]) > 0)):
		result.append(os.path.join(aggreagate_merged_dir, config["data_sample_name"] + ".spes.short_stats.txt"))
	if len(long_methods) > 0 and "long" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["long"]) > 0:
		long_read_bams = config["data_input"]["bams"]["long"]
		long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
		for base in long_read_bases:
			result.append(os.path.join(aggreagate_merged_dir, config["data_sample_name"] + ".spes." + base + "_stats.txt"))
	return result

def call_set_files():
	present = False
	if len(merge_input_rck_files()) > 0:
		present = True
	if len(long_methods) > 0 and "long" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["long"]) > 0:
		present = True
	if present:
		result = [os.path.join(aggreagate_merged_dir, config["data_sample_name"] + ".spes.rck.vcf"),
		   		  os.path.join(rck_dir, config["data_sample_name"] + ".spes.rck.adj.tsv")]
	return result

def greater_even_regex_for_number(number):
	if number < 10:
		result = "^([{a}-9]|".format(a=number) + "\\d{2}|\\d{3}|\\d{4})$"
	elif number < 100:
		tens = int(number / 10)
		deriv = number - int(number / 10) * 10
		if number % 10 == 0:
			result = "^([{tens}-9]".format(tens=tens) + "\\d|\\d{3}|\\d{4})$"
		else:
			result = "^({tens}[{deriv}-9]|".format(tens=tens, deriv=deriv)
			if number < 90:
				result += "[{tensp}-9]\\d|".format(tensp=tens+1)
			result += "\\d{3}|\\d{4})$"
	else:
		raise Exception()
	return result

def regex_extra_re_string(method_specific_re_names, regex):
	result = ""
	for re_name in method_specific_re_names:
		result += " --keep-extra-field-regex \"" + re_name.lower() + "=" + regex  + "\" "
	return result

def long_samples():
	return [base + "_sens" for base in long_read_bases]
	

def short_sample():
	return config["data_sample_name"] + "_short_sens" if len(short_methods) > 0 and (("illumina" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["illumina"]) > 0) or ("linked" in config["data_input"]["bams"] and  len(config["data_input"]["bams"]["linked"]) > 0)) else ""

def samples():
	result = []
	if len(short_sample()) > 0:
		result += [short_sample()]
	return  result + long_samples()

def long_sample_sources():
	return [os.path.join(rck_dir, base + ".sens.rck.adj.tsv") for base in long_read_bases]

def short_sample_sources():
	return os.path.join(rck_dir, config["data_sample_name"] + "_short.sens.rck.adj.tsv") if len(short_methods) > 0 and (("illumina" in config["data_input"]["bams"] and len(config["data_input"]["bams"]["illumina"]) > 0) or ("linked" in config["data_input"]["bams"] and  len(config["data_input"]["bams"]["linked"]) > 0)) else ""

def sample_sources():
	result = []
	if len(short_sample_sources()) > 0:
		result += [short_sample_sources()]
	return result + long_sample_sources


rule get_merged_all:
	input: stats_files(),
		   call_set_files()

rule get_filtered_rck_vcf:
	input:  os.path.join(rck_dir, config["data_sample_name"] + ".spes.rck.adj.tsv")
	output: os.path.join(aggreagate_merged_dir, config["data_sample_name"] + ".spes.rck.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	log:	os.path.join(aggreagate_merged_dir, "log", config["data_sample_name"] + ".spes.rck.vcf.log")
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=config["data_sample_name"] + "_call_set",
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output} &> {log}"

rule get_filtered_call_set_main_stats:
	input:  os.path.join(rck_dir, config["data_sample_name"] + ".spes.rck.adj.tsv")
	output: os.path.join(aggreagate_merged_dir, config["data_sample_name"] + ".spes.main_stats.txt")
	conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	log:	os.path.join(aggreagate_merged_dir, "log", config["data_sample_name"] + ".spes.main_stats.txt.log")
	params:
		rck_adj_stats=tools_methods["rck"]["rck_adj_stats"]["path"],
		sources_field=lambda wc: (config["data_sample_name"] + "_sens_supporting_sources").lower(),
	shell:
		"{params.rck_adj_stats} survivor-stat {input} --sources-field {params.sources_field} -o {output} &> {log}"

rule get_filtered_call_set_short_stats:
	input:  os.path.join(rck_dir, config["data_sample_name"] + ".spes.rck.adj.tsv")
	output: os.path.join(aggreagate_merged_dir, config["data_sample_name"] + ".spes.short_stats.txt")
	conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	log:	os.path.join(aggreagate_merged_dir, "log", config["data_sample_name"] + ".spes.short_stats.txt.log")
	params:
		rck_adj_stats=tools_methods["rck"]["rck_adj_stats"]["path"],
		sources_field=lambda wc: (config["data_sample_name"] + "_short_sens_supporting_sources").lower(),
	shell:
		"{params.rck_adj_stats} survivor-stat {input} --sources-field {params.sources_field} -o {output} &> {log}"

rule get_filtered_call_set_long_stats:
	input:  os.path.join(rck_dir, config["data_sample_name"] + ".spes.rck.adj.tsv")
	output: os.path.join(aggreagate_merged_dir, config["data_sample_name"] + ".spes.{long_base}_stats.txt")
	conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	log: os.path.join(aggreagate_merged_dir, "log", config["data_sample_name"] + ".spes.{long_base}_stats.txt.log")
	params:
		rck_adj_stats=tools_methods["rck"]["rck_adj_stats"]["path"],
		sources_field=lambda wc: (wc.long_base + "_sens_supporting_sources").lower(),
	shell:
		"{params.rck_adj_stats} survivor-stat {input} --sources-field {params.sources_field} -o {output} &> {log}"

rule get_filtered_call_set:
	output: os.path.join(rck_dir, config["data_sample_name"] + ".spes.rck.adj.tsv")
	input:  os.path.join(rck_dir, config["data_sample_name"] + ".sens.rck.adj.tsv")
	conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	log: 	os.path.join(rck_dir, "log", config["data_sample_name"] + ".spes.rck.adj.tsv.log")
	params:
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		re_regexes=lambda wc: " ".join(regex_extra_re_string([base + "_" + method + "_re" for method in long_methods], 
									          greater_even_regex_for_number(int(config["data_input"]["coverage"][base] * config["data_merge_spec"]["long_spec"]["min_support_fraction"])))
						   for base in long_read_bases) + " --keep-extra-field-regex \"" + config["data_sample_name"].lower() + "_sens_supporting_sources=" + config["data_sample_name"] + "_short_sens\"",
		min_size=config["data_merge_spec"]["min_len"]
	shell:
		"{params.rck_adj_process} filter {input} {params.re_regexes} --min-size {params.min_size} --size-extra-field svlen -o {output} &> {log}"


rule get_merged_sens_call_set_rck:
	output: os.path.join(rck_dir, config["data_sample_name"] + ".sens.rck.adj.tsv")
	input:  survivor_vcf=os.path.join(merged_dir, config["data_sample_name"] + ".sens.survivor.vcf"),
			rck_files=merge_input_rck_files(),
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	log: os.path.join(rck_dir, "log", config["data_sample_name"] + ".sens.rck.adj.tsv.log")
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		samples=lambda wc: ",".join(samples()),
		samples_source=lambda wc: ",".join(merge_input_rck_files()),
		suffix=lambda wc: config["data_sample_name"] + "_sens",
		chr_include=lambda wc: ("--chrs-include " + ",".join(config["data_merge_sens"]["chr_include"]["regions"])) if ("chr_include" in config["data_merge_sens"] and "regions" in config["data_merge_sens"]["chr_include"]) else "",
		chr_include_file=lambda wc: ("--chrs-include-file " + config["data_merge_sens"]["chr_include"]["file"]) if ("chr_include" in config["data_merge_sens"] and "file" in config["data_merge_sens"]["chr_include"]) else "",
		chr_exclude=lambda wc: ("--chrs-exclude " + ",".join(config["data_merge_sens"]["chr_exclude"]["regions"])) if ("chr_exclude" in config["data_merge_sens"] and "regions" in config["data_merge_sens"]["chr_exclude"]) else "",
		chr_exclude_file=lambda wc: ("--chrs-include-file " + config["data_merge_sens"]["chr_exclude"]["file"]) if ("chr_exclude" in config["data_merge_sens"] and "file" in config["data_merge_sens"]["chr_exclude"]) else "",
	shell:
		"{params.rck_adj_x2rck} survivor {input.survivor_vcf} --id-suffix {params.suffix} {params.chr_include} {params.chr_include_file} {params.chr_exclude} {params.chr_exclude_file}  --samples {params.samples} --samples-source {params.samples_source} --survivor-prefix {params.suffix} -o {output} &> {log}"

rule get_merged_sens_call_set_survivour:
	output: os.path.join(merged_dir, config["data_sample_name"] + ".sens.survivor.vcf")
	input:  os.path.join(merged_dir, config["data_sample_name"] + ".sens.survivor")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["survivor"]["conda"])
	log: os.path.join(merged_dir, "log", config["data_sample_name"] + ".sens.survivor.vcf.log")
	params:
		survivor=tools_methods["survivor"]["path"],
		max_distance=config["data_merge_sens"]["survivor"]["max_distance"],
		min_caller_cnt=1,
		sv_type_consider=0,
		sv_strands_consider=1,
		distance_estimate=0,
		min_sv_size=lambda wc: min([config["data_merge_sens"]["min_len"][key] for key in ["long", "illumina", "linked"]])
	shell:
		"{params.survivor} merge {input} {params.max_distance} {params.min_caller_cnt} {params.sv_type_consider} {params.sv_strands_consider} {params.distance_estimate} {params.min_sv_size} {output} &> {log}"


rule get_merged_sens_call_set_survivour_config:
	output: os.path.join(merged_dir, config["data_sample_name"] + ".sens.survivor")
	input:  merge_input_rck_vcf_files()
	run:
		with open(output[0], "wt") as dest:
			for file_name in input:
				print(file_name, file=dest)

include: "merge_svs_long.snakefile"
include: "merge_svs_short.snakefile"
include: "call_svs.snakefile"