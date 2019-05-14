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
long_read_bams = config["data_input"]["bams"]["long"]
long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
long_methods_regex = "(" + "|".join(long_methods) + ")"
long_bases_regex = "(" + "|".join(long_read_bases) + ")"

# print([os.path.join(rck_dir, "{base}_" + method + ".sens.rck.adj.tsv") for method in long_methods])
def expected_long_spes():
	# print("test")
	long_read_bams = config["data_input"]["bams"]["long"]
	long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
	result = [os.path.join(aggreagate_merged_dir, base + ".spes.rck.vcf") for base in long_read_bases]	
	# print(result)
	return result

def expected_long_sens():
	long_read_bams = config["data_input"]["bams"]["long"]
	long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
	result = [os.path.join(aggreagate_merged_dir, base + ".sens.rck.vcf") for base in long_read_bases]	
	return result

def greater_even_regex_for_number(number):
	if number < 10:
		result = "^([{number}-9]|\\d{2}|\\d{3}|\\d{4})$".format(number=number)
	elif number < 100:
		tens = int(number / 10)
		deriv = number - int(number / 10) * 10
		if number % 10 == 0:
			result = "^([{tens}-9]\\d|\\d{3}|\\d{4})$".format(tens=tens)
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


rule all:
	input: expected_long_spes(),
		   expected_long_sens()


rule get_long_spes:
	input: os.path.join(rck_dir, "{base}.spes.rck.adj.tsv")
	output: os.path.join(aggreagate_merged_dir, "{base," + long_bases_regex + "}.spes.rck.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_spes",
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"


rule get_long_spes_rck:
	input: os.path.join(rck_dir, "{base}.sens.rck.adj.tsv")
	output: os.path.join(rck_dir, "{base," + long_bases_regex + "}.spes.rck.adj.tsv")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		re_regexes = lambda wc: regex_extra_re_string([wc.base + "_" + method + "_re" for method in long_methods], greater_even_regex_for_number(int(config["data_input"]["coverage"][wc.base] * config["data_merge"]["long_spec"]["min_support_fraction"]))),
		min_size=config["data_merge"]["long_spec"]["min_len"]
	shell:
		"{params.rck_adj_process} filter {input} {params.re_regexes} --min-size {params.min_size} --size-extra-field svlen -o {output}"

rule get_long_sens_rck_vcf:
	input: os.path.join(rck_dir, "{base}.sens.rck.adj.tsv")
	output: os.path.join(aggreagate_merged_dir, "{base," + long_bases_regex + "}.sens.rck.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_sens",
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"


rule get_long_sens_rck:
	input: survivor=os.path.join(merged_dir, "{base}.sens.survivor.vcf"),
		   rck_files=[os.path.join(rck_dir, "{base}_" + method + ".sens.rck.adj.tsv") for method in long_methods]
	output: os.path.join(rck_dir, "{base," + long_bases_regex +"}.sens.rck.adj.tsv")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		samples=lambda wc: ",".join(wc.base + "_" + method for method in long_methods),
		samples_source=lambda wc: ",".join(os.path.join(rck_dir, wc.base + "_" + method + ".sens.rck.adj.tsv") for method in long_methods),
		suffix=lambda wc: wc.base + "_sens",
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
	shell:
		"{params.rck_adj_x2rck} survivor {input.survivor} --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} --samples-suffix-extra --samples {params.samples} --samples-source {params.samples_source} --survivor-prefix {params.suffix} -o {output}"


rule get_long_sens_survivor:
	input: 
		   vcfs=[os.path.join(rck_dir, "{base}_" + method + ".sens.rck.vcf") for method in long_methods],
		   survivor_file=os.path.join(merged_dir, "{base}.sens.survivor")
	output: os.path.join(merged_dir, "{base," + long_bases_regex +"}.sens.survivor.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["survivor"]["conda"])
	params:
		survivor=tools_methods["survivor"]["path"],
		max_distance=config["data_premerge"]["survivor"]["max_distance"],
		min_caller_cnt=1,
		sv_type_consider=0,
		sv_strands_consider=1,
		distance_estimate=0,
		min_sv_size=30
	shell:
		"{params.survivor} merge {input.survivor_file} {params.max_distance} {params.min_caller_cnt} {params.sv_type_consider} {params.sv_strands_consider} {params.distance_estimate} {params.min_sv_size} {output}" 

rule get_long_sens_survivor_config:
	output: os.path.join(merged_dir, "{base," + long_bases_regex +"}.sens.survivor")
	input: vcf=[os.path.join(rck_dir, "{base}_" + method + ".sens.rck.vcf") for method in long_methods]
	run:
		with open(output[0], "wt") as dest:
			for file_name in input.vcf:
				print(file_name, file=dest)

rule get_long_sens_vcf_for_survivor:
	input: os.path.join(rck_dir, "{base}_{method}.sens.rck.adj.tsv")
	output: os.path.join(rck_dir, "{base}_{method," + long_methods_regex + "}.sens.rck.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_" + wc.method,
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"


rule get_long_sens_rck_for_survivor:
	input: os.path.join(rck_dir, "{base}_{method}.rck.adj.tsv")
	output: os.path.join(rck_dir, "{base}_{method," + long_methods_regex + "}.sens.rck.adj.tsv")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		min_size=config["data_premerge"]["min_len"]["long"],
	shell:
		"{params.rck_adj_process} filter {input} --size-extra-field svlen --min-size {params.min_size} -o {output}"


rule get_long_initial_rck:
	output: os.path.join(rck_dir, "{base}_{method," + long_methods_regex + "}.rck.adj.tsv")
	input: os.path.join(raw_sv_calls_dir, "{base}_{method}.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		method=lambda wc: wc.method,
		suffix=lambda wc: wc.base + "_" + wc.method,
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
		sample_string=lambda wc: "--sample " + wc.base.lower() if wc.method == "pbsv" else ""
	shell:
		"{params.rck_adj_x2rck} {params.method} {input} --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} {params.sample_string} -o {output}"

