configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

include: "call_svs_sniffles.snakefile"

tools_methods = config["tools_methods"]
sniffles_output_dir = config["data_output"]["sniffles"]["dir"]
raw_sv_calls_dir = config["data_output"]["raw_sv_calls"]["dir"]
pre_merge_rel_process_dir = config["data_output"]["pre_merge_process"]["rel_dir"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]

def expected_sniffles_result_premerge_sv_files():
	long_read_bams = config["data_input"]["bams"]["long"]
	long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
	min_size = config["data_premerge"]["min_len"]["long"]
	min_support_by_base = {base: int(config["data_premerge"]["min_support_fraction"]["long"] * config["data_input"]["coverage"]["base_" + base]) for base in long_read_bases}
	return [os.path.join(pre_merge_aggregate_dir, base + "_sniffles.ms" + str(min_support_by_base[base]) + ".m" + str(min_size) + ".rck.vcf") for base in long_read_bases]

rule run_sniffles_premerge_all:
	input:
		expected_sniffles_result_premerge_sv_files()

rule copy_sniffles_premerge_vcf_to_aggregate_dir:
	input: os.path.join(sniffles_output_dir, pre_merge_rel_process_dir, "{base}" + "_sniffles.ms" + "{min_support,[^\.]}" + ".m" + "{min_size,[^\.]}" + ".rck.vcf")
	output: os.path.join(pre_merge_aggregate_dir, "{base}_sniffles.ms{min_support,[^\.]+}.m{min_size,[^\.]+}.rck.vcf")
	message: "copying sniffles premerge vcf {input} to the aggregate dir"
	shell:
	 "cp {input} {output}"

rule convert_premerge_rck_sniffles_to_sniffles_format:
	input: os.path.join(sniffles_output_dir, pre_merge_rel_process_dir, "{base}_sniffles.ms{min_support,[^\.]+}.m{min_size,[^\.]+}.rck.adj.tsv")
	output: os.path.join(sniffles_output_dir, pre_merge_rel_process_dir, "{base}_sniffles.ms{min_support,[^\.]+}.m{min_size,[^\.]+}.rck.vcf")
	message: "converting sniffles svs from RCK adjacency format after processing (file {input}) to Sniffles based VCF format file {output}"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_sniffles.ms" + wc.min_support + ".m" + wc.min_size,
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"

rule process_premerge_rck_sniffles:
	input: os.path.join(sniffles_output_dir, pre_merge_rel_process_dir, "{base}_sniffles.ms{min_support,[^\.]+}.rck.adj.tsv")
	output: os.path.join(sniffles_output_dir, pre_merge_rel_process_dir, "{base}_sniffles.ms{min_support,[^\.]+}.m{min_size,[^\.]+}.rck.adj.tsv")
	message: "processing pbsv sv calls for min size"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params: 
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		min_size=lambda wc: wc.min_size,
	shell:
		"{params.rck_adj_process} filter {input} --size-extra-field svlen --min-size {params.min_size} -o {output}"

rule convert_premerge_sniffles_to_rck:
	input: os.path.join(sniffles_output_dir, pre_merge_rel_process_dir, "{base}_sniffles.ms{min_support,[^\.]+}.vcf")
	output: os.path.join(sniffles_output_dir, pre_merge_rel_process_dir, "{base}_sniffles.ms{min_support,[^\.]}.rck.adj.tsv")
	message: "converting sniffles vcf file {input} to RCK adjacency format file {output}"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_sniffles.ms" + wc.min_support ,
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
	shell:
		"{params.rck_adj_x2rck} sniffles {input} --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} -o {output}"

rule filter_premerge_sniffles_by_min_support:
	input: os.path.join(raw_sv_calls_dir, "{base}_sniffles.vcf")
	output: os.path.join(sniffles_output_dir, pre_merge_rel_process_dir, "{base}_sniffles.ms{min_support,[^\.]+}.vcf")
	message: "filtering original sniffles sv calls to retain only those that are support but a minimum number of long reads spanning them"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["sniffles"]["support_cnt_script"]["conda"])
	params: 
		script=tools_methods["sniffles"]["support_cnt_script"]["path"],
		min_support=lambda wc: wc.min_support,
		python=tools_methods["sniffles"]["support_cnt_script"]["python"]
	shell: 
		"{params.python} {params.script} -i {input} -o {output} --min-support {params.min_support}"

