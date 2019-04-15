configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

include: "call_svs_grocsvs.snakefile"


tools_methods = config["tools_methods"]
grocsvs_output_dir = config["data_output"]["grocsvs"]["dir"]
raw_sv_calls_dir = config["data_output"]["raw_sv_calls"]["dir"]
pre_merge_rel_process_dir = config["data_output"]["pre_merge_process"]["rel_dir"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]

def expected_grocsvs_result_premerge_sv_files():
	illumina_read_bams = config["data_input"]["bams"]["illumina"]
	illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
	min_size = config["data_premerge"]["min_len"]["illumina"]
	return [os.path.join(pre_merge_aggregate_dir, base + "_grocsvs.m" + str(min_size) + "rck.vcf") for base in illumina_read_bases]

rule run_grocsvs_premerge_all:
	input:
		expected_grocsvs_result_premerge_sv_files()

rule copy_grocsvs_premerge_vcf_to_aggregate_dir:
	input: os.path.join(grocsvs_output_dir, pre_merge_rel_process_dir, "{base}" + "_grocsvs.m" + "{min_size}" + "rck.vcf")
	output: os.path.join(pre_merge_aggregate_dir, "{base}" + "_grocsvs.m" + "{min_size}" + "rck.vcf")
	message: "copying grocsvs premerge vcf {input} to the aggregate dir"
	shell:
	 "cp {input} {output}"

rule convert_premerge_rck_grocsvs_to_sniffles_format:
	input: os.path.join(grocsvs_output_dir, pre_merge_rel_process_dir, "{base}" + "_grocsvs.m" + "{min_size}" + ".rck.adj.tsv")
	output: os.path.join(grocsvs_output_dir, pre_merge_rel_process_dir, "{base}" + "_grocsvs.m" + "{min_size}" + "rck.vcf")
	message: "converting grocsvs svs from RCK adjacency format after processing (file {input}) to Sniffles based VCF format file {output}"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_grocsvs.m" + wc.min_size,
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"

rule process_premerge_rck_grocsvs:
	input: os.path.join(grocsvs_output_dir, pre_merge_rel_process_dir, "{base}" + "_grocsvs.rck.adj.tsv")
	output: os.path.join(grocsvs_output_dir, pre_merge_rel_process_dir, "{base}" + "_grocsvs.m" + "{min_size}" + ".rck.adj.tsv")
	message: "processing grocsvs sv calls for min size"
	conda: tools_methods["rck"]["conda"]
	params: 
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		min_size=lambda wc: wc.min_size,
	shell:
		"{params.rck_adj_process} filter {input} --min-size {params.min_size} -o {output}"

rule convert_premerge_grocsvs_to_rck:
	input: os.path.join(raw_sv_calls_dir, "{base}" + "_grocsvs.vcf")
	output: os.path.join(grocsvs_output_dir, pre_merge_rel_process_dir, "{base}" + "_grocsvs.rck.adj.tsv")
	message: "converting initial grocsvs vcf file {input} to RCK adjacency format file {output}"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_grocsvs",
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
	shell:
		"{params.rck_adj_x2rck} grocsv {input} --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} -o {output}"


