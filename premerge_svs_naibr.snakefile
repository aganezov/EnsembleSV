configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os


tools_methods = config["tools_methods"]
naibr_output_dir = config["data_output"]["naibr"]["dir"]
raw_sv_calls_dir = config["data_output"]["raw_sv_calls"]["dir"]
pre_merge_rel_process_dir = config["data_output"]["pre_merge_process"]["rel_dir"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]

def expected_naibr_result_premerge_sv_files():
	linked_read_bams = config["data_input"]["bams"]["linked"]
	linked_read_bases = [os.path.basename(name).split(".")[0] for name in linked_read_bams]
	min_size = config["data_premerge"]["min_len"]["linked"]
	return [os.path.join(pre_merge_aggregate_dir, base + "_naibr.m" + str(min_size) + ".rck.vcf") for base in linked_read_bases]

rule run_naibr_premerge_all:
	input:
		expected_naibr_result_premerge_sv_files()

rule copy_naibr_premerge_vcf_to_aggregate_dir:
	input: os.path.join(naibr_output_dir, pre_merge_rel_process_dir, "{base}" + "_naibr.m" + "{min_size}" + ".rck.vcf")
	output: os.path.join(pre_merge_aggregate_dir, "{base}" + "_naibr.m" + "{min_size}" + ".rck.vcf")
	message: "copying naibr premerge vcf {input} to the aggregate dir"
	shell:
	 "cp {input} {output}"

rule convert_premerge_rck_naibr_to_sniffles_format:
	input: os.path.join(naibr_output_dir, pre_merge_rel_process_dir, "{base}" + "_naibr.m" + "{min_size}" + ".rck.adj.tsv")
	output: os.path.join(naibr_output_dir, pre_merge_rel_process_dir, "{base}" + "_naibr.m" + "{min_size}" + ".rck.vcf")
	message: "converting naibr svs from RCK adjacency format after processing (file {input}) to Sniffles based VCF format file {output}"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_naibr.m" + wc.min_size,
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"

rule process_premerge_rck_naibr:
	input: os.path.join(naibr_output_dir, pre_merge_rel_process_dir, "{base}" + "_naibr.rck.adj.tsv")
	output: os.path.join(naibr_output_dir, pre_merge_rel_process_dir, "{base}" + "_naibr.m" + "{min_size}" + ".rck.adj.tsv")
	message: "processing naibr sv calls for min size"
	conda: tools_methods["rck"]["conda"]
	params: 
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		min_size=lambda wc: wc.min_size,
	shell:
		"{params.rck_adj_process} filter {input} --size-extra-field svlen --min-size {params.min_size} -o {output}"

rule convert_premerge_naibr_to_rck:
	input: os.path.join(raw_sv_calls_dir, "{base}" + "_naibr.bedpe")
	output: os.path.join(naibr_output_dir, pre_merge_rel_process_dir, "{base}" + "_naibr.rck.adj.tsv")
	message: "converting initial naibr vcf file {input} to RCK adjacency format file {output}"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_naibr",
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
	shell:
		"{params.rck_adj_x2rck} naibr {input} --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} -o {output}"



include: "call_svs_naibr.snakefile"


