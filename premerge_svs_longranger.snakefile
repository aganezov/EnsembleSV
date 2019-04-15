configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os


tools_methods = config["tools_methods"]
longranger_output_dir = config["data_output"]["longranger"]["dir"]
raw_sv_calls_dir = config["data_output"]["raw_sv_calls"]["dir"]
pre_merge_rel_process_dir = config["data_output"]["pre_merge_process"]["rel_dir"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]

def expected_longranger_result_premerge_sv_files():
	linked_read_bams = config["data_input"]["bams"]["linked"]
	linked_read_bases = [os.path.basename(name).split(".")[0] for name in linked_read_bams]
	min_size = config["data_premerge"]["min_len"]["linked"]
	return [os.path.join(pre_merge_aggregate_dir, base + "_longranger.m" + str(min_size) + ".rck.vcf") for base in linked_read_bases]

rule run_longranger_premerge_all:
	input:
		expected_longranger_result_premerge_sv_files()

rule copy_longranger_premerge_vcf_to_aggregate_dir:
	input: os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_longranger.m" + "{min_size}" + ".rck.vcf")
	output: os.path.join(pre_merge_aggregate_dir, "{base}" + "_longranger.m" + "{min_size}" + ".rck.vcf")
	message: "copying longranger premerge vcf {input} to the aggregate dir"
	shell:
	 "cp {input} {output}"

rule convert_premerge_rck_longranger_to_sniffles_format:
	input: os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_longranger.m" + "{min_size}" + ".rck.adj.tsv")
	output: os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_longranger.m" + "{min_size}" + ".rck.vcf")
	message: "converting longranger svs from RCK adjacency format after processing (file {input}) to Sniffles based VCF format file {output}"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_longranger.m" + wc.min_size,
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"

rule process_premerge_rck_longranger:
	input: os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_longranger.rck.adj.tsv")
	output: os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_longranger.m" + "{min_size}" + ".rck.adj.tsv")
	message: "processing longranger sv calls for min size"
	conda: tools_methods["rck"]["conda"]
	params: 
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		min_size=lambda wc: wc.min_size,
	shell:
		"{params.rck_adj_process} filter {input} --size-extra-field svlen --min-size {params.min_size} -o {output}"

rule concat_premerge_dels_and_large_rck_longranger:
	input:
		dels=os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_dels_longranger.rck.adj.tsv"),
		large=os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_largesvs_longranger.rck.adj.tsv")
	output: os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_longranger.rck.adj.tsv")
	message: "concatinating dels from {input.dels} and large svs from file {input.large} into {output}"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"]
	shell:
		"{params.rck_adj_process} cat {input.dels} {input.large} -o {output}"

rule convert_premerge_longranger_dels_to_rck:
	input: os.path.join(raw_sv_calls_dir, "{base}" + "_dels_longranger.vcf")
	output: os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_dels_longranger.rck.adj.tsv")
	message: "converting initial longranger dels vcf file {input} to RCK adjacency format file {output}"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_longranger_dels",
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
	shell:
		"{params.rck_adj_x2rck} longranger {input} --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} -o {output}"

rule convert_premerge_longranger_large_to_rck:
	input: os.path.join(raw_sv_calls_dir, "{base}" + "_largesvs_longranger.vcf")
	output: os.path.join(longranger_output_dir, pre_merge_rel_process_dir, "{base}" + "_largesvs_longranger.rck.adj.tsv")
	message: "converting initial longranger large svs vcf file {input} to RCK adjacency format file {output}"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_longranger_large",
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
	shell:
		"{params.rck_adj_x2rck} longranger {input} --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} -o {output}"



# include: "call_svs_longranger.snakefile"


