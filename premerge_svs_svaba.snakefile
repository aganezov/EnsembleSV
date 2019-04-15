configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os


tools_methods = config["tools_methods"]
svaba_output_dir = config["data_output"]["svaba"]["dir"]
raw_sv_calls_dir = config["data_output"]["raw_sv_calls"]["dir"]
pre_merge_rel_process_dir = config["data_output"]["pre_merge_process"]["rel_dir"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]

def expected_svaba_result_premerge_sv_files():
	illumina_read_bams = config["data_input"]["bams"]["illumina"]
	illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
	min_size = config["data_premerge"]["min_len"]["illumina"]
	return [os.path.join(pre_merge_aggregate_dir, base + "_svaba.m" + str(min_size) + ".rck.vcf") for base in illumina_read_bases]

rule run_svaba_premerge_all:
	input:
		expected_svaba_result_premerge_sv_files()

rule copy_svaba_premerge_vcf_to_aggregate_dir:
	input: os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}_svaba.m{min_size}.rck.vcf")
	output: os.path.join(pre_merge_aggregate_dir, "{base}_svaba.m{min_size,\d+}.rck.vcf")
	message: "copying svaba premerge vcf {input} to the aggregate dir"
	shell:
	 "cp {input} {output}"

rule convert_premerge_rck_svaba_to_sniffles_format:
	input: os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}" + "_svaba.m" + "{min_size}" + ".rck.adj.tsv")
	output: os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}_svaba.m{min_size,\d+}" + ".rck.vcf")
	message: "converting svaba svs from RCK adjacency format after processing (file {input}) to Sniffles based VCF format file {output}"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.base + "_svaba.m" + wc.min_size,
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output}"

rule process_premerge_rck_svaba:
	input: os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}" + "_svaba.rck.adj.tsv")
	output: os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}_svaba.m{min_size,\d+}.rck.adj.tsv")
	message: "processing svaba sv calls for min size"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params: 
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"],
		min_size=lambda wc: wc.min_size,
	shell:
		"{params.rck_adj_process} filter {input} --size-extra-field span --min-size {params.min_size} -o {output}"

rule concat_premerge_dels_and_large_rck_svaba:
	input:
		dels=os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}" + "_indels_svaba.rck.adj.tsv"),
		large=os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}" + "_sv_svaba.rck.adj.tsv")
	output: os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}_svaba.rck.adj.tsv")
	message: "concatinating indels from {input.dels} and large svs from file {input.large} into {output}"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_process=tools_methods["rck"]["rck_adj_process"]["path"]
	shell:
		"{params.rck_adj_process} cat {input.dels} {input.large} -o {output}"

rule convert_premerge_svaba_dels_to_rck:
	input: os.path.join(raw_sv_calls_dir, "{base}_indel_svaba.vcf")
	output: os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}_indels_svaba.rck.adj.tsv")
	message: "converting initial svaba dels vcf file {input} to RCK adjacency format file {output}"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_svaba_indels",
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
	shell:
		"{params.rck_adj_x2rck} svaba {input} --i-type indel --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} -o {output}"

rule convert_premerge_svaba_large_to_rck:
	input: os.path.join(raw_sv_calls_dir, "{base}_sv_svaba.vcf")
	output: os.path.join(svaba_output_dir, pre_merge_rel_process_dir, "{base}_sv_svaba.rck.adj.tsv")
	message: "converting initial svaba large svs vcf file {input} to RCK adjacency format file {output}"
	conda: os.path.join(tools_methods["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
	params:
		rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
		suffix=lambda wc: wc.base + "_svaba_sv",
		chr_include_file=config["data_premerge"]["chr_include"]["file"],
		chr_exclude=lambda wc: ",".join(config["data_premerge"]["chr_exclude"]["regions"]),
	shell:
		"{params.rck_adj_x2rck} svaba {input} --i-type sv --id-suffix {params.suffix} --chrs-include-file {params.chr_include_file} --chrs-exclude {params.chr_exclude} -o {output}"



# include: "call_svs_svaba.snakefile"


