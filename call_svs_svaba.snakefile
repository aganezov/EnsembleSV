configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

tools_methods = config["tools_methods"]
raw_sv_calls_dir = config["data_output"].get("raw_sv_calls", {}).get("dir", "./raw")
svaba_output_dir = os.path.join(raw_sv_calls_dir, config["data_output"].get("svaba", {}).get("dir", "svaba"))

def expected_svaba_result_sv_files():
	illumina_read_bams = config["data_input"]["bams"]["illumina"]
	illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
	result = []
	for base in illumina_read_bases:
		result.append(os.path.join(raw_sv_calls_dir, base + "_indel_svaba.vcf") for base in illumina_read_bases)
		result.append(os.path.join(raw_sv_calls_dir, "raw", base + "_sv_svaba.vcf") for base in illumina_read_bases)
	return result

def svaba_bases():
	return [os.path.basename(name).split(".")[0] for name in config["data_input"]["bams"]["illumina"]]

def svaba_base_to_file_bam_path():
	result = {os.path.basename(name).split(".")[0]: name for name in config["data_input"]["bams"]["illumina"]} 
	return result


rule run_svaba_all:
	input:
		expected_svaba_result_sv_files()

rule copy_svaba_vcf_to_aggregate_dir:
	input: indel=os.path.join(svaba_output_dir, "{base}" + "_svaba.svaba.indel.vcf"),
			sv=os.path.join(svaba_output_dir, "{base}" + "_svaba.svaba.sv.vcf")
	output: 
		indel=os.path.join(raw_sv_calls_dir, "{base}_indel_svaba.vcf"),
		sv=os.path.join(raw_sv_calls_dir, "{base}_sv_svaba.vcf")
	message: "copying svaba sv calls files to the aggregate raw sv calls dir"
	shell:
		"cp {input.indel} {output.indel} && cp {input.sv} {output.sv}"

rule run_svaba_workflow:
	input: 
		bam=lambda wc: svaba_base_to_file_bam_path()[wc.base],
		ref=config["data_input"]["fasta"]["ref"]
	output: indel=os.path.join(svaba_output_dir, "{base}_svaba.svaba.indel.vcf"),
			sv=os.path.join(svaba_output_dir, "{base}_svaba.svaba.sv.vcf"),
	message: "running svaba for {input}"
	threads: tools_methods["svaba"].get("threads", 16)
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["svaba"]["conda"])
	params:
		svaba=tools_methods["svaba"].get("path", "svaba"),
		name=lambda wc: wc.base + "_svaba",
		svaba_output_dir=svaba_output_dir
	shell:
		"cd {params.svaba_output_dir} || {params.svaba} run -t {input.bam} -G {input.ref} -p {threads} -a {params.name} --germline"
