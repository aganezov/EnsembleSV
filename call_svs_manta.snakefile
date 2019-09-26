configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

tools_methods = config["tools_methods"]
raw_sv_calls_dir = config["data_output"].get("raw_sv_calls", {}).get("dir", "./raw")
manta_output_dir = os.path.join(raw_sv_calls_dir, config["data_output"].get("manta", {}).get("dir", "manta"))

def expected_manta_result_sv_files():
	if "illumina" not in config["data_input"]["bams"]:
		return []
	illumina_read_bams = config["data_input"]["bams"]["illumina"]
	illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
	return [os.path.join(raw_sv_calls_dir, "raw", base + "_manta.vcf") for base in illumina_read_bases]

def manta_bases():
	return [os.path.basename(name).split(".")[0] for name in config["data_input"]["bams"]["illumina"]]

def manta_base_to_file_bam_path():
	result = {os.path.basename(name).split(".")[0]: name for name in config["data_input"]["bams"]["illumina"]} 
	return result


rule run_manta_all:
	input:
		expected_manta_result_sv_files()

rule copy_manta_vcf_to_aggregate_dir:
	input: os.path.join(manta_output_dir, "{base}" + "_manta.vcf")
	output: os.path.join(raw_sv_calls_dir, "raw", "{base}" + "_manta.vcf")
	message: "copying manta sv calls file {input} to the aggregate raw sv calls dir"
	shell:
		"cp {input} {output}"

rule run_manta_surface_result:
	input: lambda wc: os.path.join(manta_output_dir, wc.base + "_rundir", "results", "variants", "tumorSV.vcf.gz" if config["data_sample_type"] == "T" else "diploidSV.vcf.gz")
	output: os.path.join(manta_output_dir, "{base}_manta.vcf")
	params:
		vcf_file = lambda wc: os.path.join(manta_output_dir, wc.base +"_rundir", "results", "variants", "tumorSV.vcf.gz" if config["data_sample_type"] == "T" else "diploidSV.vcf.gz")
	shell:
		"gunzip -c {params.vcf_file} > {output}"

rule run_manta_workflow_normal:
	input: 
		run_scipt=os.path.join(manta_output_dir, "{base}_rundir", "runWorkflow.py"),
	output: os.path.join(manta_output_dir, "{base}_rundir", "results", "variants", "diploidSV.vcf.gz")
	message: "running manta workflow for {input}"
	threads: config["tools_methods"]["manta"].get("threads", 15)
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["manta"]["conda"])
	log: os.path.join(manta_output_dir, "{base}_rundir", "log", "diploidSV.vcf.gz.log")
	params:
		mode=tools_methods["manta"].get("mode", "local")
	shell:
		"{input.run_scipt} -j {threads} -m {params.mode}"

rule run_manta_workflow_tumor:
	input: 
		run_scipt=os.path.join(manta_output_dir, "{base}_rundir", "runWorkflow.py"),
	output: os.path.join(manta_output_dir, "{base}_rundir", "results", "variants", "tumorSV.vcf.gz")
	message: "running manta workflow for {input}"
	threads: config["tools_methods"]["manta"].get("threads", 15)
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["manta"]["conda"])
	log: os.path.join(manta_output_dir, "{base}_rundir", "log", "tumorSV.vcf.gz.log")
	params:
		mode=tools_methods["manta"].get("mode", "local")
	shell:
		"{input.run_scipt} -j {threads} -m {params.mode} &> {log}"


rule configure_manta:
	input: 
		bam=lambda wc: manta_base_to_file_bam_path()[wc.base],
		ref=config["data_input"]["fasta"]["ref"]
	output: os.path.join(manta_output_dir, "{base}_rundir", "runWorkflow.py")
	message: "running manta configuration for {input}"
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["manta"]["conda"])
	log: os.path.join(manta_output_dir, "{base}_rundir", "log", "runWorkflow.py.log")
	params:
		config_script=tools_methods["manta"].get("config_script", "configManta.py"),
		bam_flag=lambda wc: "--tumorBam" if config["data_sample_type"] == "T" else "--bam",
		prefix=lambda wc: os.path.join(manta_output_dir, wc.base)
	shell:
		"{params.config_script} {params.bam_flag} {input.bam} --referenceFasta {input.ref} --runDir {params.prefix}_rundir &> {log}"
