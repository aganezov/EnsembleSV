
configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

tools_methods = config["tools_methods"]
raw_sv_calls_dir = config["data_output"].get("raw_sv_calls", {}).get("dir", "./raw")
lumpy_output_dir = os.path.join(raw_sv_calls_dir, config["data_output"].get("lumpy", {}).get("dir", "lumpy"))


def expected_lumpy_result_sv_files():
	if "illumina" not in config["data_input"]["bams"]:
		return []
	illumina_read_bams = config["data_input"]["bams"]["illumina"]
	illumina_read_bases = [os.path.basename(name).split(".")[0] for name in illumina_read_bams]
	return [os.path.join(raw_sv_calls_dir, "raw", base + "_lumpy.vcf") for base in illumina_read_bases]

def lumpy_bases():
	return [os.path.basename(name).split(".")[0] for name in config["data_input"]["bams"]["illumina"]]

def lumpy_base_to_file_bam_path():
	return {os.path.basename(name).split(".")[0]: name for name in config["data_input"]["bams"]["illumina"]}

rule run_lumpy_all:
	input:
		expected_lumpy_result_sv_files()

rule copy_lumpy_vcf_to_aggregate_dir:
	input: os.path.join(lumpy_output_dir, "{base}" + "_lumpy.vcf")
	output: os.path.join(raw_sv_calls_dir, "raw", "{base}" + "_lumpy.vcf")
	message: "copying lumpy sv calls file {input} to the aggregate raw sv calls dir"
	shell:
		"cp {input} {output}"

rule run_lumpy_call:
	input:  splitters=os.path.join(lumpy_output_dir, "{base}" + ".splitters.sort.bam"),
			discordants=os.path.join(lumpy_output_dir, "{base}" + ".discordants.sort.bam"),
			bam=lambda wc: lumpy_base_to_file_bam_path()[wc.base]
	output: os.path.join(lumpy_output_dir, "{base}" + "_lumpy.vcf")
	message: "running lumpy on bam {input.bam} ;splitters {input.splitters} discordants {input.discordants}; result calls in {output}"
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["lumpy"]["conda"])
	log: os.path.join(lumpy_output_dir, "log", "{base}" + "_lumpy.vcf.log")
	params:
		lumpy=tools_methods["lumpy"].get("path", "lumpyexpress")
	shell: 
		"{params.lumpy} -B {input.bam} -S {input.splitters} -D {input.discordants} -o {output} &> {log}"


rule prepare_lumpy_discordant_alignment:
	input: lambda wc: lumpy_base_to_file_bam_path()[wc.base]
	output: os.path.join(lumpy_output_dir, "{base}" + ".discordants.unsorted.bam")
	message: "preparing discordant read alignments for lumpy from {input}"
	threads: 16
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["lumpy"]["conda"])
	log: os.path.join(lumpy_output_dir, "log", "{base}" + ".discordants.unsorted.bam.log")
	params: 
		samtools=tools_methods["lumpy"].get("samtools", "").get("path", "samtools")
	shell:
		"{params.samtools} view -@ {threads} -b -F 1294 -o {output} {input} &> {log}"


rule prepare_lumpy_split_alignment:
	input: lambda wc: lumpy_base_to_file_bam_path()[wc.base]
	output: os.path.join(lumpy_output_dir, "{base}" + ".splitters.unsorted.bam")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["lumpy"]["conda"])
	message: "extracting split reads from {input} to {output}"
	log: os.path.join(lumpy_output_dir, "log", "{base}" + ".splitters.unsorted.bam.log")
	params: 
		samtools=tools_methods["lumpy"].get("samtools", {}).get("path", "samtools"),
		script=tools_methods["lumpy"].get("extract_script", )
	shell:
		"{params.samtools} view -h {input} | {params.script} -i stdin | {params.samtools} view -o {output} -Sb - &> {log}"

rule prepare_lumpy_sort_alignments:
	input: os.path.join(lumpy_output_dir, "{base}" + ".{type}.unsorted.bam")
	output: os.path.join(lumpy_output_dir, "{base}" + ".{type}.sort.bam")
	message: "sorting lumpy alignment {input}"
	threads: 16
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["lumpy"]["conda"])
	log: os.path.join(lumpy_output_dir, "log", "{base}" + ".{type}.sort.bam.log")
	params:
		samtools=tools_methods["lumpy"].get("samtools", {}).get("path", "samtools")
	shell:
		"{params.samtools} sort -@ {threads} -o {output} {input} &> {log}"

