configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

tools_methods = config["tools_methods"]
raw_sv_calls_dir = config.get("data_output", {}).get("raw_sv_calls", {}).get("dir", "./raw")
naibr_output_dir = os.path.join(raw_sv_calls_dir, config["data_output"].get("naibr", {}).get("dir", "naibr"))

def expected_naibr_result_sv_files():
	if "linked" not in config["data_input"]["bams"]:
		return []
	linked_read_bams = config["data_input"]["bams"]["linked"]
	linked_read_bases = [os.path.basename(name).split(".")[0] for name in linked_read_bams]
	return [os.path.join(raw_sv_calls_dir, "raw", base + "_naibr.bedpe") for base in linked_read_bases]

def naibr_bases():
	return [os.path.basename(name).split(".")[0] for name in config["data_input"]["bams"]["linked"]]

def naibr_base_to_file_bam_path():
	return {os.path.basename(name).split(".")[0]: name for name in config["data_input"]["bams"]["linked"]}

rule run_naibr_all:
	input:
		expected_naibr_result_sv_files()

rule copy_naibr_svs_to_aggregate_dir:
	output: os.path.join(raw_sv_calls_dir, "raw", "{base}" + "_naibr.bedpe")
	input: os.path.join(naibr_output_dir, "{base}", "NAIBR_SVs.bedpe")
	message: "copying NAIBR SVs file {input} into aggregate dir"
	shell:
		"cp {input} {output}"

rule run_naibr_workflow:
	input: os.path.join(naibr_output_dir, "{base}", "{base}" + ".config")
	output: os.path.join(naibr_output_dir, "{base}", "NAIBR_SVs.bedpe")
	message: "executing NAIBR with config {input}"
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["naibr"]["conda"])
	threads: tools_methods["naibr"].get("threads", 15)
	log: os.path.join(naibr_output_dir, "{base}", "log", "NAIBR_SVs.bedpe.log")
	params:
		python=tools_methods["naibr"].get("python", {}).get("path", "python"),
		naibr_script=tools_methods["naibr"].get("script", {}).get("path", "NAIBR.py")
	shell:
		"{params.python} {params.naibr_script} {input} &> {log}"

rule create_naibr_config:
	input: lambda wc: naibr_base_to_file_bam_path()[wc.base]
	output: os.path.join(naibr_output_dir, "{base}", "{base}" + ".config")
	message: "creating naibr config for input {input}"
	params:
		min_mapq=tools_methods["naibr"].get("min_mapq", 40),
		base=lambda wc: wc.base
	threads: tools_methods["naibr"].get("threads", 16)
	run:
		with open(output[0], "wt") as destination:
			print("min_mapq={mapq}".format(mapq=params.min_mapq), file=destination)
			print("bam_file={bam}".format(bam=os.path.abspath(os.path.expanduser(input[0]))), file=destination)
			print("threads={threads}".format(threads=threads), file=destination)
			print("outdir={outdir}".format(outdir=os.path.abspath(os.path.join(naibr_output_dir, params.base))), file=destination)
