configfile: "sv_tools.yaml"
configfile: "data.yaml"

tools_methods = config["tools_methods"]
raw_sv_calls_dir = config["data_output"].get("raw_sv_calls", {}).get("dir", "./raw")
sniffles_output_dir = os.path.join(raw_sv_calls_dir, config["data_output"].get("sniffles", {}).get("dir", "sniffles"))


def expected_sniffles_result_sv_files():
	long_read_bams = config["data_input"]["bams"]["long"]
	long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
	return [os.path.join(raw_sv_calls_dir, "raw", base + "_sniffles.vcf") for base in long_read_bases]


def sniffles_bases():
	return [os.path.basename(name).split(".")[0] for name in config["data_input"]["bams"]["long"]]

def sniffles_base_to_bam_file_path():
	return {os.path.basename(name).split(".")[0]: name for name in config["data_input"]["bams"]["long"]}


rule run_sniffles_all:
	input:
		expected_sniffles_result_sv_files()

rule copy_sniffles_vcf_to_aggregate_dir:
	input: os.path.join(sniffles_output_dir, "{base}" + "_sniffles.vcf")
	output: os.path.join(raw_sv_calls_dir, "raw", "{base}" + "_sniffles.vcf")
	message: "copying sniffles sv calls file {input} to the aggregate raw sv calls dir"
	shell:
		"cp {input} {output}"

rule run_sniffles:
	input: lambda wc: sniffles_base_to_bam_file_path()[wc.base]
	output: os.path.join(sniffles_output_dir, "{base}" + "_sniffles.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["sniffles"]["conda"])
	params:
		sniffles=tools_methods["sniffles"].get("path", "sniffles"),
		min_support=tools_methods["sniffles"].get("min_support", 2),
		max_num_splits=tools_methods["sniffles"].get("max_num_splits", 10),
		max_distance=tools_methods["sniffles"].get("max_distance", 1000),
		min_length=tools_methods["sniffles"].get("min_length", 30),
		num_reads_report=tools_methods["sniffles"].get("num_reads_report", -1),
		min_seq_size=tools_methods["sniffles"].get("min_seq_size", 1000)
	threads: 16
	message: "running sniffles with {threads} threads on {input}"
	shell:
		"{params.sniffles} -m {input} -v {output} --threads {threads} --min_support {params.min_support} --max_distance {params.max_distance} --max_num_splits {params.max_num_splits} --min_length {params.min_length} --num_reads_report {params.num_reads_report} --min_seq_size {params.min_seq_size} --report_seq"