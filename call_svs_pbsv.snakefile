configfile: "sv_tools.yaml"
configfile: "data.yaml"

tools_methods = config["tools_methods"]
raw_sv_calls_dir = config["data_output"].get("raw_sv_calls", {}).get("dir", "./raw")
pbsv_output_dir = os.path.join(raw_sv_calls_dir, config["data_output"].get("pbsv", {}).get("dir", "pbsv"))


def expected_pbsv_result_sv_files():
	if "long" not in config["data_input"]["bams"]:
		return []
	long_read_bams = config["data_input"]["bams"]["long"]
	long_read_bases = [os.path.basename(name).split(".")[0] for name in long_read_bams]
	result = [os.path.join(raw_sv_calls_dir, "raw", base + "_pbsv.vcf") for base in long_read_bases]

	return result

def pbsv_bases():
	return [os.path.basename(name).split(".")[0] for name in config["data_input"]["bams"]["long"]]

def pbsv_base_to_bam_file_path():
	return {os.path.basename(name).split(".")[0]: name for name in config["data_input"]["bams"]["long"]}



rule run_pbsv:
	input:
		expected_pbsv_result_sv_files()

rule copy_pbsv_vcf_to_aggregate_dir:
	input: os.path.join(pbsv_output_dir, "{base}" + "_pbsv.vcf")
	output: os.path.join(raw_sv_calls_dir, "raw", "{base}" + "_pbsv.vcf")
	message: "copying pbsv sv calls file {input} to the aggregate raw sv calls dir"
	shell:
		"cp {input} {output}"


rule run_pbsv_call:
	input: ref=config["data_input"]["fasta"]["long_ref"],
		   signatures=os.path.join(pbsv_output_dir, "{base}_pbsv.svsig.gz")
	output: os.path.join(pbsv_output_dir, "{base}_pbsv.vcf")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["pbsv"]["conda"])
	params: 
		pbsv=tools_methods["pbsv"].get("path", "pbsv"),
	message: "running pbsv call on {input.signatures} with ref {input.ref}"
	shell:
		"{params.pbsv} call {input.ref} {input.signatures} {output}"


rule run_pbsv_discover:
	input: lambda wc: pbsv_base_to_bam_file_path()[wc.base]
	output: os.path.join(pbsv_output_dir, "{base}_pbsv.svsig.gz")
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["pbsv"]["conda"])
	params:
		pbsv=tools_methods["pbsv"].get("path", "pbsv"),
		base=lambda wc: wc.base
	message: "running pbsv discover with "
	shell:
		"{params.pbsv} discover --sample {params.base} {input} {output}"