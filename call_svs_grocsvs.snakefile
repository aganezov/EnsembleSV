configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os

tools_methods = config["tools_methods"]
raw_sv_calls_dir = config["data_output"].get("raw_sv_calls", {}).get("dir", "./raw")
grocsvs_output_dir = os.path.join(raw_sv_calls_dir, config["data_output"].get("grocsvs", {}).get("dir", "grocsvs"))

def expected_grocsvs_result_sv_files():
	if "linked" not in config["data_input"]["bams"]:
		return []
	linked_read_bams = config["data_input"]["bams"]["linked"]
	linked_read_bases = [os.path.basename(name).split(".")[0] for name in linked_read_bams]
	return [os.path.join(raw_sv_calls_dir, "raw", base + "_grocsvs.vcf") for base in linked_read_bases]

def grocsvs_bases():
	return [os.path.basename(name).split(".")[0] for name in config["data_input"]["bams"]["linked"]]

def grocsvs_base_to_bam_file_path():
	result = {os.path.basename(name).split(".")[0]: name for name in config["data_input"]["bams"]["linked"]} 
	return result


rule run_grocsvs_all:
	input:
		expected_grocsvs_result_sv_files()

rule copy_grocsvs_vcf_to_aggregate_dir:
	input: os.path.join(grocsvs_output_dir, "{base}" + "_grocsvs.vcf")
	output: os.path.join(raw_sv_calls_dir, "raw", "{base}" + "_grocsvs.vcf")
	message: "copying grosvs sv calls file {input} to the aggregate raw sv calls dir"
	shell:
		"cp {input} {output}"

rule run_grocsvs_obtain_vcf:
	input: os.path.join(grocsvs_output_dir, "{base}", "results", "PostprocessingStep", "svs.vcf")
	output: os.path.join(grocsvs_output_dir, "{base}" + "_grocsvs.vcf")
	message: "Copying internal result to the top level"
	shell:
		"cp {input} {output}"


rule run_grocsvs_workflow:
	input: os.path.join(grocsvs_output_dir, "{base}", "{base}.config.json")
	output: os.path.join(grocsvs_output_dir, "{base}", "results", "PostprocessingStep", "svs.vcf")
	message: "executing grocsvs pipeline with {input} config"
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["grocsvs"]["conda"])
	params:
		grocsvs=tools_methods["grocsvs"].get("path", "grocsvs")
	shell:
		"{params.grocsvs} {input}"

rule prepare_grocsvs_config_file:
	input: 
		bam=lambda wc: grocsvs_base_to_bam_file_path()[wc.base],
		ref=config["data_input"]["fasta"]["ref"],
		bwa_index=config["data_input"]["index"]["ref"] + ".bwt"
	output: os.path.join(grocsvs_output_dir, "{base}", "{base}.config.json")
	message: "preparing config file {output} for {input.bam}"
	threads: 16 # tools_methods.get("grocsvs", {}).get("processes", 16)
	params:
		idba_ud=tools_methods["grocsvs"].get("idba_ud", {}).get("path", "idba_ud"),
		bwa=tools_methods["grocsvs"].get("bwa", {}).get("path", "bwa"),
		samtools=tools_methods["grocsvs"].get("samtools", {}).get("path", "samtools"),
		tabix=tools_methods["grocsvs"].get("tabix", {}).get("path", "tabix"),
		bgzip=tools_methods["grocsvs"].get("bgzip", {}).get("path", "gbzip"),
		base=lambda wc: wc.base,
		cluster_type=tools_methods["grocsvs"].get("cluster_type", "multiprocessing"),
		bwa_index=config["data_input"]["index"]["ref"]
	run:
		json_config = {
			"ref_fasta": input.ref,
			"bwa_index": params.bwa_index,
			"binaries": {
				"bwa": params.bwa,
				"idba_ud": params.idba_ud,
				"samtools": params.samtools,
				"tabix": params.tabix,
				"bgzip": params.bgzip,
			},
			"samples": {
				str(params.base): [
					{
						"type": "TenXDataset",
						"bam": input.bam,
						"id": params.base + "_grocsvs",
					},
				],
			},
			"cluster_settings" : {
				"cluster_type": params.cluster_type,
				"processes": threads,
			}
		}
		with open(output[0], "wt") as destination:
			json.dump(json_config, destination)

rule prepare_ref_bwa_index:
	input: config["data_input"]["fasta"]["ref"]
	output: config["data_input"]["fasta"]["ref"] + ".bwt"
	message: "preparing bwa index for reference"
	conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["grocsvs"]["conda"])
	params:
		bwa=tools_methods["grocsvs"].get("bwa", {}).get("path", "bwa")
	shell:
		"{params.bwa} index -a bwtsw {input}"












