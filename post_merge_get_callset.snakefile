configfile: "sv_tools.yaml"
configfile: "data.yaml"

import os
from colelctions import defaultdict

tools_methods = config["tools_methods"]
merged_dir = config["data_output"]["merge"]["dir"]
pre_merge_aggregate_dir = config["data_output"]["pre_merge_process"]["dir"]


method_to_data_type = defaultdict(list) 
bases_by_data_type = defaultdict(set)
for data_type in config["tools_read_type_to_method"]:
	for method in config["tools_read_type_to_method"][data_type]:
		method_to_data_type[method].append(data_type)
	bams = config.get("data_input", {}).get("bams", {}).get(data_type, [])
	for bam in bams:
		base = os.path.basename(name).split(".")[0]
		bases_by_data_type[data_type].append(base)





class MergeGroup(object):
	def __init__(self, name, min_cnt, entries=None):
		self.name = name,
		self.min_cnt = min_cnt,
		self.entries = entries if entries is not None else {}
		self.sample_to_sample_name_dict = self.get_sample_to_sample_name_dict()

	def get_sample_to_sample_name_dict(self):
		result = {sample_name: sample_path for sample_name, sample_path in self.entries.get("hardcode", [])}
		methods = self.entries.get("methods", config["tools_enabled_methods"])
		data_types = set()
		for method in methods:
			for data_type in method_to_data_type[method]:
				data_types.add(data_type)
		defined_bases = self.entries.get("bases", set())
		method_based_bases = set()
		for data_type in data_types:
			for base in bases_by_data_type[data_type]:
				method_based_bases.add(base)
		if len(defined_bases) == 0:
			bases = method_based_bases
		else:
			bases = method_based_bases & defined_bases
		possible_files = [file_name for file_name in os.path.listdir(pre_merge_aggregate_dir)]
		for method in methods:
			for base in bases:
				min_len = config["data_premerge"]["min_len"].get(method_to_data_type[method], None)
				min_len_suffix = "" if min_len is None else ".m" + min_len
				min_support = config["data_premerge"]["min_support_fraction"].get(method_to_data_type[method], None)
				min_support_suffix = "" if min_support is None else ".ms" + min_support
				constructed_file_name = base + "_" + method + min_support_suffix + min_len_suffix + "rck.adj.tsv"

		


rule run_get_call_set_all:
	input: os.path.join(merged_dir, config["data_sample_name"] + "call_set.vcf")

rule run_convert_rck_to_vcf:
	output: os.path.join(merged_dir, "{sample_name}_callSet.vcf")
	input: os.path.join(merged_dir, "{sample_name}_callSet.rck.adj.tsv")
	message: "translating RCK formatted call set to VCF"
	conda: tools_methods["rck"]["conda"]
	params:
		rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
		dummy_clone=lambda wc: wc.sample_name + "_callSet",
	shell:
		"{params.rck_adj_rck2x} vcf-sniffles {input} {input} --dummy-clone {params.dummy_clone} -o {output}"

rule obtain_rck_call_set:
	input: os.path.join(merged_dir, "{sample_name}_callSet.yaml")
	output: os.path.join(merged_dir, "{sample_name}_callSet.rck.adj.tsv")
	message: "obtaining a merged call set"
	conda: tools_methods["rck"]["conda"]
	params: 
		call_set_script=tools_methods["rck"]["call_set_script"]["path"]
		output_dir=merged_dir,
		output_prefix=lambda wc: wc.sample_name
	shell:
		"{params.call_set_script} --config {input} --cosmic-segments {params.cosmic_segments} --output-dir {params.output_dir} --output-prefix {params.output_prefix}"

rule prepare_rck_call_set_config:
	output: os.path.join(merged_dir, "{sample_name}_callSet.yaml")
	input: 
		survivor_vcf=os.path.join(merged_dir, "{sample_name}.survivor.vcf")

