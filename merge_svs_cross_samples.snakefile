configfile: "data.yaml"
configfile: "cross_samples_data.yaml"
configfile: "sv_tools.yaml"
import os

tools_methods = config["tools_methods"]
cross_samples_output_dir = os.path.join(config["data_output"].get("cross_samples", {}).get("dir", "cross_samples"))
cross_samples_output_dir_vcf = os.path.join(cross_samples_output_dir, config["data_output"].get("cross_samples", {}).get("vcf", "vcf"))
cross_samples_output_dir_rck = os.path.join(cross_samples_output_dir, config["data_output"].get("cross_samples", {}).get("rck", "rck"))
cross_samples_output_dir_stats = os.path.join(cross_samples_output_dir, config["data_output"].get("cross_samples", {}).get("stats", "stats"))
cross_samples_output_dir_survivor = os.path.join(cross_samples_output_dir, config["data_output"].get("cross_samples", {}).get("survivor", "survivor"))
exp_name = config.get("data_experiment_name", "DEFAULT")

cross_samples_to_bases = {}
for sample, data in config["data_input"]["cross_samples"].items():
    base = data.get("EnsembleSV", {}).get("base", sample)
    cross_samples_to_bases[sample] = base

samples_regex = long_methods_regex = "(" + "|".join(cross_samples_to_bases.keys()) + ")"


def expected_rck_files():
    result = []
    for sample, base in cross_samples_to_bases.items():
        for suffix in ["spes", "sens", "unique"]:
            result.append(os.path.join(cross_samples_output_dir_rck, f"{sample}.{suffix}.rck.adj.tsv"))
    for suffix in ["spes", "sens"]:
        result.append(os.path.join(cross_samples_output_dir_rck, f"{exp_name}.{suffix}.rck.adj.tsv"))
    return result


def expected_vcf_files():
    result = []
    for sample, base in cross_samples_to_bases.items():
        for suffix in ["spes", "sens", "unique"]:
            result.append(os.path.join(cross_samples_output_dir_vcf, f"{sample}.{suffix}.rck.vcf"))
    for suffix in ["spes", "sens"]:
        result.append(os.path.join(cross_samples_output_dir_vcf, f"{exp_name}.{suffix}.rck.vcf"))
    return result


def expected_stats_files():
    result = []
    for sample, base in cross_samples_to_bases.items():
        for suffix in ["spes", "sens", "unique"]:
            for stats_type in ["svtypes"]:
                result.append(os.path.join(cross_samples_output_dir_stats, f"{sample}.{suffix}.{stats_type}.txt"))
    for suffix in ["spes", "sens"]:
        for stats_type in ["svtypes"]:
            result.append(os.path.join(cross_samples_output_dir_vcf, f"{exp_name}.{suffix}.{stats_type}.txt"))
    return result


def sample_rck_dir_path(sample):
    return os.path.join(config["data_input"]["cross_samples"][sample].get("EnsembleSV", {}).get("dir", sample + "_EnsembleSV"),
                        config["data_input"]["cross_samples"][sample].get("EnsembleSV", {}).get("rck_sub", "rck"))


def original_rck(sample, suffix="spes"):
    return os.path.join(sample_rck_dir_path(sample=sample), f'{cross_samples_to_bases[sample]}.{suffix}.rck.adj.tsv')

def all_original_rcks(suffix="spes"):
    return [original_rck(sample=sample, suffix=suffix) for sample in cross_samples_to_bases.keys()]












rule all:
    input:
         expected_rck_files(),
         expected_vcf_files(),
         expected_stats_files(),

rule rck_to_vcf: # for all SV sets
    # simple conversion to VCF from RCK
    input:  os.path.join(cross_samples_output_dir_rck, "{sample}.{suffix}.rck.adj.tsv")
    output: os.path.join(cross_samples_output_dir_vcf, "{sample}.{suffix}.rck.vcf")

rule stats: # for all SV sets
    # simple RCK powered stats
    input:  os.path.join(cross_samples_output_dir_rck, "{sample}.{suffix}.rck.adj.tsv")
    output: os.path.join(cross_samples_output_dir_stats, "{sample}.{suffix}.{svtype}.txt")

rule unique_rck:    # only for original samples
    # based on supporting_ids_field
    input:  original_rck=os.path.join(cross_samples_output_dir_rck, "{sample}.spes.rck.adj.tsv"),
            exp_sens=os.path.join(cross_samples_output_dir_rck, exp_name + ".sens.rck.adj.tsv"),
    output: os.path.join(cross_samples_output_dir_rck, "{sample," + samples_regex + "}.unique.rck.adj.tsv")

rule annotated_rck: # only for original samples
    # reading original rck and experiment merged sens rck, and recording support_ids and sources as annotations
    input:  original=lambda wildcards: original_rck(sample=wildcards.sample, suffix=wildcards.suffix),
            exp_sens=os.path.join(cross_samples_output_dir_rck, exp_name + ".sens.rck.adj.tsv"),
    output: os.path.join(cross_samples_output_dir_rck, "{sample," + samples_regex + "}.{suffix}.rck.adj.tsv")

rule spes_experiment_rck: # for overall experiment sv set only
    # reading all merged SVs and retaining only those, that have supporting sample-spes SV origins
    output: os.path.join(cross_samples_output_dir_rck, exp_name + ".spes.rck.adj.tsv")
    input: sens=os.path.join(cross_samples_output_dir_rck, exp_name + ".sens.rck.adj.tsv"),
           sample_spes=lambda wc: all_original_rcks(suffix="spes")

rule sens_experiment_rck: # for overall experiment sv set only
    # regular RCK powered survivor to rck conversion
    output: os.path.join(cross_samples_output_dir_rck, exp_name + ".sens.rck.adj.tsv")
    input:  sample_sens=lambda wc: all_original_rcks(suffix="sens"),
            survivor_vcf=os.path.join(cross_samples_output_dir_survivor, exp_name + ".sens.survivor.vcf")

rule sens_experiment_run_survivor:
    output: os.path.join(cross_samples_output_dir_survivor, exp_name + ".sens.survivor.vcf")
    input:  os.path.join(cross_samples_output_dir_survivor, exp_name + ".sens.survivor")
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["survivor"]["conda"])
    log:    os.path.join(cross_samples_output_dir_survivor, "log", exp_name + ".sens.survivor.vcf.log")
    params:
        survivor=tools_methods["survivor"]["path"],
        max_distance=config["data_merge"]["survivor"].get("max_distance", 1000),
        min_caller_cnt=1,
        sv_type_consider=config["data_merge"]["survivor"].get("svtype", 1),
        sv_strands_consider=config["data_merge"]["survivor"].get("strands", 1),
        distance_estimate=0,
        min_sv_size=config["data_merge"]["survivor"].get("max_distance", 30)
    shell:
        "{params.survivor} merge {input} {params.max_distance} {params.min_caller_cnt} {params.sv_type_consider} {params.sv_strands_consider} {params.distance_estimate} {params.min_sv_size} {output} &> {log}"


rule sens_experiment_prepare_survivor:
    output: os.path.join(cross_samples_output_dir_survivor, exp_name + ".sens.survivor")
    input:  original_sens_rcks=[os.path.join(cross_samples_output_dir_survivor, f"{sample}.sens.survivor.vcf") for sample in cross_samples_to_bases.keys()]
    run:
        with open(output[0], "wt") as dest:
            for file_name in input:
                print(file_name, file=dest)

rule sens_original_rck_survivor_suitable:
    output: os.path.join(cross_samples_output_dir_survivor, "{sample}.sens.survivor.vcf")
    input:  lambda wc: original_rck(sample=wc.sample, suffix="sens")
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    log:    os.path.join(cross_samples_output_dir_survivor, "log", "{sample}.sens.survivor.vcf.log")
    params:
        rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
        dummy_clone=lambda wc: wc.sample + "_sens",
        ref_extra="ref",
        alt_extra="alt",
    shell:
        "{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} -o {output} --ref-extra {params.ref_extra} --alt-extra {params.alt_extra} &> {log}"


rule generate_original_sample_sens_spes_rck_from_generic:
    output: os.path.join("{sample}_EnsembleSV", "rck", "{sample}.{suffix}.rck.adj.tsv")
    input:  lambda wc: config["data_input"]["cross_samples"].get(wc.sample, {}).get("rck", {}).get("path", cross_samples_to_bases[wc.sample] + "rck.adj.tsv")

rule generate_original_sample_generic_rck_from_vcf:
    output: os.path.join("{sample}_EnsembleSV", "rck", "{sample}.rck.adj.tsv")
    input:  lambda wc: config["data_input"]["cross_samples"].get(wc.sample, {}).get("vcf", {}).get("path", cross_samples_to_bases[wc.sample] + "vcf")


ruleorder: unique_rck > annotated_rck









rule spes_to_unqiue:
    output: os.path.join(cross_samples_output_dir, "{base}.spes.unique.rck.adj.tsv")
    input: spes=lambda wc: os.path.join(config["cross_samples"][wc.base]["EnsembleSV_path"], "merged_sv_calls", "rck", wc.base + ".spes.rck.adj.tsv"),
         survivor_merge=os.path.join(cross_samples_output_dir, "survivor", "{base}.vRestSens.vcf"),
    log: os.path.join(cross_samples_output_dir, "log", "{base}.spes.unique.rck.adj.tsv.log")
    conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    run:
        import vcf
        from rck.core.io import read_adjacencies_from_file, EXTERNAL_NA_ID, write_adjacencies_to_file

        unique_adj_ids = set()
        with open(input.survivor_merge, "rt") as source:
            vcf_reader = vcf.Reader(source)
            if "RNAMES" in vcf_reader.infos:
                vcf_reader.infos["RNAMES"] = vcf_reader.infos["RNAMES"]._replace(num='.')
            for record in vcf_reader:
                supp_vec = svtype = record.INFO.get("SUPP_VEC", "01")
                if supp_vec[0] == "0" or "1" in supp_vec[1:]:
                    continue
                for vcf_sample in record.samples:
                    adj_id = vcf_sample.data.ID.split(",")[0]
                    if adj_id == "NaN":
                        continue
                    unique_adj_ids.add(adj_id)
        adjs = read_adjacencies_from_file(input.spes)
        adjacencies = [adj for adj in adjs if adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased) in unique_adj_ids]
        write_adjacencies_to_file(output[0], adjacencies)

rule survivor_merge:
    input: os.path.join(cross_samples_output_dir, "survivor", "{base}.vRestSens.survivor"),
    output: os.path.join(cross_samples_output_dir, "survivor", "{base}.vRestSens.vcf"),
    params:
          survivor="SURVIVOR",
          max_distance=1000,
          min_caller_cnt=1,
          sv_type_consider=0,
          sv_strands_consider=1,
          distance_estimate=0,
          min_sv_size=30,
    shell:
         "{params.survivor} merge {input} {params.max_distance} {params.min_caller_cnt} {params.sv_type_consider} {params.sv_strands_consider} {params.distance_estimate} {params.min_sv_size} {output}"

rule survivor_merge_config:
    output: os.path.join(cross_samples_output_dir, "survivor", "{base}.vRestSens.survivor"),
    input: spes=lambda wc: os.path.join(config["cross_samples"][wc.base]["EnsembleSV_path"], "merged_sv_calls", "merged", wc.base + ".spes.rck.vcf"),
         rest_sens=lambda wc: [os.path.join(config["cross_samples"][base]["EnsembleSV_path"], "merged_sv_calls", "merged", base + ".sens.rck.vcf") for base in cross_samples_bases
                               if base != wc.base],
    run:
        with open(output[0], "wt") as dest:
            print(input.spes, file=dest)
            for file_name in input.rest_sens:
                print(file_name, file=dest)
