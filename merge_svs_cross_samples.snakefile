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
cross_samples_to_refs = {}
cross_samples_to_alts = {}
cross_samples_to_gt =   {}

for sample, data in config["data_input"]["cross_samples"].items():
    base = data.get("EnsembleSV", {}).get("base", sample)
    cross_samples_to_bases[sample] = base
    cross_samples_to_refs[sample] = data.get("ref_field", "ref")
    cross_samples_to_alts[sample] = data.get("alt_field", "alt")
    cross_samples_to_gt[sample] = data.get("gt_field", "or_gt")

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
            for stats_type in ["svtypes", "samples"]:
                result.append(os.path.join(cross_samples_output_dir_stats, f"{sample}.{suffix}.{stats_type}.txt"))
    for suffix in ["spes", "sens"]:
        for stats_type in ["svtypes", "samples"]:
            result.append(os.path.join(cross_samples_output_dir_stats, f"{exp_name}.{suffix}.{stats_type}.txt"))
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
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    log:    os.path.join(cross_samples_output_dir_vcf, "log", "{sample}.{suffix}.rck.vcf.log")
    params:
        rck_adj_rck2x=tools_methods["rck"]["rck_adj_rck2x"]["path"],
        dummy_clone=lambda wc: wc.sample + "_" + wc.suffix,
        ref_extra=lambda wc: cross_samples_to_refs[wc.sample] if wc.sample in cross_samples_to_refs else ",".join(cross_samples_to_refs.values()),
        alt_extra=lambda wc: cross_samples_to_refs[wc.sample] if wc.sample in cross_samples_to_alts else ",".join(cross_samples_to_alts.values()),
        gt_dummy_extra=lambda wc: cross_samples_to_gt[wc.sample] if wc.sample in cross_samples_to_gt else ",".join(cross_samples_to_gt.values()),
    shell:
        "{params.rck_adj_rck2x} vcf-sniffles {input} --dummy-clone {params.dummy_clone} --dummy-clone-gt-extra {params.gt_dummy_extra} -o {output} --ref-extra {params.ref_extra} --alt-extra {params.alt_extra} &> {log}"

rule svtype_stats: # for all SV sets
    # simple RCK powered stats
    input:  os.path.join(cross_samples_output_dir_rck, "{sample}.{suffix}.rck.adj.tsv")
    output: os.path.join(cross_samples_output_dir_stats, "{sample}.{suffix}.svtypes.txt")
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    log: os.path.join(cross_samples_output_dir_stats, "log", "{sample}.{suffix}.svtypes.txt.log")
    params:
        rck_adj_stats=tools_methods["rck"]["rck_adj_stats"]["path"],
    shell:
        "{params.rck_adj_stats} survivor-stat {input} --sources-field svtype -o {output} &> {log}"

rule samples_stats:
    input:  os.path.join(cross_samples_output_dir_rck, "{sample}.{suffix}.rck.adj.tsv")
    output: os.path.join(cross_samples_output_dir_stats, "{sample}.{suffix}.samples.txt")
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    log: os.path.join(cross_samples_output_dir_stats, "log", "{sample}.{suffix}.samples.txt.log")
    params:
        rck_adj_stats=tools_methods["rck"]["rck_adj_stats"]["path"],
        source_field=lambda wc: wc.sample.lower() + "_supporting_sources",
    shell:
        "{params.rck_adj_stats} survivor-stat {input} --sources-field {params.source_field} -o {output} &> {log}"

rule unique_rck:    # only for original samples
    # based on supporting_ids_field
    input:  original=os.path.join(cross_samples_output_dir_rck, "{sample}.spes.rck.adj.tsv"),
            exp_sens=os.path.join(cross_samples_output_dir_rck, exp_name + ".sens.rck.adj.tsv"),
    output: os.path.join(cross_samples_output_dir_rck, "{sample," + samples_regex + "}.unique.rck.adj.tsv")
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    log:    os.path.join(cross_samples_output_dir_rck, "log", "{sample," + samples_regex + "}.unique.rck.adj.tsv.log")
    run:
        import sys
        with open(log[0], "w") as stdout:
            sys.stdout = stdout
            sys.stderr = stdout
            from collections import defaultdict
            from rck.core.io import EXTERNAL_NA_ID, write_adjacencies_to_file, stream_adjacencies_from_source
            exp_sens_ids_to_sources = defaultdict(list)
            with open(input.exp_sens, "rt") as source:
                exp_sens = stream_adjacencies_from_source(source=source)
                for adj in exp_sens:
                    supporting_ids = adj.extra.get(exp_name.lower() + "_supporting_source_ids").split(",")
                    supporting_sources = adj.extra.get(exp_name.lower() + "_supporting_sources").split(",")
                    for aid_id in supporting_ids:
                        exp_sens_ids_to_sources[aid_id] = supporting_sources
            with open(input.original, "rt") as source:
                original = stream_adjacencies_from_source(source=source)
                unique = (adj for adj in original if len(exp_sens_ids_to_sources.get(adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased), [])) == 1)
                write_adjacencies_to_file(file_name=output[0], adjacencies=unique)

rule annotated_rck: # only for original samples
    # reading original rck and experiment merged sens rck, and recording support_ids and sources as annotations
    input:  original=lambda wildcards: original_rck(sample=wildcards.sample, suffix=wildcards.suffix),
            exp_sens=os.path.join(cross_samples_output_dir_rck, exp_name + ".sens.rck.adj.tsv"),
    output: os.path.join(cross_samples_output_dir_rck, "{sample," + samples_regex + "}.{suffix}.rck.adj.tsv")
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    log:    os.path.join(cross_samples_output_dir_rck, "log", "{sample," + samples_regex + "}.{suffix}.rck.adj.tsv")
    params:
        ann_prefix=lambda wc: wc.sample + "_" + wc.suffix + "_" + exp_name,
    run:
        import sys
        from collections import defaultdict
        import itertools
        from rck.core.io import EXTERNAL_NA_ID, write_adjacencies_to_file, stream_adjacencies_from_source

        def stream_annotated_adjs(input_adjacencies, merged_with_dict):
            for adj in input_adjacencies:
                merged_with = merged_with_dict.get(adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased), [])
                if len(merged_with) > 0:
                    adj.extra[params.ann_prefix + "_merged_with"] = ",".join(set(merged_with))
                yield adj

        with open(log[0], "w") as stdout:
            sys.stdout = stdout
            sys.stderr = stdout
            exp_sens_supporting_ids_pairing = defaultdict(list)
            with open(input.exp_sens, "rt") as source:
                exp_sens = stream_adjacencies_from_source(source=source)
                for adj in exp_sens:
                    supporting_ids = adj.extra.get(exp_name.lower() + "_supporting_source_ids").split(",")
                    for id1, id2 in itertools.permutations(supporting_ids, r=2):
                        exp_sens_supporting_ids_pairing[id1].append(id2)
            with open(input.original, "rt") as source:
                original = stream_adjacencies_from_source(source=source)
                annotated_original = stream_annotated_adjs(input_adjacencies=original, merged_with_dict=exp_sens_supporting_ids_pairing)
                write_adjacencies_to_file(file_name=output[0], adjacencies=annotated_original, sort_adjacencies=False)


rule spes_experiment_rck: # for overall experiment sv set only
    # reading all merged SVs and retaining only those, that have supporting sample-spes SV origins
    output: os.path.join(cross_samples_output_dir_rck, exp_name + ".spes.rck.adj.tsv")
    input: sens=os.path.join(cross_samples_output_dir_rck, exp_name + ".sens.rck.adj.tsv"),
           sample_spes=lambda wc: all_original_rcks(suffix="spes")
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    log:    os.path.join(cross_samples_output_dir_rck, "log", exp_name + ".spes.rck.adj.tsv")
    run:
        import sys

        def spes_filter(adjacencies, spes_original_ids):
            for adj in adjacencies:
                supporting_ids = adj.extra.get(exp_name.lower() + "_supporting_source_ids").split(",")
                for sup_id in supporting_ids:
                    if sup_id in spes_original_ids:
                        yield adj

        with open(log[0], "w") as stdout:
            sys.stdout = stdout
            sys.stderr = stdout
            from rck.core.io import EXTERNAL_NA_ID, write_adjacencies_to_file
            spes_original_adjs_ids = set()
            for sample_file in input.sample_spes:
                with open(sample_file, "rt") as source:
                    adjacencies = stream_adjacencies_from_source(source=source)
                    for adj in adjacencies:
                        spes_original_adjs_ids.add(adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased))
            with open(input.sens, "rt") as source:
                exp_sensitive_adjacencies = stream_adjacencies_from_source(source=source)
                exp_spes_adj = spes_filter(adjacencies=exp_sensitive_adjacencies, spes_original_ids=spes_original_adjs_ids)
                write_adjacencies_to_file(file_name=output[0], adjacencies=exp_spes_adj)

rule sens_experiment_rck: # for overall experiment sv set only
    # regular RCK powered survivor to rck conversion
    output: os.path.join(cross_samples_output_dir_rck, exp_name + ".sens.rck.adj.tsv")
    input:  sample_sens=lambda wc: all_original_rcks(suffix="sens"),
            survivor_vcf=os.path.join(cross_samples_output_dir_survivor, exp_name + ".sens.survivor.vcf")
    conda:  os.path.join(config["tools_methods_conda_dir"], tools_methods["rck"]["conda"])
    log:    os.path.join(cross_samples_output_dir_rck, "log", exp_name + ".sens.rck.adj.tsv.log")
    params:
        rck_adj_x2rck=tools_methods["rck"]["rck_adj_x2rck"]["path"],
        samples=lambda wc: ",".join(sample + "_sens" for sample in sorted(cross_samples_to_bases.keys())),
        samples_source=lambda wc: ",".join([original_rck(sample=sample, suffix="sens") for sample in sorted(cross_samples_to_bases.keys())]),
        suffix=lambda wc: exp_name.lower(),
        sample_name_prefix=lambda wc: "--samples-suffix-extra" if config["data_merge"].get("sample_name_prefix", False) else ""
    shell:
        "{params.rck_adj_x2rck} survivor {input.survivor_vcf} --id-suffix {params.suffix} --samples {params.samples} {params.sample_name_prefix} --samples-source {params.samples_source} --survivor-prefix {params.suffix} -o {output} &> {log}"

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
        min_sv_size=config["data_merge"]["survivor"].get("min_length", 30)
    shell:
        "{params.survivor} merge {input} {params.max_distance} {params.min_caller_cnt} {params.sv_type_consider} {params.sv_strands_consider} {params.distance_estimate} {params.min_sv_size} {output} &> {log}"


rule sens_experiment_prepare_survivor:
    output: os.path.join(cross_samples_output_dir_survivor, exp_name + ".sens.survivor")
    input:  original_sens_rcks=[os.path.join(cross_samples_output_dir_survivor, f"{sample}.sens.survivor.vcf") for sample in cross_samples_to_bases.keys()]
    log:    os.path.join(cross_samples_output_dir_survivor, "log", exp_name + ".sens.survivor.log")
    run:
        import sys
        with open(log[0], "w") as stdout:
            sys.stdout = stdout
            sys.stderr = stdout
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
