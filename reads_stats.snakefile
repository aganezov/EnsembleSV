configfile: "sv_tools.yaml"
configfile: "data.yaml"

tools_methods = config["tools_methods"]
reads_stats_output_dir = os.path.join(config["data_output"].get("reads_stats", {}).get("dir", "reads_stats"))

long_read_bams = config["data_input"].get("bams", {}).get("long", [])
long_reads_bams_by_bases = {os.path.basename(name).split(".")[0]: name for name in long_read_bams}
long_read_bases = list(long_reads_bams_by_bases.keys())


def reads_alignments_stats():
    result = []
    for long_read_base in long_read_bases:
        result.append(os.path.join(reads_stats_output_dir, long_read_base + ".reads_stats.alignment.txt"))
    return result


def reads_query_stats():
    reads_long_bases = list(config["data_input"].get("reads", {}).keys())
    result = []
    for long_read_base in long_read_bases:
        if long_read_base not in reads_long_bases:
            continue
        result.append(os.path.join(reads_stats_output_dir, long_read_base + ".reads_stats.query.txt"))
    return result


rule all:
    input: reads_alignments_stats(),
         reads_query_stats()

rule alignment_stats:
    input: lambda wc: long_reads_bams_by_bases[wc.base]
    output: os.path.join(reads_stats_output_dir, "{base}.reads_stats.alignment.txt")
    conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["reads_stats"]["conda"])
    log: os.path.join(reads_stats_output_dir, "log", "{base}.reads_stats.alignment.txt.log")
    params:
          median_flag=lambda wc: "--no-median" if not tools_methods["reads_stats"]["median"] else "",
          N50_flag=lambda wc: "--no-N50" if not tools_methods["reads_stats"]["N50"] else "",
          min_truncated_length=tools_methods["reads_stats"]["min_truncated_length"],
          max_truncated_length=tools_methods["reads_stats"]["max_truncated_length"],
          python=lambda wc: tools_methods["reads_stats"]["python"],
          script=tools_methods["reads_stats"]["path"],
    shell:
          "{params.python} {params.script} alignment {input} {params.median_flag} {params.N50_flag} --truncation-min {params.min_truncated_length} --truncation-max {params.max_truncated_length} -o {output} &> {log}"

rule query_stats:
    input: lambda wc: config["data_input"]["reads"][wc.base]
    output: os.path.join(reads_stats_output_dir, "{base}.reads_stats.query.txt")
    conda: os.path.join(config["tools_methods_conda_dir"], tools_methods["reads_stats"]["conda"])
    log: os.path.join(reads_stats_output_dir, "log", "{base}.reads_stats.query.txt.log")
    params:
          median_flag=lambda wc: "--no-median" if not tools_methods["reads_stats"]["median"] else "",
          N50_frag=lambda wc: "--no-N50" if not tools_methods["reads_stats"]["N50"] else "",
          min_truncated_length=tools_methods["reads_stats"]["min_truncated_length"],
          max_truncated_length=tools_methods["reads_stats"]["max_truncated_length"],
          python=lambda wc: tools_methods["reads_stats"]["python"],
          script=tools_methods["reads_stats"]["path"],
    shell:
          "{params.python} {params.script} fastq {input} {params.median_flag} {params.N50_flag} --truncation-min {params.min_truncated_length} --truncation-max {params.max_truncated_length} -o {output} &> {log}"
