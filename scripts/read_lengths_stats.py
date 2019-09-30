import argparse
import bisect
import statistics
import sys

import pysam


def get_size_bins(bins_strs):
    result = [-3000000000]
    for str_value in bins_strs:
        result.append(int(float(str_value)))
    result.append(3000000000)
    return result


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="commands", dest="command")
    shared_parser = argparse.ArgumentParser(add_help=False)
    shared_parser.add_argument("--bins", type=str, default="1500,3000,4500,6000,7500,9000,10500,12000,13500,15000,16500,18000,19500,21000,25000,30000,35000,40000,50000,60000")
    shared_parser.add_argument("--no-median", action="store_false", dest="median")
    shared_parser.add_argument("--no-N50", action="store_false", dest="compute_N50")
    shared_parser.add_argument("--truncation-min", type=int, default=1500, dest="trunc_min")
    shared_parser.add_argument("--truncation-max", type=int, default=1000000000, dest="trunc_max")
    shared_parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    subparsers.required = True
    fastq_parser = subparsers.add_parser("fastq", parents=[shared_parser])
    fastq_parser.add_argument("fastq", nargs="+")
    alignment_parser = subparsers.add_parser("alignment", parents=[shared_parser])
    alignment_parser.add_argument("alignment", nargs="+")
    alignment_parser.add_argument("--format", choices=["bam", "sam", "cram"], default="bam")
    args = parser.parse_args()
    bins = get_size_bins(bins_strs=args.bins.split(","))
    if args.command == "fastq":
        all_lengths = []
        lengths = {bin_size: 0 for bin_size in bins}
        for source in args.fastq:
            with pysam.FastxFile(source) as fh:
                for entry in fh:
                    length = len(entry.sequence)
                    all_lengths.append(length)
                    bin_size = bins[bisect.bisect_right(bins, length)]
                    lengths[bin_size] += 1
        print("bin size", "raw", sep=",", file=args.output)
        for bin_size in bins:
            print(bin_size, lengths.get(bin_size, 0), sep=",", file=args.output)
        truncated_lengths = [length for length in all_lengths if args.trunc_min <= length <= args.trunc_max]
        if args.median:
            print(file=args.output)
            print("median length (all):", statistics.median(all_lengths), file=args.output)
            print("median length (truncated):", statistics.median(truncated_lengths), file=args.output)
        if args.compute_N50:
            print(file=args.output)
            total_length = sum(all_lengths)
            total_truncated_length = sum(truncated_lengths)
            all_lengths = sorted(all_lengths, reverse=True)
            truncated_lengths = sorted(truncated_lengths, reverse=True)
            com_sum = 0
            for length in all_lengths:
                com_sum += length
                if com_sum >= total_length / 2:
                    print("N50 (all):", length, file=args.output)
                    break
            com_sum = 0
            for length in truncated_lengths:
                com_sum += length
                if com_sum >= total_truncated_length / 2:
                    print("N50 (truncated):", length, file=args.output)
    elif args.command == "alignment":
        all_query_lengths = []
        all_alignment_lengths = []
        query_lengths = {bin_size: 0 for bin_size in bins}
        alignment_lengths = {bin_size: 0 for bin_size in bins}
        format = "r"
        if args.format == "bam":
            format += "b"
        elif args.format == "cram":
            format += "c"
        processed_reads = set()
        for source in args.alignment:
            with pysam.AlignmentFile(source, format) as i_stream:
                for entry in i_stream:
                    query_id = entry.query_name
                    if query_id in processed_reads:
                        continue
                    processed_reads.add(query_id)
                    query_length = entry.query_length
                    all_query_lengths.append(query_length)
                    alignment_length = entry.query_alignment_length
                    all_alignment_lengths.append(alignment_length)
                    query_length_bin = bins[bisect.bisect_right(bins, query_length)]
                    alignment_length_bin = bins[bisect.bisect_right(bins, alignment_length)]
                    query_lengths[query_length_bin] += 1
                    alignment_lengths[alignment_length_bin] += 1
        print("bin size", "raw-aligned", "aligned", sep=",", file=args.output)
        for bin_size in bins:
            print(bin_size, query_lengths.get(bin_size, 0), alignment_lengths.get(bin_size, 0), sep=",", file=args.output)
        truncated_query_lengths = [length for length in all_query_lengths if args.trunc_min <= length <= args.trunc_max]
        truncated_alignment_lengths = [length for length in all_alignment_lengths if args.trunc_min <= length <= args.trunc_max]
        if args.median:
            print(file=args.output)
            print("median query length (all):", statistics.median(all_query_lengths), file=args.output)
            print("median query length (truncated):", statistics.median(truncated_query_lengths), file=args.output)
            print("median query length (all):", statistics.median(all_alignment_lengths), file=args.output)
            print("median query length (truncated):", statistics.median(truncated_alignment_lengths), file=args.output)
        if args.compute_N50:
            print(file=args.output)
            all_query_length = sum(all_query_lengths)
            all_query_lengths = sorted(all_query_lengths, reverse=True)
            truncated_query_length = sum(truncated_query_lengths)
            truncated_query_lengths = sorted(truncated_query_lengths, reverse=True)
            all_alignment_length = sum(alignment_lengths)
            all_alignment_lengths = sorted(alignment_lengths, reverse=True)
            truncated_alignment_lengths = sorted(truncated_alignment_lengths, reverse=True)
            truncated_alignment_length = sum(truncated_alignment_lengths)
            for name, total_length, lengths in [("all query", all_query_length, all_query_lengths),
                                                ("truncated query", truncated_query_length, truncated_query_lengths),
                                                ("all alignment", all_alignment_length, all_alignment_lengths),
                                                ("truncated alignment", truncated_alignment_length, truncated_alignment_lengths)]:
                com_sum = 0
                for length in lengths:
                    com_sum += length
                    if com_sum >= total_length / 2:
                        print("N50 ({name}):".format(name=name), length, file=args.output)
                    break


if __name__ == "__main__":
    main()
