import argparse
import sys
from collections import defaultdict

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=argparse.FileType("rt"), default=sys.stdin)
    parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("--min-support", required=True, type=int)
    args = parser.parse_args()

    cnt = defaultdict(int)
    for line_cnt, line in enumerate(args.input):
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            print(line, file=args.output)
            continue
        sv_type_entry = line.split("\t")[7].split(";")[0].split("=")
        if len(sv_type_entry) < 2:
            continue
        sv_type = sv_type_entry[1]
        if sv_type not in ["DUP", "DEL", "INS", "INV", "BND"]:
            continue
        info = line.split("\t")[9].split(":")
        f,s = info[0].split("/")
        fs, ss = info[1].split(",")
        fs, ss = int(fs), int(ss)
        ts = 0
        if f != "0":
            ts += fs
        if s != "0":
            ts += ss
        cnt[int(info[2])] += ts
        if ts >= args.min_support:
            print(line, file=args.output)

    # cnt = dict(cnt)
    # for key in sorted(cnt.keys()):
    #      print(key, "=", cnt[key])
