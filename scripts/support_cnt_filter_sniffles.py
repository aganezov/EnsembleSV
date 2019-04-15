import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=argparse.FileType("rt"), default=sys.stdin)
    parser.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("--min-support", required=True, type=int)
    args = parser.parse_args()

    for line in args.input:
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            print(line, file=args.output)
            continue
        info = line.split("\t")[7].split(";")
        for entry in info:
            e_data = entry.split("=")
            if e_data[0] == "RE" and int(e_data[1]) >= args.min_support:
                print(line, file=args.output)
                break
