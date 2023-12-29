import argparse
import sys


def main(args):
    if args.output is None:
        fout = sys.stdout
    else:
        fout = open(args.output, 'w')

    with open(args.input, 'r') as fin:
        line = fin.readline()
        words = line.split()
        while len(words) == 0:
            line = fin.readline()
            words = line.split()
        ntax = int(words[0])

    nkeep = ntax + 1
    nskip = nkeep * (args.gnum - 1)

    with open(args.input, 'r') as fin:
        ns = 0
        while ns < nskip:
            line = fin.readline()
            if len(line.split()) > 0:
                ns += 1

        nk = 0
        while nk < nkeep:
            line = fin.readline()
            if len(line.split()) > 0:
                fout.write(line)
                nk += 1

    if not args.output is None:
        fout.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input with stacked gene alignments in phylip format (not interleaved)",
                        required=True)

    parser.add_argument("-n", "--gnum", type=int,
                        help="Gene number",
                        required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output file",
                        required=False)

    main(parser.parse_args())

