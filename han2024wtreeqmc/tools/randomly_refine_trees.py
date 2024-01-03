import argparse
import tswrapper 
import treeswift
import sys


def main(args):
    # Count number of trees
    ntree = 0
    with open(args.input, 'r') as fin:
        for line in fin:
            if len(line.strip()) > 0:
                ntree += 1

    fin = open(args.input, 'r')

    for i in range(ntree):
        # Read tree from list
        line = fin.readline()
        while len(line.strip()) == 0:
            line = fin.readline()
        tree = treeswift.read_tree(line, "newick")

        # Randomly refine and write to newick string
        tswrapper.randomly_refine(tree, args.seed)
        sys.stdout.write("%s\n" % tswrapper.newick(tree))

    fin.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input containing tree list",
                        required=True)

    parser.add_argument("-s", "--seed", type=int, default=12345,
                        help="Seed for random number generator",
                        required=False)

    main(parser.parse_args())

