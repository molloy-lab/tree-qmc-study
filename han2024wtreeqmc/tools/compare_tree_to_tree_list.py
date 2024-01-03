import argparse
from compare_two_trees import compare_two_trees
import treeswift
import sys


def main(args):
    if args.prefix is None:
        prefix = ""
    else:
        prefix = str(args.prefix + ",")

    # Read tree 1
    with open(args.tree, 'r') as fin:
        line = fin.readline()
        while len(line.strip()) == 0:
            line = fin.readline()
        tree1 = treeswift.read_tree(line, "newick")

    # Count number of genes / lines
    ngen = 0
    with open(args.treelist, 'r') as fin:
        for line in fin:
            if len(line.strip()) > 0:
                ngen += 1

    with open(args.treelist, 'r') as fin:
        for gene in range(ngen):
            # Read tree 2 from list
            line = fin.readline()
            while len(line.strip()) == 0:
                line = fin.readline()
            tree2 = treeswift.read_tree(line, "newick")

            # Compare two trees
            info = compare_two_trees(tree1, tree2)
            [n1, n2, nl, i1, i2, fn, fp] = info

            sys.stdout.write("%s%d,%d,%d,%d,%d,%d,%d,%d\n" \
                % (prefix, gene + 1, n1, n2, nl, i1, i2, fn, fp))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", type=str,
                        help="Input file containing tree list 1", required=True)
    parser.add_argument("-l", "--treelist", type=str,
                        help="Input file containing tree list 2", required=True)
    parser.add_argument("-p", "--prefix", type=str,
                        help="Append prefix to each row of CSV", required=False)

    main(parser.parse_args())

