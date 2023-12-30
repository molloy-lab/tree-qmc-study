import argparse
from compare_two_trees import compare_two_trees
import treeswift
import sys


def main(args):
    # Get prefix for CSV line
    if args.prefix is None:
        prefix = ""
    else:
        prefix = str(args.prefix + ',')

    # Count number of genes / lines
    ngen = 0
    with open(args.treelist1, 'r') as fin1:
        for line in fin1:
            if len(line.strip()) > 0:
                ngen += 1

    fin1 = open(args.treelist1, 'r')
    fin2 = open(args.treelist2, 'r')

    for gene in range(ngen):
        # Read tree from list 1
        line1 = fin1.readline()
        while len(line1.strip()) == 0:
            line1 = fin1.readline()
        tree1 = treeswift.read_tree(line1, "newick")

        # Read tree from list 2
        line2 = fin2.readline()
        while len(line2.strip()) == 0:
            line2 = fin2.readline()
        tree2 = treeswift.read_tree(line2, "newick")

        # Compare two trees
        info = compare_two_trees(tree1, tree2)
        [n1, n2, nl, i1, i2, fn, fp] = info

        #if fn > 0: 
        #    sys.stderr.write("Found %d false negatives" % fn)

        #if fp > 0: 
        #    sys.stderr.write("Found %d false positives" % fp)

        sys.stdout.write("%s%d,%d,%d,%d,%d,%d,%d,%d\n" \
            % (prefix, gene + 1, n1, n2, nl, i1, i2, fn, fp))

    fin1.close()
    fin2.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-l1", "--treelist1", type=str,
                        help="Input containing tree list 1",
                        required=True)
    parser.add_argument("-l2", "--treelist2", type=str,
                        help="Input containing tree list 2",
                        required=True)
    parser.add_argument("-p", "--prefix", type=str,
                        help="Add prefix to each row of CSV",
                        required=False)

    main(parser.parse_args())
