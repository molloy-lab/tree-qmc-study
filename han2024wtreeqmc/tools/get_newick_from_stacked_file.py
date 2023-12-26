import argparse
import sys
import treeswift


def main(args):
    if args.output is None:
        fout = sys.stdout
    else:
        fout = open(args.output, 'w')

    nskip = args.gnum - 1

    with open(args.input, 'r') as fin:
        for ns in range(nskip):
            line = fin.readline()

        line = fin.readline()

    tree = treeswift.read_tree(line, schema='newick')

    for node in tree.traverse_preorder():
        node.edge_length = None
        if not node.is_leaf():
            node.label = None

    fout.write(tree.newick())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input with stacked gene trees in newick format",
                        required=True)

    parser.add_argument("-n", "--gnum", type=int,
                        help="Gene number",
                        required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output file",
                        required=False)

    main(parser.parse_args())

