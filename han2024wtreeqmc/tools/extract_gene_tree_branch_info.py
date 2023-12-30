import argparse
import treeswift
from tswrapper import newick
import sys


def get_leaf_set_string(tree):
    labels = []
    for leaf in tree.traverse_leaves():
        labels.append(leaf.label)
    labels = sorted(labels)
    leaf_set_str = ','.join(labels)
    return leaf_set_str


def get_branch_info(tree):
    info = {}

    leaves = sorted([leaf.label for leaf in tree.traverse_leaves()])
    root = leaves[0]

    leaf_set = set(leaves)

    for node in tree.traverse_postorder():
        if node.is_leaf():
            pass
        elif node.is_root():
            pass
        else:
            leaves_below = sorted([leaf.label for leaf in node.traverse_leaves()])
            if leaves_below[0] == root:
                leaf_set_below = leaf_set.difference(set(leaves_below))
                leaves_below = sorted(list(leaf_set_below))

            branch = ','.join(leaves_below)

            label = node.label
            if label is None:
                # FastTree-2 outputs polytomies for identical sequences,
                # IQTree randomly refines polytomies and sets support to None
                label = "NA"
                node.label = '0.000'
                node.edge_length = 0.0

            info[branch] = [str(node.edge_length), label]

    return info


def main(args):
    if args.prefix is None:
        prefix = ""
    else:
        prefix = str(args.prefix + ',')

    ngen = 0
    with open(args.truetrees, 'r') as fin_true:
        for line in fin_true:
            if len(line.strip()) > 0:
                ngen += 1
    
    fin_true = open(args.truetrees, 'r')
    fin_esti = open(args.estitrees, 'r')
    if not args.output is None:
        fout = open(args.output, 'w')

    for gene in range(ngen):
        # Read true tree
        line = fin_true.readline()
        while len(line.strip()) == 0:
            line = fin_true.readline()
        true_tree = treeswift.read_tree(line, "newick")

        # Read estimated tree
        line = fin_esti.readline()
        while len(line.strip()) == 0:
            line = fin_esti.readline()
        esti_tree = treeswift.read_tree(line, "newick")

        # Check trees are on the same leaf sets
        if get_leaf_set_string(true_tree) != \
           get_leaf_set_string(esti_tree):
            sys.exit("True and estimated trees are on different leaf sets!")

        # Get branch info
        true_info = get_branch_info(true_tree)
        esti_info = get_branch_info(esti_tree)

        true_brln_list = []
        esti_brln_list = []
        esti_supp_list = []

        # Write branch info
        branches = sorted([key for key in esti_info.keys()])

        for branch in branches:
            [esti_brln, esti_supp] = esti_info[branch]

            try:
                [true_brln, true_supp] = true_info[branch]
            except KeyError:
                true_brln = "NA"

            true_brln_list.append(true_brln)
            esti_brln_list.append(esti_brln)
            esti_supp_list.append(esti_supp)

        true_brln_str = ','.join(true_brln_list)
        esti_brln_str = ','.join(esti_brln_list)
        esti_supp_str = ','.join(esti_supp_list)

        sys.stdout.write("%s%d,\"%s\",\"%s\",\"%s\"\n" \
            % (prefix, gene + 1, true_brln_str, esti_brln_str, esti_supp_str))

        # Write cleaned up output tree
        if not args.output is None:
            fout.write(newick(esti_tree))
            fout.write('\n')

    fin_true.close()
    fin_esti.close()
    if not args.output is None:
        fout.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--truetrees", type=str,
                        help="True gene trees",
                        required=True)
    parser.add_argument("-e", "--estitrees", type=str,
                        help="Estimated gene trees",
                        required=True)
    parser.add_argument("-p", "--prefix", type=str,
                        help="Append prefix to each row of CSV",
                        required=False)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file containing fixed gene trees",
                        required=False)
    
    main(parser.parse_args())
