"""
Grab branch info

Written by Erin K. Molloy (github: ekmolloy)
"""
import argparse
import treeswift
import sys


def parse_branch_label(label):
    data = {}

    try:
        xs = label.split(';')
    except TypeError:
        return data

    for i, x in enumerate(xs):
        try:
            key, value = x.split('=')

            if i == 0:
                key = key.replace("\'[", '')
            if i == len(xs) - 1:
                value = value.replace("]\'", '')

            data[key] = float(value)
        except ValueError:
            return data

    return data


def main(args):
    with open(args.input, 'r') as fi:
        line = fi.readline()
        tree = treeswift.read_tree(line, schema='newick')

    taxa = set([x.label for x in tree.traverse_leaves()])

    for node in tree.traverse_preorder():
        bipA = set([x.label for x in node.traverse_leaves()])
        bipB = taxa.difference(bipA)

        if (len(bipA) > 1) and (len(bipB) > 1):
            bipA = ','.join(sorted(list(bipA)))
            bipB = ','.join(sorted(list(bipB)))
            tmp = sorted([bipA, bipB])

            elen = node.get_edge_length()
            if elen is not None:
                brln = str("%f" % float(elen))
            else:
                brln = "NA"

            if node.label is not None:
                label = parse_branch_label(node.label)
                if len(label) > 0:
                    q1 = str("%f" % label["q1"])
                    q2 = str("%f" % label["q2"])
                    q3 = str("%f" % label["q3"])
                    qc = str("%f" % label["QC"])
                    en = str("%f" % label["EN"])
                    pp = str("%f" % label["pp1"])
                    info = ','.join([pp, en, q1, q2, q3, qc])
                else:
                    info = node.label + ",NA,NA,NA,NA,NA"

            else:
                info = "NA,NA,NA,NA,NA,NA"

            sys.stdout.write("%s,\"%s\",\"%s\",%s,%s\n" %
                             (args.prefix, tmp[0], tmp[1], brln, info))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input newick string",
                        required=True)

    parser.add_argument("-p", "--prefix", type=str,
                        help="Prefix",
                        required=True)

    main(parser.parse_args())

