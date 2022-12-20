"""
This file is used prepare quartets for wQMC.
Copyright (c) 2022 Erin K. Molloy
All rights reserved.
License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import argparse
import dendropy
import os
import sys




def main(args):
    # Read tree
    taxa = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(path=args.tree,
                             schema='newick',
                             rooting='force-unrooted',
                             taxon_namespace=taxa)

    # Read taxa
    xlst = []
    with open(args.taxa, 'r') as fi:
        for line in fi:
            xlst.append(line.strip())

    # Relabel taxa
    for x in taxa:
        ind = int(x.label)
        x.label = xlst[ind]

    # Write relabeled tree
    with open(args.output, 'w') as fo:
        fo.write(tree.as_string(schema="newick")[5:])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", type=str,
                        help="Input file containing tree",
                        required=True)
    parser.add_argument("-x", "--taxa", type=str,
                        help="Input file containing list of taxa (one per line)",
                        required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file",
                        required=True)

    main(parser.parse_args())

