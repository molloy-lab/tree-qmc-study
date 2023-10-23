"""
This file is used prepare quartets for wQMC.
Copyright (c) 2022 Erin K. Molloy
All rights reserved.
License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import argparse
import os
import sys


def add_taxon(xmap, xlst, old_x):
    tot = len(xlst)
    try:
        new_x = xmap[old_x]
    except KeyError:
        xmap[old_x] = tot
        xlst.append(old_x)
        new_x = xmap[old_x]
    return new_x


def write_taxon_map(xlst, outf):
    with open(outf, 'w') as fo:
        for x in xlst:
            fo.write(str("%s\n" % x))


def main(args):
    xmap = {}
    xlst = []

    outq = args.output + ".qrt"
    outx = args.output + "_taxon_map.txt"

    with open(outq, 'w') as fo, \
         open(args.input, 'r') as fi:

        for line in fi:
            [old_q, weight] = line.split('; ')
           
            [old_qab, old_qcd] = old_q.split('),(')
            [old_qa, old_qb] = old_qab.split(',')
            [old_qc, old_qd] = old_qcd.split(',')

            old_qa = old_qa.replace('((', '')
            old_qd = old_qd.replace('))', '')

            new_qa = add_taxon(xmap, xlst, old_qa)
            new_qb = add_taxon(xmap, xlst, old_qb)
            new_qc = add_taxon(xmap, xlst, old_qc)
            new_qd = add_taxon(xmap, xlst, old_qd)

            if new_qa < new_qb:
                left = new_qa
                new_qab = str("%d,%d" % (new_qa, new_qb))
            else:
                left = new_qb
                new_qab = str("%d,%d" % (new_qb, new_qa))

            if new_qc < new_qd:
                rght = new_qc
                new_qcd = str("%d,%d" % (new_qc, new_qd))
            else:
                rght = new_qd
                new_qcd = str("%d,%d" % (new_qd, new_qc))

            if left < rght:
                new_q = str("%s|%s:%s" % (new_qab, new_qcd, weight))
            else:
                new_q = str("%s|%s:%s" % (new_qcd, new_qab, weight))

            fo.write("%s" % new_q)

    write_taxon_map(xlst, outx)    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file containing quartets",
                        required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output prefix",
                        required=True)

    main(parser.parse_args())

