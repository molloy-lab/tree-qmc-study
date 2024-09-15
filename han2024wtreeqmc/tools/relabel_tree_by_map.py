import argparse
import treeswift
import sys


def get_newick_string(tree):
    """
    Creates a newick string with  branch lengths in decimal notation.
    This function was taken from Phybase.py, see
    github.com/ngcrawford/CloudForest/blob/master/cloudforest/phybase.py

    Parameters
    ----------
    tree : treeswift tree

    Returns
    -------
    new_tree : string
               newick representation of tree
    """
    from decimal import Decimal

    str_newick_tree = tree.newick()

    colon_s = 0
    comma_back_paren_s = 0
    num = ''
    new_tree = ''

    for count, char in enumerate(str_newick_tree):
        if char == ':':
            colon_s = count
            continue
        if char in (')', ','):
            comma_back_paren_s = 1
            try:
                tmp = float(num)
                num = '%1.12f' % Decimal(num)
                new_tree += ":" + num
            except ValueError:
                pass
            colon_s = 0
            num = ''
        if colon_s != 0:
            num = num + char
        if colon_s == 0:
            new_tree += char

    new_tree = new_tree.strip('\'').strip('\"').strip('\'')
    if new_tree[-1] != ';':
        new_tree = new_tree + ';'
    new_tree = new_tree + '\n'

    return new_tree


def read_taxon_name_map(fname):
    nmap = {}
    with open(fname, 'r') as f:
        for line in f:
            labels = line.split()

            if len(labels) != 2:
                sys.exit("Bad name map!\n")

            [slab, elab] = labels
            slab = slab.strip()
            elab = elab.strip()
            #elab = elab.replace("UUU", "_")

            nmap[slab] = elab
    return nmap


def relabel_taxa(tree, nmap):
    for l in tree.traverse_leaves():
        l.label = nmap[l.label]


def main(args):
    # Read taxon name map
    nmap = read_taxon_name_map(args.taxon_name_map)

    with open(args.input, 'r') as fi, \
         open(args.output, 'w') as fo:
        for line in fi:
            tmp = line.strip()
            if tmp.find(';') > -1:
                tree = treeswift.read_tree(tmp, schema='newick')
                relabel_taxa(tree, nmap)

                if args.removelt4:
                    taxa = [x for x in tree.traverse_leaves()]
                    if len(taxa) > 3:
                        fo.write(get_newick_string(tree))
                else:
                    fo.write(get_newick_string(tree))
            else:
                sys.stdou.write("WARNING: Skipping line %s because may not be a tree!\n" % line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input file", required=True)
    parser.add_argument("-x", "--taxon_name_map", type=str,
                        help="Taxon name map file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output file", required=True)
    parser.add_argument("--removelt4", action="store_true")

    main(parser.parse_args())


