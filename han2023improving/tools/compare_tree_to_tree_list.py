import argparse
import dendropy
from compare_two_trees import compare_trees


def main(args):
    if args.prefix is None:
        p = ""
    else:
        p = str(args.prefix + ",")

    with open(args.output, 'a+') as fo, \
         open(args.tree, 'r') as f1, \
         open(args.treelist, 'r') as f2:

        l1 = f1.readline();

        i = 1
        for l2 in f2:
            taxa = dendropy.TaxonNamespace()
            tre1 = dendropy.Tree.get(string=l1,
                                     schema='newick',
                                     rooting='force-unrooted',
                                     taxon_namespace=taxa)

            tre2 = dendropy.Tree.get(string=l2,
                                     schema='newick',
                                     rooting='force-unrooted',
                                     taxon_namespace=taxa)

            [nl, ei1, ei2, fn, fp, rf] = compare_trees(tre1, tre2)
            if rf == "NA":
                fo.write('%s%d,%d,%d,%d,%s,%s,%s\n' % \
                         (p, i, nl, ei1, ei2, fn, fp, rf))
            else:
                fo.write('%s%d,%d,%d,%d,%d,%d,%1.6f\n' % \
                         (p, i, nl, ei1, ei2, fn, fp, rf))

            i = i + 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--tree", type=str,
                        help="Input file containing tree list 1", required=True)
    parser.add_argument("-l", "--treelist", type=str,
                        help="Input file containing tree list 2", required=True)
    parser.add_argument("-p", "--prefix", type=str,
                        help="Append prefix to each row of CSV", required=False)
    parser.add_argument("-o", "--output", type=str,
                        help="Output CSV file", required=True)

    main(parser.parse_args())

