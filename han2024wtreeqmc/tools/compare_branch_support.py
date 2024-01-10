import argparse
import numpy
import treeswift
import sys


def get_leaf_set(tree):
    return set([leaf.label for leaf in tree.traverse_leaves()])


def get_astral_support_value(label):
    if label[0] != '[':
        sys.exit("Unable to parse!")
    else:
        label = label[1:-1].split(';')
        pp1 = numpy.nan
        q1 = numpy.nan
        for x in label:
            [key, val] = x.split('=')
            if key == "pp1":
                #pp1 = float(val)
                pp1 = val
            elif key == "q1":
                #q1 = float(val)
                q1 = val
        return [pp1, q1]
    sys.exit("Unable to parse!")


def get_branch_set(tree):
    bset = {}

    leaves = sorted([leaf.label for leaf in tree.traverse_leaves()])
    nl = len(leaves)
    root = leaves[0]  # Randomly pick root

    leaf_set = set(leaves)

    for node in tree.traverse_postorder():
        if node.is_leaf():
            pass
        elif node.is_root():
            pass
        else:
            clade = sorted([leaf.label for leaf in node.traverse_leaves()])
            if clade[0] == root:
                clade = sorted(list(leaf_set.difference(set(clade))))

            nb = len(clade)
            if (nb > 1) and (nb < nl - 1):
                cstr = ','.join(clade)
                bset[cstr] = {}
                bset[cstr]["supp"] = get_astral_support_value(node.label)
                bset[cstr]["elen"] = node.edge_length
    return bset


def process_false_negatives(ttre_bset, etre1_bset, etre2_bset):
    fn_tt = [] # first holds pp1, second holds q1
    fn_e1 = []
    fn_e2 = []
    count = 0

    for br in ttre_bset:
        # Process tree 1
        e1_miss = True
        try:
            e1_supp = etre1_bset[br]["supp"]
            e1_miss = False
        except KeyError:
            e1_supp = ["NA", "NA"]

        # Process tree 2
        e2_supp = ["NA", "NA"]
        e2_miss = True
        try:
            e2_supp= etre2_bset[br]["supp"]
            e2_miss = False
        except KeyError:
            pass

        # Add false negative info
        if e1_miss or e2_miss:
            count += 1
            # Save information if FP in et1 OR et2
            fn_tt.append(ttre_bset[br]["supp"])
            fn_e1.append(e1_supp)
            fn_e2.append(e2_supp)

    if not count:
        return "\"\",\"\",\"\",\"\",\"\",\"\"" 

    fn_tt = numpy.array(fn_tt, dtype=object)
    fn_e1 = numpy.array(fn_e1, dtype=object)
    fn_e2 = numpy.array(fn_e2, dtype=object)

    data = []
    for i in range(fn_tt.shape[1]):
        data.append('\"' + ','.join(fn_tt[:,i]) + '\"')
        data.append('\"' + ','.join(fn_e1[:,i]) + '\"')
        data.append('\"' + ','.join(fn_e2[:,i]) + '\"')

    #return '\n'.join(data)
    return ','.join(data)


def process_false_positives(ttre_bset, etre1_bset, etre2_bset):
    fp_e1 = []
    fp_e2 = []
    count = 0

    brs = set(etre1_bset.keys()).union(etre2_bset.keys())

    for br in brs:
        try:
            # True positive branch so do nothing
            tmp = ttre_bset[br]
        except KeyError:
            # False positive branch so grab support
            count += 1
            try:
                fp_e1.append(etre1_bset[br]["supp"])
            except KeyError:
                fp_e1.append(["NA", "NA"])

            try:
                fp_e2.append(etre2_bset[br]["supp"])
            except KeyError:
                fp_e2.append(["NA", "NA"])

    if not count:
        return "\"\",\"\",\"\",\"\""

    fp_e1 = numpy.array(fp_e1, dtype=object)
    fp_e2 = numpy.array(fp_e2, dtype=object)

    data = []
    for i in range(fp_e1.shape[1]):
        data.append('\"' + ','.join(fp_e1[:,i]) + '\"')
        data.append('\"' + ','.join(fp_e2[:,i]) + '\"')

    #return '\n'.join(data)
    return ','.join(data)


def process_true_positives(ttre_bset, etre1_bset, etre2_bset):
    tp_tt = [] # first holds pp1, second holds q1
    tp_e1 = []
    tp_e2 = []
    count = 0

    for br in ttre_bset:
        # Process tree 1
        e1_miss = True
        try:
            e1_supp = etre1_bset[br]["supp"]
            e1_miss = False
        except KeyError:
            e1_supp = ["NA", "NA"]

        # Process tree 2
        e2_supp = ["NA", "NA"]
        e2_miss = True
        try:
            e2_supp= etre2_bset[br]["supp"]
            e2_miss = False
        except KeyError:
            pass

        # Add false negative info
        if (not e1_miss) or (not e2_miss):
            count += 1
            # Save information if TP in et1 OR et2
            tp_tt.append(ttre_bset[br]["supp"])
            tp_e1.append(e1_supp)
            tp_e2.append(e2_supp)

    if not count:
        return "\"\",\"\",\"\",\"\",\"\",\"\"" 

    tp_tt = numpy.array(tp_tt, dtype=object)
    tp_e1 = numpy.array(tp_e1, dtype=object)
    tp_e2 = numpy.array(tp_e2, dtype=object)

    data = []
    for i in range(tp_tt.shape[1]):
        data.append('\"' + ','.join(tp_tt[:,i]) + '\"')
        data.append('\"' + ','.join(tp_e1[:,i]) + '\"')
        data.append('\"' + ','.join(tp_e2[:,i]) + '\"')

    #return '\n'.join(data)
    return ','.join(data)


def compare_branch_support(ttre, etre1, etre2):
    # Suppress unifurcations
    ttre.suppress_unifurcations()
    etre2.suppress_unifurcations()
    etre2.suppress_unifurcations()

    # Restrict trees to same leaf set
    ttre_lset = get_leaf_set(ttre)
    etre1_lset = get_leaf_set(etre1)
    etre2_lset = get_leaf_set(etre2)

    keepset = ttre_lset.intersection(etre1_lset.intersection(etre2_lset))

    ttre = ttre.extract_tree_with(keepset)
    etre1 = etre1.extract_tree_with(keepset)
    etre2 = etre2.extract_tree_with(keepset)

    # Get branch set for trees
    ttre_bset = get_branch_set(ttre)
    etre1_bset = get_branch_set(etre1)
    etre2_bset = get_branch_set(etre2)

    fn_data = process_false_negatives(ttre_bset, etre1_bset, etre2_bset)
    fp_data = process_false_positives(ttre_bset, etre1_bset, etre2_bset)
    tp_data = process_true_positives(ttre_bset, etre1_bset, etre2_bset)

    return [fn_data, fp_data, tp_data]


def main(args):
    # Get prefix of CSV line
    if args.prefix is None:
        prefix = ""
    else:
        prefix = str(args.prefix + ',')

    # Read trees
    with open(args.truetree, 'r') as fin:
        line = fin.readline()
        ttre = treeswift.read_tree(line, "newick")

    # Read trees
    with open(args.estitree1, 'r') as fin:
        line = fin.readline()
        etre1 = treeswift.read_tree(line, "newick")

    with open(args.estitree2, 'r') as fin:
        line = fin.readline()
        etre2 = treeswift.read_tree(line, "newick")

    # Compare two trees
    [fn, fp, tp] = compare_branch_support(ttre, etre1, etre2)

    # Write CSV to standard output
    # PREFIX,
    # FN_TT_PP, FN_E1_PP, FN_E2_PP
    # FN_TT_QS, FN_E1_QS, FN_E2_QS
    # FP_E1_PP, FP_E2_PP
    # FP_E1_QS, FP_E2_QS
    # TP_TT_PP, TP_E1_PP, TP_E2_PP
    # TP_TT_QS, TP_E1_QS, TP_E2_QS

    #sys.stdout.write("%s\n%s\n\n%s\n\n%s\n" % (prefix, fn, fp, tp))
    sys.stdout.write("%s%s,%s,%s\n" % (prefix, fn, fp, tp))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t", "--truetree", type=str,
                        help="Input containing true tree with ASTRAL -u 2 branch support",
                        required=True)

    parser.add_argument("-e1", "--estitree1", type=str,
                        help="Input containing estimated tree with ASTRAL -u 2 branch support",
                        required=True)

    parser.add_argument("-e2", "--estitree2", type=str,
                        help="Input containing estimated tree with ASTRAL -u 2 branch support",
                        required=True)

    parser.add_argument("-p", "--prefix", type=str,
                        help="Add prefix to CSV output",
                        required=False)

    main(parser.parse_args())
