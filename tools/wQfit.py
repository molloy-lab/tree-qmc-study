"""
This file was written by Yunheng Han on Feb 13, 2023.
This file was edited by Erin Molloy on Feb 13, 2023
"""
import argparse
import sys

def join_(a, b, c, d):
    q0 = a + ',' + b if a < b else b + ',' + a
    q1 = c + ',' + d if c < d else d + ',' + c
    q = '((' + q0 + '),(' + q1 +'));' if q0 < q1 else '((' + q1 + '),(' + q0 +'));'
    return q

def join(q):
    return join_(q[0], q[1], q[2], q[3])

def split(q):
    s = []
    for st in q.split(','):
        s.append(st.replace('(', '').replace(')', '').replace(';', ''))
    return s

def merge(l):
    s = ''
    for label in l:
        if s == '':
            s = label
        else:
            s = s + ',' + label
    return s

def combine(l):
    s = []
    for li in l:
        s.extend(li)
    s.sort()
    return s

class Tree:
    class Node:
        def __init__(self, label = ''):
            self.label = label
            self.children = []
            self.parent = None
            self.depth = -1
            self.info = ''

    def __init__(self, newick):
        i = len(newick) - 1
        while newick[i] != ')': i -= 1
        self.root = self.build_tree(newick[:i+1])
        self.get_depth(self.root, 0)

    def get_depth(self, node, depth):
        node.depth = depth
        for child in node.children:
            self.get_depth(child, depth + 1)

    def build_tree(self, newick):
        if len(newick) == 0 or newick[0] != '(':
            s = newick.split(':')
            root = self.Node(s[0])
            if len(s) > 1: root.info = s[1]
        else:
            root = self.Node()
            j = 0
            k = 1
            for i in range(len(newick)):
                if newick[i] == '(': j += 1
                if newick[i] == ')': j -= 1
                if newick[i] == ',' and j == 1:
                    root.children.append(self.build_tree(newick[k : i]))
                    k = i + 1
            i = len(newick) - 1
            while newick[i] != ')': i -= 1
            root.info = newick[i + 1:len(newick)]
            root.children.append(self.build_tree(newick[k : i]))
            for child in root.children:
                child.parent = root
        return root

    def lca(self, a, b):
        ca = a
        cb = b
        while ca.depth > cb.depth: ca = ca.parent
        while cb.depth > ca.depth: cb = cb.parent
        while ca != cb:
            ca = ca.parent
            cb = cb.parent
        return ca
    
    def get_quartet(self, leaves):
        lcas = []
        for i in range(4):
            for j in range(i + 1, 4):
                lcas.append([self.lca(leaves[i], leaves[j]), i, j])
        def skey(a): return - a[0].depth
        lcas.sort(key=skey)
        if lcas[0][0] == lcas[1][0]: return ''
        a = leaves[lcas[0][1]].label
        b = leaves[lcas[0][2]].label
        c = ''
        d = ''
        for i in range(4):
            if not leaves[i].label in [a, b]:
                c = leaves[i].label
        for i in range(4):
            if not leaves[i].label in [a, b, c]:
                d = leaves[i].label
        return join([a, b, c, d])
    
    def get_quartets(self):
        s = {}
        leaves = self.get_leaves(self.root)
        for a in range(len(leaves)):
            for b in range(a + 1, len(leaves)):
                for c in range(b + 1, len(leaves)):
                    for d in range(c + 1, len(leaves)):
                        q = self.get_quartet([leaves[a], leaves[b], leaves[c], leaves[d]])
                        if q != '':
                            if not q in s: s[q] = 0.0
                            s[q] += 1.0
        return s

    def bipartition(self, node):
        labels = self.get_labels()
        subset0 = self.get_labels_(node)
        subset0.sort()
        labelset = set(subset0)
        subset1 = []
        for label in labels:
            if not label in labelset:
                subset1.append(label)
        subset1.sort()
        return [subset0, subset1]

    def quardpartition(self, node):
        if node.parent == None:
            return []
        else:
            labels = self.get_labels()
            labelset = set()
            subsets0 = []
            for child in node.parent.children:
                if child != node:
                    subset = self.get_labels_(child)
                    for label in subset:
                        labelset.add(label)
                    if len(subset) > 0: 
                        subsets0.append(subset)
            subsets1 = []
            for child in node.children:
                subset = self.get_labels_(child)
                for label in subset:
                    labelset.add(label)
                if len(subset) > 0:
                    subsets1.append(subset)
            subset = []
            for label in labels:
                if not label in labelset:
                    subset.append(label)
            if len(subset) > 0:
                subsets0.append(subset)
            return [subsets0, subsets1]

    def get_quardpartitions(self):
        nodes = self.get_nodes(self.root)
        s = []
        for node in nodes:
            qup = self.quardpartition(node)
            if len(qup) == 0: continue
            s.append([qup, node])
        return s

    def get_bipartitions(self):
        nodes = self.get_nodes(self.root)
        s = []
        for node in nodes:
            if node == self.root: continue
            if len(self.root.children) == 2 and node == self.root.children[0]: continue
            if len(node.children) == 0: continue
            s.append([self.bipartition(node), node])
        return s

    def get_nodes(self, root):
        s = [root]
        if len(root.children) != 0:
            for child in root.children:
                s.extend(self.get_nodes(child))
        return s
    
    def get_leaves(self, root):
        s = []
        if len(root.children) != 0:
            for child in root.children:
                s.extend(self.get_leaves(child))
        else:
            s.append(root)
        return s

    def to_string(self, full = True):
        return self.display_tree(self.root, full) + ';'

    def display_tree(self, root, full):
        s = ''
        info = root.info if full else ''
        if len(root.children) == 0:
            s = root.label
            if info != '': s += ':'
        else:
            for child in root.children:
                s = s + ',' + self.display_tree(child, full)
            s = s[1:]
            s = '(' + s + ')'
        return s + info

    def get_labels(self):
        return self.get_labels_(self.root)

    def get_labels_(self, root):
        if len(root.children) == 0:
            return [root.label]
        else:
            s = []
            for child in root.children:
                s.extend(self.get_labels_(child))
            return s


def compute_wqfit(quartets_t1, quartets_t2, qweights, taxa):
    s12 = s11 = s22 = 0.0

    ntax = len(taxa)

    for a in range(ntax):
        for b in range(a + 1, ntax):
            for c in range(b + 1, ntax):
                for d in range(c + 1, ntax):
                    q0 = join([taxa[a], taxa[b], taxa[c], taxa[d]])
                    q1 = join([taxa[a], taxa[c], taxa[b], taxa[d]])
                    q2 = join([taxa[a], taxa[d], taxa[b], taxa[c]])
                    
                    qws = [qweights[q] if q in qweights else 0.0 for q in [q0, q1, q2]]
                    tot = sum(qws)

                    normw = [qw / tot for qw in qws]
                    
                    i1 = i2 = -1
                    w1 = w2 = 0
                    
                    if q0 in quartets_t1: 
                        i1 = 0
                        w1 = normw[0]
                    elif q1 in quartets_t1: 
                        i1 = 1
                        w1 = normw[1]
                    elif q2 in quartets_t1: 
                        i1 = 2
                        w1 = normw[2]
                    
                    if q0 in quartets_t2: 
                        i2 = 0
                        w2 = normw[0]
                    elif q1 in quartets_t2: 
                        i2 = 1
                        w2 = normw[1]
                    elif q2 in quartets_t2: 
                        i2 = 2
                        w2 = normw[2]
                    
                    d12 = 2 if i1 == i2 else -1
                    
                    s12 += d12 * w1 * w2
                    s11 += 2 * w1 * w1
                    s22 += 2 * w2 * w2
    
    #print(s12, s11, s22)

    wqfit = (2.0 * s12) / float(s11 + s22)

    return wqfit


def main(args):
    # Read tree 1 (typically true tree)
    with open(args.tree1, 'r') as f:
        newick = f.readlines()[0]
    tree1 = Tree(newick)
    quartets_t1 = tree1.get_quartets()

    # Read tree 2 (typically estimated tree)
    with open(args.tree2, 'r') as f:
        newick = f.readlines()[0]
    tree2 = Tree(newick)
    quartets_t2 = tree2.get_quartets()

    # Find shared taxon set
    tax1 = [node.label for node in tree1.get_leaves(tree1.root)]
    tax2 = [node.label for node in tree2.get_leaves(tree2.root)]
    taxa = list(set().union(tax1, tax2))

    # Read weighted quartets
    qweights = {}
    with open(args.qweights, 'r') as f:
        for line in f.readlines():
            s = line.split()
            qweights[join(split(s[0]))] = int(s[1])

    # Compute wQfit score defined in Avni et al. (2015)
    wqfit = compute_wqfit(quartets_t1, quartets_t2, qweights, taxa)

    sys.stdout.write("%1.12f" % wqfit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t1", "--tree1", type=str,
                        help="Input file containing tree 1",
                        required=True)
    parser.add_argument("-t2", "--tree2", type=str,
                        help="Input file containing tree 2",
                        required=True)
    parser.add_argument("-qw", "--qweights", type=str,
                        help="Input file containing weighted quartets in wQFM format",
                        required=True)

    main(parser.parse_args())
