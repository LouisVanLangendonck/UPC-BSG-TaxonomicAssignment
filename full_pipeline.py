import argparse
import networkx as nx
import os
import time

def get_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-idump", "--inpdump",
                        help="Location of nodes.dmp file",
                        type=str, default="./nodes.dmp")
    p.add_argument("-iread", "--inpread",
                        help="Location of nodes.dmp file",
                        type=str, default="./sample.inp")
    return p

def find_pruning_node(tree, leaf):
    if leaf != 1:
        for pred in tree.predecessors(leaf):
            if tree.out_degree(leaf) > 1:
                return []
            else:
                return [leaf] + find_pruning_node(tree, pred)
    else:
        return []

def prune_tree(tree):
    to_delete = set()
    for node in tree.nodes():
        if tree.out_degree(node) == 0:
            if tree.nodes[node]['rank'] != 'species':
                for node_to_delete in find_pruning_node(tree=tree, leaf=node):
                    to_delete.add(node_to_delete)
    for node in to_delete:
        tree.remove_node(node)
    return tree

def make_tree(input_file,species_only_leaves = False):
    G = nx.DiGraph()
    corr_ranks = ['superkingdom','phylum','class','order','family','genus','species']
    nodes_to_delete = set()

    with open(input_file, 'r') as node_file:
        print('Reading input file & Building Tree...')
        for nodes in node_file.read().split("\t|\n")[:-1]: #Iter Nodes
            tax_id, parent_tax_id, rank_str = nodes.split("\t|\t")[0:3]
            tax_id = int(tax_id)
            parent_tax_id = int(parent_tax_id)
            if tax_id != 1:
                if rank_str not in corr_ranks:
                    G.add_node(tax_id)
                    nodes_to_delete.add(tax_id)
                else:
                    G.add_node(tax_id, rank = rank_str)
                G.add_edge(parent_tax_id, tax_id)

    for node in nodes_to_delete:
        for preds in G.predecessors(node):
            for success in G.successors(node):
                G.add_edge(preds, success)
        G.remove_node(node)

    if species_only_leaves:
        print('Pruning Non-Species leaves...')
        final_res = prune_tree(G)
    print("Tree conversion completed!")
    return final_res

def find_skeleton(input_graph, reads):
    order = ['genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'root']
    res_graph = nx.DiGraph()
    input_graph.nodes[1]['rank'] = 'root'

    nodes_to_analyze = reads
    for level in order:
        predecessor_set = set()
        for node in nodes_to_analyze:
            for pred in input_graph.predecessors(node):
                if input_graph.nodes[pred]['rank'] == level:
                    res_graph.add_edge(pred, node)
                    predecessor_set.add(pred)
                else:
                    predecessor_set.add(node)
        nodes_to_analyze = list(predecessor_set)

    to_delete = set()
    for node in res_graph.nodes():
        if res_graph.out_degree(node) == 1 and res_graph.in_degree(node) == 1:
            to_delete.add(node)

    for node in list(to_delete):
        for pred in res_graph.predecessors(node):
            for suc in res_graph.successors(node):
                res_graph.add_edge(pred, suc)
        res_graph.remove_node(node)
    return res_graph

def nr_of_leaf_succ(graph, node):
    children = 0
    for suc in graph.successors(node):
        children += nr_of_leaf_succ(graph, suc)
    graph.nodes[node]['nr_of_leaf_succ'] = children
    if children == 0:
        return 1
    return children

def best_F(major_tree, skel_tree, iteration):
    res = {}
    for node in skel_tree.nodes():
        if skel_tree.in_degree(node) == 0:
            root = node
            nr_of_leaf_succ(skel_tree, root)
    root_skel_succ = skel_tree.nodes[root]['nr_of_leaf_succ']
    for node in skel_tree.nodes():
        if skel_tree.out_degree(node) != 0:
            node_skel_succ = skel_tree.nodes[node]['nr_of_leaf_succ']
            TP = node_skel_succ
            FP = major_tree.nodes[node]['nr_of_leaf_succ'] - TP
            FN = root_skel_succ - node_skel_succ
            prec = TP / (TP + FP)
            rec = TP / (TP + FN)
            F = 2 / ((1/prec) + (1/rec))
            res[node] = F
    print('Most fitting mapping of read {} is tax-id: {} (with F-score: {})'.format(iteration, max(res, key=res.get), max(res.values())))

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    local = os.getcwd()
    
    start_time = time.time()
    important_tree = make_tree(args.inpdump, species_only_leaves=True)
    print('Adding successor information on taxonomy tree...')

    for node in important_tree.nodes():
        if important_tree.in_degree(node) == 0:
            root_of_important_tree = node
    nr_of_leaf_succ(important_tree, root_of_important_tree)

    with open(args.inpread, 'r') as inp:
        it = 0
        print('Reading {} file'.format(args.inpread))
        for lines in inp.read().splitlines():
            it += 1
            leaves = lines.split(' ')[1:]
            leaves = [int(leave) for leave in leaves]
            important_skel = find_skeleton(input_graph=important_tree, reads=leaves)
            best_F(important_tree, important_skel, it)
        
        if it == 0:
            print('No skeleton trees found. Check naming conventions or specifiy input directory...')
        else:
            print('All taxonomy assignements completed!')
        print('Total program execution time: {}s'.format(time.time()-start_time))
