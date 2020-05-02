import numpy as np
import scipy
import sys
import itertools
from multiprocessing import Pool
import pickle
import random
import networkx as nx
sys.path.insert(0, "../../../lib")# add the library folder to the path I look for modules
import component as co
import  directed_random_graph # this is the module  I wrote in /lib

def out_nonlinear_remove(BG, Genes, TFs, genes_rm):
    '''Return the subgraph in the SCC, and then removes TFs whose in degree is not the same as initially '''
    BG_filtered = BG.copy()
    BG_tf = dict(BG.in_degree(TFs))  # dictionary of in degree for TF in the bipartite graph
    BG_filtered.remove_nodes_from(genes_rm)
    count = 0
    while True:
        Gc_nodes = co.out_component(BG_filtered)
        Gc = BG_filtered.subgraph(Gc_nodes).copy()
        if len(Gc_nodes) == len(BG_filtered):
            print("number of deletion iterations:", count)
            return BG_filtered, count
        count += 1

        Gc_G = set(Gc.nodes()).intersection(Genes)  # set of genes in the giant cluster
        Gc_tf = dict(Gc.in_degree(TFs))  # dictionary of in degree for TF in the SCC
        filtered = {v for v, k in Gc_tf.items() if
                    BG_tf[v] == k}  # filtered dictionary of SCC if TF has same indegree as initially
        if len(filtered) == 0:
            return nx.empty_graph(), count
        BG_filtered = BG_filtered.subgraph(set(filtered).union(Gc_G)).copy()



def filter_out(BG):
    '''It extract the out subgraph'''
    nodes_out = co.out_component(BG)
    Gc = BG.subgraph(nodes_out).copy()
    Genes = [n for n, d in Gc.nodes(data=True) if d['bipartite'] == 0]
    N1 = len(Genes)
    TFs = [n for n, d in Gc.nodes(data=True) if d['bipartite'] == 1]
    N2 = len(TFs)
    in_deg = scipy.sparse.csr_matrix(list(dict(Gc.in_degree(TFs)).values()))
    return Genes, TFs, in_deg


def create_graph_and_remove(N1,N2,d,c,perc_fraction):
    _,BG = directed_random_graph.create_bipartite(N1,N2,d,c)
    Genes,TFs,in_deg=filter_out(BG)
    print('Starting from ',len(Genes)+len(TFs), 'nodes in out component')
    genes_rm = random.sample(Genes,int(perc_fraction*N1))
    print(' we remove ',len(genes_rm),' genes initially')
    BG_filtered,count = out_nonlinear_remove(BG,Genes,TFs,genes_rm)
    print('we are left with:',len(BG_filtered.nodes()), 'nodes in non-linear out-component')
    return len(set(BG_filtered.nodes()).intersection(Genes)),len(set(BG_filtered.nodes()).intersection(TFs))
def save_obj(obj, name ):
    with open('dic-'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
def main():
    N1=10000
    N2=10000
    c = 2
    perc_fraction=0.05
    ds=np.append(np.linspace(0.01,1.5,10),np.linspace(1.51,2.5,20))
    pool = Pool()
    data=pool.starmap(create_graph_and_remove,itertools.product ([N1],[N2],ds,[c],[perc_fraction]))
    dic={'data':data,'ds':ds,'c':c,'perc_fraction':perc_fraction,'N1':N1}
    pool.close()
    save_obj(dic,'perc-typeI')


if __name__ == '__main__':
    main()
