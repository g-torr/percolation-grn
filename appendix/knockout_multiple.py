from pylab import *
import networkx as nx
from networkx.algorithms import bipartite
import itertools
import warnings
import numpy.random
# from scipy.special import factorial
import itertools
import scipy
import sys
import pickle
import argparse
from multiprocessing import Pool

sys.path.insert(0, "../lib")  # add the library folder to the path I look for modules
from directed_random_graph import *  # this is the function I wrote
from pathlib import Path
folder_name='knockout-cascade_multiple'
Path(folder_name).mkdir(parents=True, exist_ok=True)# create folder if it doesn't exist

def save_obj(obj, name):
    with open(folder_name+'/dic-' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open('dic-' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def remove_dynamics(eta, xi, in_deg, node, process=0):
    N1 = eta.shape[0]

    def time_step_on(n, n_start):
        m = np.where((n * eta != in_deg).toarray(), 0, 1)
        n = scipy.sparse.csr_matrix(m) * xi
        n = np.where(n.toarray() > 0, 1, 0)
        n[0][n_start == 0] = 0  # enforce elimination constrain of Genes initially eliminated
        n = scipy.sparse.csr_matrix(n)
        return n, m  # it is in the array format

    t = 0
    n_start = np.ones(N1)
    n_start[node] = 0
    n = scipy.sparse.csr_matrix(n_start)
    N_iter = 1000
    while t < N_iter:
        n_new, m = time_step_on(n, n_start)
        if n_new.nnz == n.nnz:
            if all(n_new.indices == n.indices):  # exit if no improvements occur
                break
        n = n_new
        t += 1
    return n.toarray()[0], m[0]


def out_statistics(eta, xi, in_deg, node):
    n, m = remove_dynamics(eta, xi, in_deg, node)
    return [np.count_nonzero(n > 0) / len(n), np.count_nonzero(m > 0) / len(m)]

'''
Not using it
def filter_SCC(BG, Genes):
    #It extract the SCC subgraph
    Gc = maximum_strongly_connected_component_subgraph(BG)
    Gc_G = set(Gc.nodes()).intersection(Genes)  # set of genes in the giant cluster
    Genes = [n for n, d in Gc.nodes(data=True) if d['bipartite'] == 0]
    N1 = len(Genes)
    TFs = [n for n, d in Gc.nodes(data=True) if d['bipartite'] == 1]
    N2 = len(TFs)
    eta = bipartite.biadjacency_matrix(Gc, Genes)  # from genes  to TF
    xi = bipartite.biadjacency_matrix(Gc, TFs)  # from  TF to  genes
    in_deg = scipy.sparse.csr_matrix(list(dict(Gc.in_degree(TFs)).values()))
    return Genes, TFs, eta, xi, in_deg
'''

def create_graph_and_remove(N1, N2, d, c):
    _, BG = create_bipartite(N1, N2, d, c)
    Genes = [n for n, d in BG.nodes(data=True) if d['bipartite'] == 0]
    TFs = [n for n, d in BG.nodes(data=True) if d['bipartite'] == 1]
    pool = Pool()
    eta = bipartite.biadjacency_matrix(BG, Genes)  # from genes  to TF
    xi = bipartite.biadjacency_matrix(BG, TFs)  # from  TF to  genes
    in_deg = scipy.sparse.csr_matrix(list(dict(BG.in_degree(TFs)).values()))
    data = pool.starmap(out_statistics, itertools.product([eta], [xi], [in_deg], range(len(Genes))))
    pool.close()
    return data


def main():
    parser = argparse.ArgumentParser(
        description='Single node removal is performed starting by each of the N1 nodes.  Simulation is repeated Nrep times. Code gives the fraction of nodes in out-component\n Usage: python percolation.py --d 1 --c 1  -N1 10000 -N2 10000 --Nrep 1000')
    parser.add_argument('--d', type=float, default=1,
                        help="Out degree of gene to be explored for percolation. \n Defaul value set to  1 ")
    parser.add_argument('--c', type=float, default=1, help="Out degree of TF c, float \n Defaul value set to 1")
    parser.add_argument("-N1", help="Number of genes", type=int, const=10000, default=10000, nargs='?')
    parser.add_argument("-N2", help="Number of TF", type=int, const=10000, default=10000, nargs='?')
    parser.add_argument("-Nrep", help="Number of repetition", type=int, const=100, default=100, nargs='?')
    # parser.add_argument("--threads", type=int, default=-1, help="Number of cores to be used")
    args = parser.parse_args()
    c = args.c
    d = args.d
    N1 = args.N1  # number genes
    N2 = args.N2  # number TFs
    Nrep = args.Nrep
    # cores = args.threads

    g = []
    for i in range(Nrep):
        data = array(create_graph_and_remove(N1, N2, d, c))
        x = np.array(data)[:, 0]
        g += [np.round((1 - x[x > 0.5]) * N1)]  # this is the number of nodes in the out-component
    description='g is the number of nodes that do not belong to out-component as a result of each elimination experiment. Rows of g are obtained for different networks,each row has dimension N1 and contains the result for each of the node removal experiment'
    dic = {"d": d, "c": c, "N1": N1, "N2": N2, "g": g, 'Nrep': Nrep,'description': description}
    save_obj(dic,'c:'+str(c)+' d:'+str(d))

if __name__ == '__main__':
    main()