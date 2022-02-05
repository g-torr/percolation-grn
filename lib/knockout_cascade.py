import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
import itertools
import warnings
# from scipy.special import factorial
import itertools
import scipy
import sys
from multiprocessing import Pool

sys.path.insert(0, "lib")  # add the library folder to the path I look for modules
from directed_random_graph import *  # this is the function I wrote

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
