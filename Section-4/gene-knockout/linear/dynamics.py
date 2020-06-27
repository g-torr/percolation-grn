from pylab import *
import networkx as nx
from networkx.algorithms import bipartite
import itertools
import warnings
import scipy
import numpy.random
import argparse
import sys
import os
sys.path.insert(0, "../../../lib")  # add the library folder to the path I look for modules
import pickle
from directed_random_graph import *
from multiprocessing import Pool


def create_graph(cout, dout, N1, N2):
    _,BG = create_bipartite(N1, N2, dout, cout)
    Genes = [n for n, d in BG.nodes(data=True) if d['bipartite'] == 0]
    TFs = [n for n, d in BG.nodes(data=True) if d['bipartite'] == 1]
    eta = bipartite.biadjacency_matrix(BG, Genes)  # from top_nodes to bottom nodes
    xi = bipartite.biadjacency_matrix(BG, TFs, format="csc")  # from  bottom nodes to  top_nodes
    J = (eta * scipy.sparse.diags(1 / array(list(dict(BG.in_degree(TFs)).values())), format="csr") * xi)
    return J


def save_obj(obj,dout,cout):
    name='cout:'+str(cout)+'dout:'+str(dout)+'.pkl'
    with open('dic-' + name , 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)





# --------This is the simulation------


def replics_parallel(J):
    pool = Pool()
    '''Simulation at fixed T for different replicas
    Initial condition is chosen to be the same as for cavity'''
    N1 = J.shape[0]


    data = []
    data = pool.starmap(dynamics_light_parallel, itertools.product([J], range(N1)))
    # for replica in range(N_replic):
    #        data+=[dynamics_light(J,psi_init,T)]
    pool.close()
    return data


def dynamics_light_parallel(J,node):

    N1 = J.shape[0]
    N_iterations =100
    n_start = ones(N1)
    n_start[node] = 0
    n = scipy.sparse.csr_matrix(n_start)
    t = 1
    while t < N_iterations:
        a = (n * J).toarray()
        n_new = where(a > 0, 1, 0)
        n_new[0][n_start == 0] = 0  # enforce elimination constrain of Genes initially eliminated
        n_new = scipy.sparse.csr_matrix(n_new)
        if n_new.nnz == n.nnz:
            if all(n_new.indices == n.indices):  # exit if no improvements occur
                break
        n = n_new
        t += 1
    return n.toarray()[0]





def load_input(args):
    N1 = args.N1  # number genes
    N2 = args.N2  # number TFs
    dout = args.dout
    cout = args.cout
    return N1, N2, cout, dout



def main():
    parser = argparse.ArgumentParser(
        description='Simulation of magnetisation for spin glass. Many replica considered. Use the argument --create_graph if you do not want to load the matrix of couplings from the dictionary. Data will be saved in the data folder. Read the Readme.md for more info ')
    parser.add_argument("-N1", help="Number of genes", type=int, const=1000, default=1000, nargs='?')
    parser.add_argument("-N2", help="Number of TF", type=int, const=1000, default=1000, nargs='?')
    parser.add_argument('--dout', type=int, default=2, help=" average out degree of genes. Default set to 2")
    parser.add_argument('--cout', type=int, default=2, help=" average out degree of TF. Default set to 2")
    args = parser.parse_args()
    N1, N2, cout, dout = load_input(args)
    J = create_graph(cout, dout, N1, N2)
    simulation = array(replics_parallel(J))
    magn_simulation = [replic for replic in simulation]
    dic = {"dout": dout, "cout": cout, "N1": N1, "N2": N2, "J": J,  "magn_simulation": magn_simulation}
    '''while True:
        flag = input("save the results? [y/n]")
        if (flag == "y" or flag == "yes"):
            save_obj(dic, "states:cin" + str(cin)+"din:" + str(din))
            break
        elif (flag == "n" or flag == "not"):
            print("I do not save the result")
            break
        else:
            "I do not understand your answer, type yes or not"
    '''
    save_obj(dic, dout,cout)


if __name__ == '__main__':
    main()
