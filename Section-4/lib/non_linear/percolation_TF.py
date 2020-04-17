from pylab import *
from networkx.algorithms import bipartite
from multiprocessing import Pool
import sys
import scipy
sys.path.insert(0, "../")# add the library folder to the path I look for modules
import numpy as np
import itertools
from directed_random_graph import *
sys.path.insert(0, "../")# add the library folder to the path I look for modules

def initialise_genes(m_start,xi):
    n=scipy.sparse.csr_matrix(m_start)*xi
    n=np.where(n.toarray()>0,1,0)
    n=scipy.sparse.csr_matrix(n)
    return  n
def dynamics_light_perc(N1,N2,d,c,N_iter,m_start,TFs_initially_on):
    '''It enforce initial nodes to remain off '''
    def time_step_on(n):
        m=np.where((n*eta!=in_deg).toarray(),0,1)
        m[0][m_start==0]=0# enforce elimination constrain of TF initially eliminated
        n=scipy.sparse.csr_matrix(m)*xi
        n=np.where(n.toarray()>0,1,0)
        n=scipy.sparse.csr_matrix(n)
        return n,m # it is in the array format
    def time_step_off(n):
        m=np.where((n*eta!=in_deg).toarray(),0,1)
        m[0][m_start==1]=1# enforce retain of TF initially activated
        n=scipy.sparse.csr_matrix(m)*xi
        n=np.where(n.toarray()>0,1,0)
        n=scipy.sparse.csr_matrix(n)
        return n,m # it is in the array format

    Genes = arange(0,N1)
    TFs = arange(N1,N1+N2)
    _,F=create_bipartite(N1,N2,d,c)
    eta=bipartite.biadjacency_matrix(F,Genes)#from genes  to TF
    xi=bipartite.biadjacency_matrix(F,TFs)#from  TF to  genes
    in_deg=scipy.sparse.csr_matrix(list(dict(F.in_degree(TFs)).values()))
    t=0
    n=initialise_genes(m_start,xi)
    m_start=m_start.toarray()[0]
    if TFs_initially_on:
        time_step=time_step_on
    else:
        time_step=time_step_off
    while t< N_iter:
        n_new,m=time_step(n)
        if n_new.nnz==n.nnz:
           if all(n_new.indices==n.indices):# exit if no improvements occur
               return np.mean(n.toarray()), np.mean(m)
        n=n_new
        t+=1
    return np.mean(n.toarray()),np.mean(m)
def percolation_parallel(ds,cs,N1,N2,cores,TFs_initially_on=True,N_iter=1000,perc_fraction=0.99):
    if cores < 1:
        pool = Pool()
    else:
        pool = Pool(processes=cores)
    if TFs_initially_on:
        m_start=scipy.sparse.random(1,N2,density=perc_fraction)
    else:
        m_start=scipy.sparse.random(1,N2,density=1-perc_fraction)
    m_start.data[:] = 1 # turn on the genes accordingly to the specified protocol
    data=pool.starmap(dynamics_light_perc,itertools.product([N1],[N2],ds,cs,[N_iter],[m_start],[TFs_initially_on]))
    data=reshape(data,(len(ds),len(cs),2))
    pool.close()
    return data


