from pylab import *
import networkx as nx
from networkx.algorithms import bipartite
import itertools
import warnings
import numpy.random
#from scipy.special import factorial
import itertools
import scipy
import sys
import pickle
import argparse
from multiprocessing import Pool
#from scipy.stats import poisson
sys.path.insert(0, "../../../../lib")# add the library folder to the path I look for modules
sys.path.insert(0, "../lib")# add the library folder to the path I look for modules
from directed_random_graph import *# this is the function I wrote
import non_linear.percolation_Genes as nl


''' 
#This run the dynamics initial condition, but it does not remove the node, which could be activated in later stages 
def dynamics_light(N1,N2,d,c,N_iter,m_start): 
    #It does not store the entire trajectory, but only the manetisation in the final state 
    def time_step(n): 
        m=where((n*eta!=in_deg).toarray(),0,1) 
        n=scipy.sparse.csr_matrix(m)*xi 
        n=where(n.toarray()>0,1,0) 
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
    while t< N_iter: 
        n,m=time_step(n) 
        t+=1 
    return mean(n.toarray()),mean(m) 
'''

def heatmap(magn,ds,cs,alpha,Genes_initially_on,histeresys=False):
    figure()
    imshow(magn[:,:,0],extent=(cs[0],cs[-1],ds[0],ds[-1]),origin="lower", aspect="auto",cmap='RdYlBu')
    colorbar()
    xlabel("c")
    ylabel("d")
    ylim(ds[0],ds[-1])
    xlim(cs[0],cs[-1])
    if histeresys:
        plot(cs,alpha *log(1 + cs*alpha))# for d larger then this,the histeresys begins
        plot(cs,(-1 +exp(cs*alpha) )*alpha) # for d larger then this,the histeresys terminate
        title("Gene  magnetisation,genes initially on - off ")
    elif Genes_initially_on:
        plot(cs,(-1 +exp(cs*alpha) )*alpha) # for d larger then this,the state with n=1 becomes instable
        title("Gene  magnetisation,TF initially on ")
    else:
        plot(cs,alpha *log(1 + cs*alpha))# for d larger then this,the state with  n=0 genes becomes stable, and n=1 remains stable
        title("Gene  magnetisation,genes initially off ")
def save_obj(obj, name ):
    with open('./dic-'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
def load_obj(name ):
    with open('./dic-' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def main():
    parser = argparse.ArgumentParser(description='Save the colormap of the magnetisation for different network realisation. Two sets of simulations are performed for the given set of parameter. The first one starts with TFS initially on means that 99% of TF are on. Second with TFS initially off which means that only 5 of TF are on.  \n Usage: python percolation.py --ds 0.01 1 10--cs 0.01 1 10 -N1 10000 -N2 10000 ')
    parser.add_argument('--ds', type=float,nargs='*', default=[0.01,1,10],help="[cmin,cmax, Npoints]. gives the list of ds to be explored for percolation. \n Defaul value set to 0 1 10" )
    parser.add_argument('--cs', type=float,nargs='*', default=[0.01,1,10],help="[cmin,cmax, Npoints]. gives the list of cs to be explored for percolation.  \n Defaul value set to 0.01 1 10" )
    parser.add_argument("-N1", help="Number of genes",type=int,const=10000, default=10000,nargs='?')
    parser.add_argument("-N2", help="Number of TF",type=int,const=10000, default=10000,nargs='?')
    parser.add_argument("--threads", type=int, default=-1, help="Number of cores to be used")
    parser.add_argument("--perc_fraction",default=0.95,type=float, help="percolation fraction p, 1 is unperturbed")
    args = parser.parse_args()
    cores = args.threads
    cs= linspace((args.cs)[0],(args.cs)[1],round((args.cs)[2]))
    ds= linspace((args.ds)[0],(args.ds)[1],round((args.ds)[2]))
    perc_fraction=args.perc_fraction
    N1=args.N1 # number genes
    N2=args.N2 # number TFs
    data_on= nl.percolation_parallel(ds,cs,N1,N2,cores,Genes_initially_on=True,perc_fraction=perc_fraction)
    data_off= nl.percolation_parallel(ds,cs,N1,N2,cores,Genes_initially_on=False,perc_fraction=perc_fraction)
    #heatmap(data_on,ds,cs,N2/N1,Genes_initially_on=True)
    #heatmap(data_off,ds,cs,N2/N1,Genes_initially_on=False)
    #heatmap(data_on -data_off,ds,cs,N2/N1,Genes_initially_on=None,histeresys=True)
    #show()
    dic={"ds":ds,"cs":cs,"N1":N1,"N2":N2,"data_on":data_on,"data_off":data_off,'perc_fraction':perc_fraction }
    '''while True: 
        flag=input("save the results? [y/n]") 
        if (flag=="y" or flag=="yes"): 
            save_obj(dic,"magnetisation") 
            break 
        elif (flag=="n" or flag=="not"): 
            print("I do not save the result") 
            break 
        else: 
            "I do not understand your answer, type yes or not" 
    '''
    save_obj(dic,"magnetisation")
if __name__ == '__main__':
    main()
