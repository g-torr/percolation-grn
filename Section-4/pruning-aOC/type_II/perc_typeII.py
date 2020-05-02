import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import random
import pickle
import argparse
import itertools
from multiprocessing import Pool
import scipy
import sys
sys.path.insert(0, "../../../lib")# add the library folder to the path I look for modules
#import configurational_model
import configurational_model_regulatory
import pruning
import component as co



def generate_degree_seq(gamma, N):
    kseq = np.ceil( np.random.pareto(gamma, N))
    cond = kseq > N
    while any(cond):
        temp_seq = np.ceil( np.random.pareto(gamma, np.count_nonzero(cond)))
        kseq[cond] = temp_seq
        cond = kseq > N
    return np.array(kseq, dtype=int)






def create_graph_and_remove(N1, N2, gamma_G, c, perc_fraction):
    aseq = generate_degree_seq(gamma_G, N1)
    bseq = np.random.poisson(1 + c, N2)
    if sum(bseq) < N1:
        raise ValueError('degree sequence non compatible for regulations')
    #generate the graph
    BG = configurational_model_regulatory.configuration_model_mirroring(aseq, bseq, nx.MultiDiGraph())
    Genes, TFs, in_deg = pruning.filter_out(BG)
    genes_rm = random.sample(Genes, int(perc_fraction * N1))
    BG_filtered, count = pruning.out_nonlinear_remove(BG, Genes, TFs, genes_rm)
    print('we are left with:', len(BG_filtered.nodes()), 'nodes in non-linear out-component')
    return len(set(BG_filtered.nodes()).intersection(Genes))/len(Genes), len(set(BG_filtered.nodes()).intersection(TFs))/len(TFs), np.mean(
        aseq), np.mean(bseq)
def save_obj(obj):
    with open('./dic-powerlaw-both.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def main():
   parser = argparse.ArgumentParser('generates the fraction of nodes in aOC after elimination of 5 per cent of nodes')

   parser.add_argument('--threads',type= int,default=None,help = 'number of parallel process to execute. Default behaviour is all cores availbles')
   args=parser.parse_args()
   threads= args.threads

   N1=300000
   N2=300000
   perc_fraction = 0.05
   print('After removal, network contains ',N1*(1-perc_fraction), 'genes and ',N2, 'TFs')
   gamma_Gs=np.append(np.linspace(1.9,2.6,10),[3,4,5])
   c = 0.4  # average in degree of genes is c+1
   pool=Pool(threads)
   data=pool.starmap(create_graph_and_remove,itertools.product([N1],[N2],gamma_Gs,[c],[perc_fraction]))
   print('g, t, <dout>,<cout>, powerlaw exponent')
   print(np.c_[data,gamma_Gs])
   pool.close()
   dic={'data':data,'N1':N1,'gamma_Gs':gamma_Gs,'c':c,'descr':'data is a matrix containing 4 columns formatted as follows:\n1st: fraction of genes in aOC,2nd fraction of TFs in aOC,3rd mean out-degree of genes,4th mean out-degree of TFs'}
   save_obj(dic)
if __name__ == '__main__':
    main()
