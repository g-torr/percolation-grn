import numpy as np
import argparse 
import sys
import pickle
sys.path.insert(0, "../lib")  # add the library folder to the path I look for modules
import knockout_cascade
def save_obj(obj, name ):
    with open('dic-'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def main():
    parser = argparse.ArgumentParser(
        description='Single node removal is performed starting by each of the N1 nodes.  Simulation is repeated once for each network. Code gives the fraction of nodes in out-component\n Usage: python percolation.py --d 1 --c 1  -N1 10000 -N2 10000 --Nrep 1000')
    parser.add_argument("--ds",type = float, nargs = '*', default = [0.05, 4, 40],help = "[dmin,dmax,Npoints]. Out degree of gene to be explored for percolation. \n Defaul value set to  [0.05, 1.1, 2] ")
    parser.add_argument('--cs',type = float, nargs = '*', default = [0.05, 2, 30],help = "[cmin,cmax,cpoints]. Out degree of TFs to be explored for percolation. \n Defaul value set to  [0.05, 1.1, 2] ")
   
    parser.add_argument("-N1", help="Number of genes", type=int, const=10000, default=10000, nargs='?')
    parser.add_argument("-N2", help="Number of TF", type=int, const=10000, default=10000, nargs='?')
    # parser.add_argument("--threads", type=int, default=-1, help="Number of cores to be used")
    args = parser.parse_args()
    cs = np.linspace((args.cs)[0], (args.cs)[1], int((args.cs)[2]))
    ds = np.linspace((args.ds)[0], (args.ds)[1], int((args.ds)[2]))
    N1 = args.N1  # number genes
    N2 = args.N2  # number TFs
    # cores = args.threads
    lambda_heat =[]
    for d in ds:
        for c in cs:
            data=np.array(knockout_cascade.create_graph_and_remove(N1,N2,d,c))
            lambda_heat+=[np.count_nonzero(data[:,0]>0.5)/len(data)]
    
    dic = {'lambda':np.array(lambda_heat).reshape((len(ds),len(cs))),'ds':ds,'cs':cs,'N1':N1,'N2':N2}
    save_obj(dic,'lambda_heatmap')
if __name__ == '__main__':
    main()