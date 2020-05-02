from pylab import *
sys.path.insert(0, "../lib/")# add the library folder to the path I look for modules
from directed_random_graph import *# this is the function I wrote
from multiprocessing import Pool
import itertools
import pickle
import component
def save_obj(obj):
    with open('dic-SCC.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
def SCC_nonlinear(BG,Genes,TFs,process=0):
    '''Return the subgraph in the SCC, and then removes TFs whose in degree is not the same as initially '''
    BG_filtered=BG.copy()
    BG_tf=dict(BG.in_degree(TFs))# dictionary of in degree for TF in the bipartite graph
    count=0
    while True:
        Gc= component.maximum_strongly_connected_component_subgraph(BG_filtered)
        #print(len(Gc),len(BG_filtered))
        #print(BG_filtered.nodes())
        if len(Gc)==len(BG_filtered):
            print("number of deletion iterations:",count)
            return BG_filtered,count
        count+=1
        
        Gc_G=set(Gc.nodes()).intersection(Genes)# set of genes in the giant cluster
        Gc_tf=dict(Gc.in_degree(TFs))# dictionary of in degree for TF in the SCC
        filtered={v for v,k in Gc_tf.items() if BG_tf[v]==k} # filtered dictionary of SCC if TF has same indegree as initially
        #BG_filtered.remove_nodes_from([n for n in BG_filtered if n not in set(filtered).union(Gc_G)])#remove TFs with wrond in-degree and Genes not in SCC
        BG_filtered=BG_filtered.subgraph(set(filtered).union(Gc_G)).copy()
        #print(set(filtered),Gc_G)
def main():
    N1=50000
    N2=50000
    c=0.5
    data=[]
    ds=linspace(0.01,3,40)
    for d in ds:
        _,BG=create_bipartite(N1,N2,d,c)
        TFs=arange(N1,N1+N2)
        Genes=arange(0,N1)
        SG,count=SCC_nonlinear(BG,Genes,TFs)
        data+=[[len(set(SG.nodes()).intersection(Genes)),len(set(SG.nodes()).intersection(TFs))]]
    data=array(data)
    dic = {'c':c,'N1':N1,'N2':N2,'ds':ds,'data':data}
    save_obj(dic)
if __name__ == '__main__':
    main()