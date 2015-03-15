from scipy.spatial.distance import squareform
import numpy as np

def similarSamples(D, names, linkage='complete', subcluster=10, subclusterCenter='Sabkha'):
    import scipy.cluster.hierarchy as sch
    import numpy as np
    # Compute dendrogram.
    Y = sch.linkage(D, method=linkage)
    Z1 = sch.dendrogram(Y, orientation='right')
    idx = Z1['leaves']
    sim = list(np.array(names)[idx])
    pos = sim.index(subclusterCenter)
    D = squareform(D)
    nrange = idx[max(pos-subcluster,0): pos+subcluster]    
    return D[nrange, nrange], sim[max(pos-subcluster,0): pos+subcluster]

def readtxt(textfile, delim='\t'):
    matrix = [line.rstrip().split(delim) for line in open(textfile)]
    data = np.array([map(float, row[1:]) for row in matrix[1:]])
    headers = matrix[0][1:]
    col0 = [row[0] for row in matrix[1:]]
    assert col0==headers
    return data, headers
    
if __name__ == "__main__":
    textfile = '/home/handreas/Data/Beta_Diversity_VM/Beta/binary_jaccard_otu_table.txt'
    data, headers = readtxt(textfile)
    print similarSamples(squareform(data[:100,:100]), headers[:100], subclusterCenter="e44b8e811f6da5b9f082452be250aa3f")
