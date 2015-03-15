from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import numpy as np
import pdb


def simpleHeatmap(D):        
    Y = sch.linkage(squareform(D), method='average')
    fig = plt.figure(figsize=(8,8))
    dendrogram = fig.add_axes([0.01,0.1,0.3,0.8])
    Z = sch.dendrogram(Y, orientation='right')
    idx = Z['leaves'] ## 
    dendrogram.set_yticklabels(idx)
    
    D = D[idx,:][:,idx]    
    matrix = fig.add_axes([0.37,0.1,0.6,0.8])
    im = matrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
    matrix.set_xticklabels([])
    matrix.set_yticklabels([])
    plt.show()

def collapseIdenticalLabels(b):
    a = []
    prevEl = ''
    for e in b:
        if e.split("/")[0] != prevEl.split("/")[0]:
            a.append(e)
        else:
            a.append("")
        prevEl = e
    return a

def collapseIdenticalLabels2(b):
    a = []
    prevEl = ''
    for e in b:
        if e != prevEl:
            a.append(e)
        else:
            a.append("")
        prevEl = e
    return a

def heatmap(D, names, linkage1='single', linkage2='single', collapseLabels=None, filename=None):
    import matplotlib.pyplot as plt 
    import numpy as np
    plt.xkcd()
    # Compute and plot first dendrogram.
    fig = plt.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.01,0.1,0.15,0.6])
    Y = sch.linkage(D, method=linkage1)
    Z1 = sch.dendrogram(Y, orientation='right')
    idx1 = Z1['leaves']
    #print idx1
    #print np.array(names)[idx1]
    ax1.set_xticks([])
    #ax1.set_yticks([])
    labels = np.array(names)[idx1]
    if collapseLabels is not None:
        labels = collapseLabels(labels)
    ax1.set_yticklabels(labels)

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.37,0.71,0.6,0.2])
    Y = sch.linkage(D, method=linkage2)
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    D = squareform(D)

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.37,0.1,0.6,0.6])
    idx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    if filename:
        plt.savefig(filename)
    else:
        plt.show()

def readtxt(textfile, delim='\t'):
    matrix = [line.rstrip().split(delim) for line in open(textfile)]
    data = np.array([map(float, row[1:]) for row in matrix[1:]])
    headers = matrix[0][1:]
    col0 = [row[0] for row in matrix[1:]]
    assert col0==headers
    return data, headers

def similarSamples(D, names, linkage='complete', subcluster=10, subclusterCenter='Sabkha'):
    """selects the relevant part of a large heatmap, returns data and headers around a cluster-focus
       for  'Zooming in'"""
    # Compute dendrogram.
    Y = sch.linkage(D, method=linkage)
    Z1 = sch.dendrogram(Y, no_plot=True)
    idx = Z1['leaves']
    sim = list(np.array(names)[idx])
    pos = sim.index(subclusterCenter)
    D = squareform(D)
    nrange = idx[max(pos-subcluster,0): pos+subcluster]    
    return D[nrange,:][:,nrange], sim[max(pos-subcluster,0): pos+subcluster]

if __name__ == "__main__":
    from dunn import flatcluster, flatclusterDunn
    import sys
    #textfile = '/home/zain/Downloads/weighted_unifrac_crabgut.txt' #'/home/handreas/Data/Beta_Diversity_VM/Beta/binary_jaccard_otu_table.txt'
    textfile = '/data/EarthMicrobiomeProject/BetaDiversity/weighted_unifrac_allfreshwater_r2000_red3sig.txt'
    if len(sys.argv)>1:
        textfile = sys.argv[1]
    data, headers = readtxt(textfile)
    fc = flatclusterDunn(squareform(data))
    #heatmap(squareform(data), headers) ## one corner of the dist-mat

