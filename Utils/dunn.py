"""
Tools for assessing clustering quality.
See Wikipedia: Dunn index
"""

from scipy.spatial.distance import squareform, pdist
import scipy.cluster.hierarchy as sch
import numpy as np
#from heatmap import simpleHeatmap

def clusterDistance(cluster1, cluster2, DM):
    """implemented as average distance"""
    return DM[cluster1][:,cluster2].mean()
    
def dunn(clustering, DM):
    clusters = len(clustering)
    maxIntra = max([clusterDistance(c, c, DM) for c in clustering])
    minInter = np.inf
    for cl1 in range(1, clusters):
        for cl2 in range(cl1):
            minInter = min(clusterDistance(clustering[cl1], clustering[cl2], DM), minInter)
    return minInter/maxIntra

def compactness(cluster, clustering, DM):
    homogeneity = clusterDistance(cluster, cluster, DM)
    separation = np.mean([clusterDistance(cluster, otherCluster, DM) for otherCluster in clustering])
    return homogeneity, separation

def flatcluster(DM, threshold):
    Y = sch.linkage(DM, method="average")
    n = len(Y)+1
    clusterdict = dict([(el, [el]) for el in range(n)])
    for cidx, cluster in enumerate(Y):        
        if cluster[2] < threshold:
            cl1, cl2 = map(int, cluster[:2])
            clusterdict[cidx+n] = clusterdict[cl1] + clusterdict[cl2]
            del clusterdict[cl1]
            del clusterdict[cl2]
    return clusterdict.values()

def flatclusterDunn(D, save=None):
    DM = squareform(D)
    Y = sch.linkage(DM, method="average")
    if save: dump.dump(Y, "%s_linkage.npy" % save)
    n = len(Y)+1
    clusterdict = dict([(el, [el]) for el in range(n)])
    bestClustering, bestDunn = None, 0
    for cidx, cluster in enumerate(Y):        
        cl1, cl2 = map(int, cluster[:2])
        clusterdict[cidx+n] = clusterdict[cl1] + clusterdict[cl2]
        del clusterdict[cl1]
        del clusterdict[cl2]
        clustering = clusterdict.values()
        if not len(clustering) in [1,n]:
            currentDunn = dunn(clustering, DM)
            if save:
                dump.dump((clustering, currentDunn), "%s_%05d_%04d.pcl" % (save, cidx, len(clustering)))
            #print "Dunn(%s):"%clustering, currentDunn
            if currentDunn > bestDunn:
                bestDunn = currentDunn
                bestClustering = clustering ## maybe copy
    return bestClustering, bestDunn

def subclusterRows(Y, idx):
    """ recursively retrieves relevant rows for a specified subcluster
    Input: Y   -  a clustering matrix 4 x (n-1)
           idx -  a cluster idx according to the cluster naming induced by sch.linkage,
           eg. if we have n=5, 0-4 are the singletons, 5 refers to the cluster created in row 0 in Y,
                                                       6 refers to the cluster created in row 1 in Y,
                                                       ...
                                                       10 (2n-2)  ---''----------''-----in row 4 (n-1) in Y.
    """
    if idx > len(Y):
        rowIdx = idx-(len(Y)+1)
        row = Y[rowIdx]
        cl1, cl2 = map(int, row[:2])
        r1, m1 = subclusterRows(Y, cl1)
        r2, m2 = subclusterRows(Y, cl2)
        return r1 + r2 + [rowIdx], m1 + m2
    else:
        return [], [idx]

def subcluster(Y, idx):
    rows, members = subclusterRows(Y, idx)
    rows.sort()
    Y1 = Y[rows]
    mapping = dict(zip(sorted(members), range(len(members))))
    for i, row in enumerate(rows):
        oldClusterIdx = row + len(Y)+1
        newClusterIdx =   i + len(rows) + 1
        mapping[oldClusterIdx] = newClusterIdx
    for row in Y1:
        row[0] = mapping[row[0]]
        row[1] = mapping[row[1]]
    revMapping = dict([(v,k) for k,v in mapping.items()])
    return Y1, mapping, revMapping
        
if __name__ == "__main__":
    import matplotlib.pyplot as plt 
    #D = [1,5,6,0.9,4,5,0.1,1,4.1,5.1]
    names = ['Null', 'Eins', 'Zwei', 'Drei', 'Vier', 'Fuenf']
    D =  [1.1, 5, 6, 0.9, 3.9, 4.9, 0.2, 1, 4.1, 5.1]
    D = pdist([[0],[1.1],[5],[6],[0.9],[1]])
    DM = squareform(D)
    Y = sch.linkage(DM, method="complete")
    Y1, mapping, revMapping = subcluster(Y, 9)
    print Y
    print Y1
    Z = sch.dendrogram(Y1, leaf_label_func=lambda x: "%s-%s-%s"%(x,revMapping[x],names[revMapping[x]]))
    print Z['leaves']
    #link_color_func=lambda k: colors[k]
    plt.show()
    #fig = plt.figure(figsize=(8,8))
    #Z = sch.dendrogram(Y, truncate_mode="lastp", p=5, show_contracted=True)
    #plt.show()
    #print flatcluster(D, 1.3)    
    #print flatclusterDunn(D)
    
    #clustering = [[0,1,4], [2,3]]
    #print dunn(clustering, D)
    #clustering = [[0,1], [4], [2,3]]
    #print dunn(clustering, D)
    #clustering = [[0],[1], [4], [2],[3]]
    #print dunn(clustering, D)

"""
3 [  1.           2.           0.4472136    2.        ]
4 [  0.           3.           2.30867928   3.        ]

5 [  1.           4.           0.4472136    2.        ]
7 [  0.           5.           2.30867928   3.        ]
    
5 [  1.           4.           0.4472136    2.        ]
6 [  2.           3.           2.23606798   2.        ]
7 [  0.           5.           2.30867928   3.        ]
8 [  6.           7.          10.95810203   5.        ]]    

    The dendrogram illustrates how each cluster is
    composed by drawing a U-shaped link between a non-singleton
    cluster and its children. The height of the top of the U-link is
    the distance between its children clusters. It is also the
    cophenetic distance between original observations in the two
    children clusters. It is expected that the distances in Z[:,2] be
    monotonic, otherwise crossings appear in the dendrogram.
    
    Parameters
    ----------
    Z : ndarray
        The linkage matrix encoding the hierarchical clustering to
        render as a dendrogram. See the ``linkage`` function for more
        information on the format of ``Z``.
    p : int, optional
        The ``p`` parameter for ``truncate_mode``.
    truncate_mode : str, optional
        The dendrogram can be hard to read when the original
        observation matrix from which the linkage is derived is
        large. Truncation is used to condense the dendrogram. There
        are several modes:
    
        * None/'none': no truncation is performed (Default)
        * 'lastp': the last ``p`` non-singleton formed in the linkage
          are the only non-leaf nodes in the linkage; they correspond
          to to rows ``Z[n-p-2:end]`` in ``Z``. All other
          non-singleton clusters are contracted into leaf nodes.
        * 'mlab': This corresponds to MATLAB(TM) behavior. (not
          implemented yet)
        * 'level'/'mtica': no more than ``p`` levels of the
          dendrogram tree are displayed. This corresponds to
          Mathematica(TM) behavior.
    
    color_threshold : double, optional
        For brevity, let :math:`t` be the ``color_threshold``.
        Colors all the descendent links below a cluster node
        :math:`k` the same color if :math:`k` is the first node below
        the cut threshold :math:`t`. All links connecting nodes with
        distances greater than or equal to the threshold are colored
        blue. If :math:`t` is less than or equal to zero, all nodes
        are colored blue. If ``color_threshold`` is ``None`` or
        'default', corresponding with MATLAB(TM) behavior, the
        threshold is set to ``0.7*max(Z[:,2])``.
"""
