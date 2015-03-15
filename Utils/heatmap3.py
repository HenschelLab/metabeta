"""
More feature rich version than originial heatmap.py:

can take a list of colorvectors which will display additional features for the clustered elements:
used by betadiversity plots, where query samples, their Ecosystem and envoclass will be colored
"""

from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import numpy as np
import pdb

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

def annotatedDendrogram(Y, names, colorVectors):
    """clustering already done"""
    import pylab
    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.01,0.02,0.5,0.96])
    Z = sch.dendrogram(Y, orientation='right')
    order = Z['leaves']
    labels = np.array(names)[order] if names else []
    labelOffset = 0.22 if names else 0
    barwidth = (1 - (0.57 + labelOffset))/len(colorVectors)
    ax1.set_xticks([])    
    ax1.set_yticklabels(labels)
    for idx, colorVector in enumerate(colorVectors):
        colorVector = np.array(colorVector)[order]
        
        axc = fig.add_axes([0.53 + labelOffset+idx*barwidth,0.02,barwidth*.9,0.96])
        y_pos = np.arange(len(colorVector))
        values = np.ones(len(colorVector))
        widths = np.ones(len(colorVector))
        axc.set_xticks([])
        axc.set_yticks([])
        axc.barh(y_pos, values, widths,linewidth=0, color=colorVector)
        
def annotatedDendrogram2(Y, names, colorVectors): ## colorVectors are allowed to be multidimensional
    """clustering already done"""
    import pylab
    ecosystemColors = {'Plant': 'DarkGreen', 'Geothermal': 'SaddleBrown', 'Soil': 'Gold', 'Biofilm': 'SlateGray', 'Animal/Human': 'DarkViolet', 'Freshwater': 'b', 'Marine': 'Cyan', 'Anthropogenic': 'DarkOrange', 'Air': 'AliceBlue', 'Hypersaline':'r'}
    ecosystemColors1 = ecosystemColors.values()
    ecosystemsIndex = dict([(eco,idx) for idx,eco in enumerate(ecosystemColors.keys())])
    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.01,0.02,0.5,0.96])
    Z = sch.dendrogram(Y, orientation='right')
    order = Z['leaves']
    labels = np.array(names)[order] if names else []
    labelOffset = 0.22 if names else 0
    barwidth = (1 - (0.57 + labelOffset))/len(colorVectors)
    ax1.set_xticks([])    
    ax1.set_yticklabels(labels)
    for idx, colorVector in enumerate(colorVectors):
        colorVector = np.array(colorVector)[order]
        
        axc = fig.add_axes([0.53 + labelOffset+idx*barwidth,0.02,barwidth*.9,0.96])
        bottom=np.zeros(len(colorVector))
        axc.set_xticks([])
        axc.set_yticks([])
        ## TO BE TESTED!
        for colorIdx, colorVectorColumn in enumerate(colorVector):
            y_pos = np.arange(len(colorVector))
            values = np.ones(len(colorVector))
            widths = np.ones(len(colorVector))
            axc.barh(y_pos, colorVectorColumn, widths, linewidth=0, color=ecosystemColors1[colorIdx])
            bottom += colorVectorColumn
            
def annotatedDendrogramInfomatrix(Y, names, infomatrix, dendrogramColors, colorVectors=[], vectorColors=[]): ## colorVectors are allowed to be multidimensional
    """clustering already done
    vectorColors and colorVectors should be of equal length
    vectorColors expected to be in some fixed order ['Blue', 'Marine', 'Gold', ...]
    colorVectors are matrices of the form [0.5, 0, 0, 0, 0.5, 0, 0, ... ] (same length as the corresponging vector color)
    creates a horizontal bar diagram for these vectors for each dendrogram element (sample)
    """
    import pylab
    import matplotlib.cm as cm
    
    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.01,0.02,0.5,0.96])
    Z = sch.dendrogram(Y, orientation='right', link_color_func=lambda k: dendrogramColors.get(k, 'k'))
    order = Z['leaves']
    labels = np.array(names)[order] if names else []
    labelOffset = 0.22 if names else 0
    matrixwidth = (1 - (0.7 + labelOffset))
    barwidth = (1 - (0.87 + labelOffset))/len(colorVectors) if colorVectors else 0 ## could be negative, if names, so fix that later!
    ax1.set_xticks([])    
    ax1.set_yticklabels(labels)
    axm = fig.add_axes([0.67 + labelOffset, 0.02, matrixwidth*.95,0.96])
    axm.imshow(infomatrix[order][::-1], interpolation='nearest', cmap=cm.get_cmap("Greys"), aspect='auto')
    axm.set_xticks([])    
    axm.set_yticklabels([])
    
def annotatedDendrogramInfomatrices(Y, names, infomatrices, colormaps, matrixwidths, dendrogramColors): 
    """clustering already done
    creates a vertical bar diagram for these vectors for each dendrogram element (sample)
    """
    import pylab
    import matplotlib.cm as cm
    
    fig = pylab.figure(figsize=(8,8))
    ax1 = fig.add_axes([0.00,0.00,0.5,1])
    Z = sch.dendrogram(Y, orientation='right', link_color_func=lambda k: dendrogramColors.get(k, 'k'))
    order = Z['leaves']
    labels = np.array(names)[order] if names else []
    labelOffset = 0.22 if names else 0
    
    ax1.set_xticks([])    
    ax1.set_yticklabels(labels)
    matrixPos = 0.525 + labelOffset
    for idx, (infomatrix, colormap, matrixwidth) in enumerate(zip(infomatrices, colormaps, matrixwidths)):
        matrixwidth = (1 - (0.53 + labelOffset))*matrixwidth
        axm = fig.add_axes([matrixPos, 0.0, matrixwidth,1.0])
        matrixPos += 0.001 + matrixwidth
        
        axm.imshow(infomatrix[order][::-1], interpolation='nearest', cmap=cm.get_cmap(colormap), aspect='auto') # "Greys"
        axm.set_xticks([])    
        axm.set_yticklabels([])
    return order

def annotatedDendrogramInfomatrices2(Y, names, infomatrices, colormaps, matrixwidths, dendrogramColors): 
    """clustering already done
    creates a horizontal bar diagram for these vectors for each dendrogram element (sample)
    """
    import pylab
    import matplotlib.cm as cm
    
    fig = pylab.figure(figsize=(8,1))
    ax1 = fig.add_axes([0.00,0.25,1,0.75])
    Z = sch.dendrogram(Y, orientation='top', link_color_func=lambda k: dendrogramColors.get(k, 'k'))
    order = Z['leaves']
    labels = np.array(names)[order] if names else []
    labelOffset = 0.22 if names else 0
    
    ax1.set_yticks([])    
    ax1.set_xticklabels(labels)
    matrixPos = 1 - (0.525 + labelOffset)
    for idx, (infomatrix, colormap, matrixwidth) in enumerate(zip(infomatrices, colormaps, matrixwidths)):
        matrixwidth = (1 - (0.53 + labelOffset))*matrixwidth
        axm = fig.add_axes([0.0, matrixPos-matrixwidth, 1.0, matrixwidth])
        matrixPos -= 0.001 + matrixwidth
        try:
            axm.imshow(infomatrix[:,order], interpolation='nearest', cmap=cm.get_cmap(colormap), aspect='auto') # "Greys"
        except:
            pdb.set_trace()
        axm.set_yticks([])    
        axm.set_xticklabels([])
    return order
    

       
def heatmap(D, names, linkage1='single', linkage2='single', collapseLabels=None, filename=None, colorVectors=None):
    import pylab
    import numpy as np
    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(8,8))
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
    #ax1.set_yticklabels(labels)

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
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    for idx, colorVector in enumerate(colorVectors):
        colorVector = np.array(colorVector)[idx1] ## sort according to dendrogram!
        axc = fig.add_axes([0.18+idx*0.06,0.1,0.05,0.6])
        y_pos = np.arange(len(colorVector))
        values = np.ones(len(colorVector))
        widths = np.ones(len(colorVector))
        axc.set_xticks([])
        axc.set_yticks([])
        axc.barh(y_pos, values, widths,linewidth=0, color=colorVector)
        
    if filename:
        pylab.savefig(filename)
    else:
        pylab.show()

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
    textfile = '/home/zain/Downloads/weighted_unifrac_crabgut.txt' #'/home/handreas/Data/Beta_Diversity_VM/Beta/binary_jaccard_otu_table.txt'
    data, headers = readtxt(textfile)
    #ndata, nheaders = similarSamples(squareform(data), headers, subclusterCenter="e44b8e811f6da5b9f082452be250aa3f")
    colorVector = [{'B':'r', 'S':'b'}[h.split(".")[-1][0]] for h in headers]
    colorVector2= [{'F':'y', 'M':'g'}[h.split(".")[-1][1]] for h in headers]
    Y = sch.linkage(squareform(data))
    annotatedDendrogram(Y, [], colorVectors=[colorVector, colorVector2])
    #heatmap(squareform(data), headers, colorVectors=[colorVector,colorVector2,colorVector]) ## one corner of the dist-mat
    #heatmap(squareform(ndata), nheaders) 
