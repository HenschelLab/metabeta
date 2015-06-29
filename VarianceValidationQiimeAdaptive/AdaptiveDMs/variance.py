import glob
import numpy as np
from scipy.spatial.distance import squareform
from scipy.stats import ttest_ind, ttest_rel
import pdb
from pylab import *

ardir =  "/data/Projects/EcoPhyl/Results/ValidationAdaptive/AdaptiveRarefaction/AdaptiveDMs"
qrdir = "/data/Projects/EcoPhyl/Results/ValidationAdaptive/Qiime/RarifiedDMs"

def getVariances(pattern):
    oldnames = None
    m3 = []
    for armatrix in glob.glob(pattern):
        m = np.genfromtxt(armatrix, skip_header=1, dtype=None)
        names =  [row[0] for row in  m]
        assert not oldnames or names == oldnames, "sample names not in the same order!!!"
        oldnames = names
        m0 = [list(row)[1:] for row in  m] 
        m3.append(m0)
    m3 = np.array(m3)
    return squareform(m3.std(axis=0))

arv = getVariances("%s/*.csv"%ardir)
qrv = getVariances("%s/*.txt"%qrdir)

print ttest_ind(arv, qrv)
print ttest_rel(arv, qrv)

figure()
boxplot([arv, qrv])
ylabel("Variance for corresponding pairwise sample distance")
xticks([1,2], ["Adaptive", "Static"])

show()
