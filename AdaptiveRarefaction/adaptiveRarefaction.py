"""
Adaptive rarefaction for phylogenetic Beta-diversity
Author: Andreas Henschel
License: GPL v3.0\n\n

Given a number of microbial community samples (OTU vs sample Biom table) of possibly strongly varying sampling depth, all individual sample pairs are rarefied to a sample size necessary only for that pair, then weighted UniFrac is called on that pair only.
The result is an entry to the total beta-diversity matrix, which at the end is outputted in csv (tab-delimited) format, adhering to Qiime's way of storing distance matrices, such that it can be further processed with Qiime scripts (e.g. principal_coordinates.py)
Built-in jack-knifing: if number of repetitions are provided, random subsampling is done repeatedly and a series of distance matrices is dumped in the provided output directory.
I use the entire greengenes97 tree, unifrac trims the tree to the relevant part. Note: OTUs must match the tree labels.

Alternatively, instead of a biom table, community sample data can be provided in other formats:
1. pickled/numpy files for each sample, the first containing a list of OTU names, the second (in same order a numpy array containing abundances).
The
2. As initally employed through the database, we also provide programmatic access to the SQL tables storing abundance information,
If a SQL server is hosted locally or elsewhere, you would need to change the access details.
Not parallel. See also parallel version, because this turned out to be very slow for large amounts of samples.

Wraps around FastUniFrac (as provided in cogent)
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2797552/

Dependencies:
Requires PyCogent, Biom v2.1 (http://biom-format.org/), NumPy

TODO:
"""

import numpy as np
from cogent.parse.tree import DndParser
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent.maths.unifrac.fast_tree import UniFracTreeNode
from collections import Counter
import glob, os
import dump
from biom.parse import parse_biom_table

class Sample:
    def __init__(self, sampleID, datadir):
        ## from predefined dumps
        self.sampleID = sampleID
        self.otus = dump.load("%s/%s_otus.pcl" %(datadir, sampleID))
        self.seqCounts = np.load("%s/%s_sc.npy" %(datadir, sampleID))
        self.size = self.seqCounts.sum()
        self.seqP = self.seqCounts/float(self.seqCounts.sum())
    def subsample(self, size):
        return Counter(np.random.choice(self.otus, size, replace=True, p=self.seqP))

class BiomSample(Sample):
    def __init__(self, sampleID, s, observationIDs):
        self.sampleID = sampleID
        self.otus = observationIDs[s.indices]
        self.seqCounts = s.data
        self.size = self.seqCounts.sum()
        self.seqP = self.seqCounts/float(self.seqCounts.sum())        
        
def makeOTUdict(otu, subsample1, subsample2, sample1, sample2):
    otuDict = {}
    if otu in subsample1.keys(): otuDict[sample1.sampleID] = subsample1[otu]
    if otu in subsample2.keys(): otuDict[sample2.sampleID] = subsample2[otu]
    return otuDict

def unifrac2(sample1, sample2, tree, repetitions=1, subsampleSize='auto'):
    distances = []
    if subsampleSize == 'auto':
        subsampleSize = int(min(sample1.size, sample2.size)*.8) # 80% of the smaller sample
    for i in range(repetitions):
        subsample1  = sample1.subsample(subsampleSize) 
        subsample2  = sample2.subsample(subsampleSize) 
        allOtus = set(subsample1.keys()).union(subsample2.keys())
        envs = dict([(otu, makeOTUdict(otu, subsample1, subsample2, sample1, sample2)) for otu in allOtus])
        #pdb.set_trace()
        res = fast_unifrac(tree, envs, weighted=True)
        try:
            distances.append(res['distance_matrix'][0][0,1])
        except:
            pdb.set_trace()
    print subsampleSize, np.array(distances).mean()
    return np.array(distances)


if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-b", "--biom", help= "OTU biom table containing abundances of OTUs for a number of samples")
    parser.add_argument("-o", "--outdir", help= "directory where the output (distance matrices) should go")
    parser.add_argument("-t", "--treefile", help= "Phylogenetic tree in Newick format")
    parser.add_argument("-r", "--repetitions", help= "Number of repetitions (multiple adaptive rarefaction)", type=int, default=1)
    args = parser.parse_args()
    if not args.biom: print "No biom table provided"
    if not args.outdir: print "No output directory provided"
    if not args.treefile: print "No tree file provided"
    if not (args.biom and args.outdir and args.treefile): sys.exit(1)
    
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir) ## if this fails, at least an early death (could do with some error msg)
    #outdir = "/data/Projects/EcoPhyl/Results/AdaptiveRarefaction"
    #treefile = "/data/GreenGenesSQL/97_otus.tree"    
    tree = DndParser(open(args.treefile), UniFracTreeNode)
    #biom = "/home/zain/Projects/KnowYourEnv/Data/BiomTables/somesoils.biom"
    try:
        otutable = parse_biom_table(open(args.biom, "U"))
    except:
        raise Exception("Can't parse Biom table %s" % args.biom)    
    sampleIDs = otutable._sample_index.keys()
    ## Select all sample IDs of samples to be compared to each other
    D = [np.zeros((len(sampleIDs), len(sampleIDs))) for i in range(args.repetitions)]
    samples = [BiomSample(sampleID, otutable.data(sampleID, axis='sample', dense=False), otutable._observation_ids) for sampleID in sampleIDs]

    for i in range(1, len(sampleIDs)):
        for j in range(i):
            try:
                results = unifrac2(samples[i], samples[j], tree, repetitions=len(D))
                for k in range(len(D)):
                    D[k][i,j] = results[k]
            except ValueError as e:
                print e
                for k in range(len(D)):
                    D[k][i,j] = np.nan

    delim = "\t"
    header = delim + delim.join(sampleIDs)

    ## writing output matrixes to csv files
        
    for e, m in enumerate(D):
        m += m.T 
        w = open("%s/%s_DM%03d.csv" % (args.outdir, os.path.splitext(os.path.basename(args.biom))[0], e), "w")    
        print >> w, header    
        for ri, row in enumerate(D[e]):
            print >> w, sampleIDs[ri] + "\t" + "\t".join(map(str, row))
        w.close()
        #np.savetxt(, np.hstack((rows, m)), delimiter=delim, header=header, comments='')
#np.save("%s/%s_DM" % (smallbiom[:-5], outdir), D)
