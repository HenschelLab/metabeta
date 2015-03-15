"""
Posthoc cluster analysis and visualization.
Author: Andreas Henschel
License: GPL v3.0\n\n

After hierarchical clustering (upgma.py) of beta-diversity matrices (adaptiveRarefaction.py), we here analyse bottom-up all cluster-formations (clusterings), if they reconstitute (more or less)
ecosystems and sub-ecosystems as defined by EnvO by testing every cluster from the dendrogram for enriched EnvO terms.
This also includes possible over-arching parent concepts (a cluster might be formed of samples from feces, rumen, cecum, all of which have the parent class 'gut'.

Create dictionary (globalSampleDict): sample-id: Sample class

Input
Sample usage (Linux):

## for all samples of at least size 2000 (you can download the necessary samplenames/distance matrix/linkage files from ecophyl.info/html/DataDownload)
/usr/bin/python posthocEnvoEnrichmentTest.py -d ../Data/Clusterings/allsamples_02000__sampleNames.npy -l ../Data/Clusterings/upgma_linkage_samples_02000__sampleNames.npy -n ../Data/Clusterings/upgma_linkage_samples_02000__sampleNames.pcl -c ../Data/clusterLocal.rc -o ../Results

## for all Plant samples
/usr/bin/python posthocEnvoEnrichmentTest.py -l ../Data/Clusterings/Plant_unifrac_UPGMA.npy -d ../Data/Clusterings/Plant_unifrac.npy -n ../Data/Clusterings/Plant_unifrac.pcl -c ../Data/clusterLocal.rc -o /tmp/Posthoc

Warning: do not use large beta-diversity clusterings with the remote database server, many queries are required!

Output:
* A csv spread sheet is generated contain Clustering statistics (compactness in terms of homogeneity and separation)
* HTML web page including EnvO-colored Dendrogram including identified 'enriched clusters'
with tooltips (onMouseOver) describing samples
[Problem with large datasets: either less than one pixel per sample or very large dendrogram]


Dependencies:
  Python 2.7
  MySQL with ServerMicroBiome database needs to be installed and configured (see below).
  It will be possible (albeit slow) to run this script through the remote MySQL server at ecophyl.info. See the provided cluster.rc file.
  If a local MySQL is possible (recommended), adapt the cluster.rc file accordingly.
  Environmental Ontology (EnvO) files (envo5Terms.pcl, text descriptions of EnvO terms in Python's pickle formatted)
                                       envo5.pcl - relevant part of EnvO extracted from envo.obo, EnvO in OBO format, available on obofoundry.org)
  networkx for EnvO graph related stuff (topology_sort, general to specific ordering, traversals...)

  Also, for HTML generation, imagemagick's identify is required to be available. Probably only works on Linux.

## The script also imports a number of self made modules, make sure they are in the PYTHONPATH! (Utils directory)
dump, dunn, heatmap3, envoTools

TODO:
Fisher's exact test could be optimized:
the numbers for the contingency table could be inferred from the joined subclusters, rather then calculating this new every time
I.e., the cluster object should provide that information

"""

from matplotlib.colors import cnames, hex2color
import sys, os, glob
import numpy as np
import MySQLdb
from MySQLdb.cursors import DictCursor
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import networkx as nx
from collections import Counter, defaultdict ## requires Python 2.7
import pdb
import colorsys
from scipy.stats import fisher_exact

## self made modules, make sure they are in the PYTHONPATH
import dump
import dunn
from heatmap3 import annotatedDendrogram, annotatedDendrogramInfomatrix, annotatedDendrogramInfomatrices, annotatedDendrogramInfomatrices2
from envoTools import Ontology


areaHTML = '''<area class="hotspot" shape="rect" coords="%s,%s,%s,%s" 
onmouseover="tooltip.show('%s');"
onmouseout="tooltip.hide()" />\n'''

headerHTML = '''<html>
<head>
<script type="text/javascript" language="javascript" src="tooltip.js"></script>
<link rel="stylesheet" type="text/css" href="style.css" />
</head>
<body>
<H2>Dendrogram of Microbial Communities - Ecosystem: %s </H2>
<P>Hierarchical clustering (UPGMA) was used to produce a comprehensive mapping of microbial communities belonging to a specified ecosystem.</P>

<STRONG>Hover over the samples to get more sample information (Sample ID, study, assigned EnvO-Terms, assigned Ecosystems, sample description, and isolation source).</STRONG><BR><BR>
The color bars below the dendrogram indicate: study (random coloring, allows to visualize samples belonging to the same study), Envo coloring, Ecosystem coloring.
Note that the latter two indicate, when a sample is assigned several EnvO terms or ecosystems, respectively.
<img src="%s" alt="Dendrogram" usemap="#dbcanmap" border="0"/>
<map name="dbcanmap">'''

def generateHTML(image, sampleInfo, xsize, ysize, title="", sf=1):
    htmlfile = "%s/HTML/environmentDendrogram_%s_%s.html" % (outdir, environment, minSampleSize)
    print "Producing", htmlfile
    f = open(htmlfile, "w")
    print >> f, headerHTML%(title, image)
    sf = xsize/float(len(sampleInfo))
    for idx, info in enumerate(sampleInfo):
        print >> f, areaHTML % (1+int(idx*sf), int(0.2*ysize), 1+int((idx+1)*sf), int(0.8*ysize), info)
    print >> f, "</map></body></html>"
    f.close()


class Cluster:
    """This class describes a cluster as a list of Sample objects that 
    
    """
    def __init__(self, sampleIDs, mostly=0.0):
        self.samples = [globalSampleDict[sampleID] for sampleID in sampleIDs]
        self.enrichedIn = []
        self.mostly = mostly

    def len(self):
        return len(self.samples)
    def studies(self):
        """Returns the study Id the samples stem from"""
        return set([sample.study for sample in self.samples])
    def analyze(self, verbose=True, doFisher=True):
        """check if a cluster has a predominant portion from one ecosystem/envo-term
        and vice versa, most of this envo class is in the cluster
        makes also sure to report decently sized envo groups only.
        """
        envoCounter = Counter()
        for sample in self.samples:
            envoCounter.update(sample.envoParentIDs)

        ## candidates for potentially dominant EnvO categories
        ## I.e., dont do Fisher for all, many would be considered "significant"
        self.dominantEnvos = [(envoID, count) for (envoID, count) in envoCounter.items() if count>len(self.samples)*self.mostly]
        ## if most of an envo class is contained in a cluster: report!
        ## i.e if cluster.dominantEnvos > 90% of global distribution - report
        for envoID, envoCount in self.dominantEnvos:
            if envoCount > len(self.samples)*self.mostly: ## redundant(?)          
                ## getting data for 2x2 contingency table for Fisher's test # example numbers are for EnvO term rhizosphere
                pvalue = np.NaN                
                if doFisher:
                    notThisEnvo_cluster = len(self.samples) - envoCount              # nr of samples without specific EnvoID in current cluster #6     
                    rest = len(names) - len(self.samples)                            # nr of samples outside cluster # 1891
                    envo_outsideCluster = len(envo2samples[envoID]) - envoCount      # how many times envo term occurs outside current cluster #462 - 196 = 266
                    notThisEnvo_outsideCluster = rest - envo_outsideCluster          # 1891 - 266 = 1625
                    try:
                        odds, pvalue = fisher_exact([[envoCount, notThisEnvo_cluster], [envo_outsideCluster, notThisEnvo_outsideCluster]])
                    except ValueError:
                        pvalue = np.NaN
                precision = float(envoCount)/len(self.samples)
                recall = float(envoCount)/len(envo2samples[envoID])
                fmeasure = 2*precision*recall/(precision + recall)                
                #if (fmeasure > 0.5) and envo.G.node[envoID]['daglevel'] > 0 and envoID not in ['ENVO:00002009'] and len(self.studies()) > 1: ## includes "enriched" clusters
                #if envo.G.node[envoID]['daglevel'] > 0 and envoID not in ['ENVO:00002009'] and len(self.studies()) > 1 and self.len() > 10:  ## just "compact clusters with many studies"
                #    hitinfo = Hitinfo(envoID, envoCount, self.len(), len(envo2samples[envoID]), len(self.studies()),
                #                      envo.G.node[envoID]['daglevel'], envo.G.node[envoID]['name'], pvalue, precision, recall, fmeasure)
                hitinfo = Hitinfo(envoID, envoCount, self.len(), len(envo2samples[envoID]), len(self.studies()),
                                  envo.G.node[envoID]['daglevel'], envo.G.node[envoID]['name'], pvalue, precision, recall, fmeasure)
                self.enrichedIn.append(hitinfo) # (pvalue, envoID)                   
        
class Hitinfo:
    """Information container for an enriched cluster that is considered a hit by some of the measures (predominently F-measure).
       Output can be given in HTML as well.
       Bitvector is used to visualize the cluster members in the dendrogram (the bitvector should result in a continuous bar, stretching over the members).
       __cmp__ function to simply do a ranking of Hitinfo objects based on F-measure 
       """
    def __init__(self, envoID, count, clusterSize, envoTotal, supportStudies, daglevel, envoTerm, pvalue, precision, recall, fmeasure):
        self.envoID=envoID
        self.count=count
        self.clusterSize=clusterSize
        self.envoTotal=envoTotal
        self.supportStudies=supportStudies
        self.daglevel=daglevel
        self.envoTerm=envoTerm
        self.logPvalue=-np.log(pvalue) if not (np.isnan(pvalue) or pvalue==0) else np.NaN
        self.pvalue = pvalue
        self.precision=precision
        self.recall=recall
        self.fmeasure=fmeasure
    def addInfo(self, compactness, dunnIndex=0):
        self.homogeneity, self.separation = compactness
        self.dunnIndex = dunnIndex
    def pprint(self, output=sys.stdout, html=False):
        ## 
        formattedString = "%s\t%s\t%s\t%s\t%s\t%02d %-50s\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (self.envoID, self.count, self.clusterSize, self.envoTotal, self.supportStudies, self.daglevel, self.envoTerm[:50], self.precision, self.recall,
                                                                                          self.fmeasure, self.pvalue, self.logPvalue, self.homogeneity, self.separation, self.dunnIndex)
        if not html:
            print >> output, formattedString
        else:
            print >> output, "<TR><TD>%s</TD></TR>" % ("</TD><TD>".join(formattedString.split("\t")))
            
    def generateBitvector(self, length, cluster):
        self.bitvector = np.zeros(length)
        self.bitvector[cluster] = self.fmeasure #self.logPvalue
    def __cmp__(self, other):
        return cmp(self.fmeasure, other.fmeasure)
class Sample:
    """Sample class:
    Describes a sample in terms of study, EnvO annotation, etc. Needed for enrichment analysis and visualization (eg. tooltips) 
    Makes heavy use of SQL information, so make sure you have either a local copy of the database or a fast internet connection and a small dataset, when using the public server remotely.    
    """
    def __init__(self, sampleID, study):
        self.sampleID = sampleID
        if study is None:
            curs.execute("SELECT study FROM samples_unified WHERE sample_event_ID='%s'" % sampleID)
            self.study = curs.fetchone()["study"]
        else:
            self.study = study
        self.envoParentIDs = set()
        self.envoChildrenIDs = set()
        self.generalInfo()
    def generalInfo(self):
        curs.execute("SELECT * FROM samples_EnvO_annotation_unified NATURAL JOIN samples_unified NATURAL JOIN envoColors WHERE sample_event_ID =  '%s'"% self.sampleID)
        results = curs.fetchall()
        self.ecosystem = set([rec["ecosystem"] for rec in results])
        self.isolationSource = set([rec["isolation_source"] for rec in results])
        self.envoTerm = set([rec["OntologyTerm"] for rec in results])
        self.title = set([rec["title"] for rec in results])
        self.tooltip = "<STRONG>%s</STRONG><BR>%s - %s<BR>%s<BR>%s<BR>%s" % (self.sampleID, self.study, "/".join(self.envoTerm), "/".join(self.ecosystem), "/".join(self.title),  "/".join(self.isolationSource))
                     
    def lookupEnvos(self): 
        query = "SELECT OntologyID FROM samples_EnvO_annotation_unified WHERE sample_event_ID='%s' AND ontology='envo'"%(self.sampleID)        
        curs.execute(query)
        result = curs.fetchall()
        if not result:
            print >> sys.stderr, "Warning: no ontology annotation for sample", self.sampleID
        for rec in result:
            envoID = rec["OntologyID"]
            self.addEnvoParents(envoID)
    def colors(self, criteria="color"):
        def expand(color):
            return {'b':'blue'}.get(color, color) ## extend for other single letter abbrevs, so far only b is used
        query = "SELECT color, ecocolor FROM samples_EnvO_annotation_unified NATURAL JOIN envoColors WHERE sample_event_ID='%s' AND ontology='envo'"%(self.sampleID)
        curs.execute(query)
        colors = sorted(list(set([rec[criteria] for rec in curs.fetchall()])))
        if criteria=="ecocolor":
            colors = [cnames[expand(color.lower())] for color in colors]
        lc = len(colors)
        colorVector = np.zeros((60,3))
        if lc:
            colorSegments = np.array_split(colorVector, lc)
            for color, colorSegment in zip(colors, colorSegments):
                colorSegment[:] = hex2color(color)
            return np.vstack(colorSegments).T
        else:
            print "Warning: %s has no assigned colors" % self.sampleID
            return colorVector.T
 
    def addEnvoParents(self, envoID):
        self.envoParentIDs.add(envoID)
        for parent in envo.G.successors(envoID):
            self.addEnvoParents(parent)
    def addEnvoChildren(self, envoID):
        self.envoChildrenIDs.add(envoID)
        for child in envo.G.predecessors(envoID):
            self.addEnvoChildren(child)
    def updateEnvo2samplesDict(self):
        for envoID in self.envoParentIDs:
            envo2samples[envoID].append(self.sampleID)        

class LeveledDAG:
    """class that helps to deal with the depth of a term inside the EnvO multitree, to account for how specific/general a term is.
    Instantiates with a given graph (usually envo.L, where envo is an Ontology object, see envoTools)"""
    def __init__(self, H):
        self.H = H
        self.toplevel = [e for e in H.nodes() if H.out_degree(e)==0]
        for topnode in self.toplevel:
            #self.H.node[topnode]['daglevel'] = 0
            self.addDAGlevels(topnode, 0)
    def addDAGlevels(self, node, level):
        """Recursive procedure for descending top-down into the tree. If a node is encountered that already has an attribute 'daglevel',
        then the smaller value is used"""
        self.H.node[node]['daglevel'] = min(self.H.node[node].get('daglevel', 999999), level)
        for predecessor in self.H.predecessors(node):
            self.addDAGlevels(predecessor, level+1)

def lookupEnvoIdColor(envoID):
    """Every EnvO term has a dedicated color that has been precalculated and stored in the database. See colorEnvograph.py"""
    curs.execute("SELECT color, ecocolor FROM `envoColors` WHERE OntologyID='%s'" % envoID)
    colors = list(set(rec["color"] for rec in curs.fetchall()))
    if colors:
        return colors[0] # just take the first assigned color - there could be more, will have to be ignored for the moment...
    else:
        print "warning: no color found for Envo ID", envoID, envo.G.node[envoID]['name']
        return 'Gray'


def figureEnvoSamples():
    """print out EnvO-sample association statistics:
       How many samples have been associated with each Envo-term.
       Envo-terms are sorted topologically"""
    H = envo.G.subgraph(envo2samples.keys())
    H.add_edge('ENVO:00000894', 'ENVO:00000428') ## retrofit: polar biome is a biome (is obsolete in EnvO, hence disconnected)
    h = nx.topological_sort(H)
    nrSamples, envoTerms = zip(*[(len(envo2samples[e]), envo.G.node[e]['name']) for e in h if len(envo2samples[e])>100])
    fig = plt.figure(figsize=(20,20))
    
    ax = fig.add_subplot(111)
    ax.bar(np.arange(len(envoTerms)), np.log(nrSamples))
    ax.set_xticks(np.arange(len(envoTerms)))
    ax.set_xticklabels(envoTerms, ha='right', rotation=45, size='x-small')
    fig.savefig("/tmp/envo2Samples.png")
    plt.show()

def colorize(vector, colormap="Greys"):# gist_rainbow"):
    import matplotlib.cm as cm
    cmap = cm.get_cmap(colormap)
    cmap.set_bad('w',1.)
    return map(cmap, vector)
    #return ["#%02X%02X%02X" % tuple(255*np.array(colorsys.hsv_to_rgb(i/2., .7, 1))) for i in vector]
def recursiveDendrogramColoring(Y, idx, color, overwrite=False):
    """Coloring branches of the dendrogram"""
    if idx<Y.shape[0] + 1: return
    if dendrogramColors.has_key(idx) and not overwrite: return ## leave sub-cluster in the previous coloring
    dendrogramColors[idx] = color
    cl1, cl2 = map(int, Y[idx - (Y.shape[0] + 1)][:2])
    recursiveDendrogramColoring(Y, cl1, color)
    recursiveDendrogramColoring(Y, cl2, color)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--config", help= "configuration file for database access, see cluster.rc", default="../Data/cluster.rc")
    parser.add_argument("-l", "--linkage", help= "linkage file, npy format, generated by adaptiveRarefaction.py")
    parser.add_argument("-d", "--distancematrix", help= "distance matrix, npy format, generated by adaptiveRarefaction.py")
    parser.add_argument("-n", "--samplenames", help= "list of sample names, must be in same order as in distance matrix, pcl format, generated by adaptiveRarefaction.py")
    parser.add_argument("-e", "--ontodir", help= "data directory that contains EnvO files as generated by oboParser.py", default="../EnvoAnnotation/data/ontology/")
    parser.add_argument("-o", "--outdir", help= "output directory")
    parser.add_argument("--html", help= "output directory", action="store_false")
    #parser.add_argument("-r", "--repetitions", help= "Number of repetitions (multiple adaptive rarefaction)", type=int, default=1)
    args = parser.parse_args()
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    #if not args.outdir: print "No output directory provided"; sys.exit(1)
    #outdir = "/data/Projects/EcoPhyl/Results"

    ## MySQL - read out parameters db, host, user, passwd from a file 
    ## Very simple parsing of MySQL config file, deals with the opt_name=value assignments, syntax exactly as in
    ## MySQL docs: http://dev.mysql.com/doc/refman/5.1/en/option-files.html
    parameterfile = args.config
    if not os.path.exists(parameterfile):
        raise Exception("No MySQL paramater file (default: cluster.rc) provided")
    for line in open(parameterfile):
        if "=" in line and not line.strip().startswith("#"): # skip non-assignment lines, comments
            try:
                exec(line.strip())
            except:
                print "ignoring", line.strip(), "in parameter file", parameterfile
    if not set('database server user password'.split()).issubset(vars().keys()):
        raise Exception ("Not all MySQL access data provided in cluster.rc")
    ## or just hardcode server access in here:
    conn = MySQLdb.connect(db=database, host=server, user=user, passwd=password)
    curs = conn.cursor(DictCursor)
    
    # Load ontology (see EnvO tools)
    ## Ontology files: pickled list of terms
    envo = Ontology("%s/envo5Terms.pcl" % args.ontodir, "%s/envo5.pcl" % args.ontodir, "envo")
    envoL = LeveledDAG(envo.G)
    
    ## For each environment, there was a beta diversity distance matrix produced using scripts in AdaptiveRarefaction
    ## some of them where so large, they had to be calculated on the cluster and then have been copied 
    #envNamesFiles = glob.glob("%s/*_unifrac.pcl" % args.datadir)[2:]
    
    try:
        #sampleNamesFile = "%s/upgma_linkage_samples_%05d_%s_sampleNames.pcl" % (datadir, minSampleSize, environment)        
        ## Create global sample dict (can also be loaded from precalculated pickled files), as well as a dictionary Envo -> list of samples (envo2samples)
        globalSampleDict = {}
        envo2samples = defaultdict(list)
        sampleIDs = dump.load(args.samplenames)
        for sampleID in sampleIDs: 
            s = Sample(sampleID, None)
            s.lookupEnvos()
            s.updateEnvo2samplesDict()
            globalSampleDict[sampleID] = s
        for en, samples in envo2samples.items():
            envo2samples[en] = list(set(samples))        

        # generate colors for studies
        studies = list(set([s.study for s in globalSampleDict.values()]))
        hues = np.linspace(0, 0.9, len(studies))
        htmlColors = ["#%02X%02X%02X" % tuple(255*np.array(colorsys.hsv_to_rgb(i, .7, 1))) for i in hues]
        studyColorDict = dict(zip(studies, htmlColors))
        dendrogramColors = {}

        threshold1, threshold2 = .4, .45 ## Thresholds for dendrogram branch-length to be considered, between 0 and 1,
        ## avoids cluster enrichment test too close to the bottom (i.e. samples too similar, agglomerate them first a bit more)
        ## or too close to top, i.e. too general categories, they are often "enriched" in 'environmental feature' etc, as all samples have that annotation

        ## Read in input:
        ## open previous clustering created with sch.linkage on a comprehensive beta-diversity matrix
        Y = np.load(args.linkage)
        dm = np.load(args.distancematrix)
        names = np.array(dump.load(args.samplenames))
        ## Could do with some sanity checks!
        assert dm.shape[0] == dm.shape[1] == Y.shape[0] + 1 == names.shape[0], Exception("Linkage, Distance Matrix and Sample name list do not match in size! Check input!")
        assert Y.shape[1] == 4, Exception("%s is not a linkage matrix, wrong dimensions" % args.linkage)
        ###################################
        ## Posthoc clustering algorithm  ##
        ###################################
        ## The "heart" of our analysis:
        ## given the scipy linkage matrix, start bottom up to form clusters and check with each merge, whether the newly formed cluster is enriched in an EnvO category
        ## produces "bitvectors" to be placed under the dendrogram that indicates "enriched cluster"-membership (also coloring the dendrogram branches automatically, like in
        ## Lozupone/Knight, PNAS 2007, "Global patterns...", Fig 4
        n = len(Y)+1
        clusterdict = dict([(el, [el]) for el in range(n)]) # the actual clustering works with indices that match the numbering in Y, start with a dictionary, where each element is its own cluster
        hitinfos = {} ## contains results
        hitinfosClassic = [] ## contains results for classic clustering, just looking at compactness
        bitvectors = []
        for cidx, cluster in enumerate(Y):
            if cluster[2] < threshold2:
                cl1, cl2 = map(int, cluster[:2])
                clusterdict[cidx+n] = clusterdict[cl1] + clusterdict[cl2] ## merging subclusters
                if len(clusterdict[cidx+n]) > 10 and threshold1 < cluster[2]: ## every cluster that is tested for enrichment should be of a certain size
                    c = Cluster(names[clusterdict[cidx+n]]) ## creates a cluster object with the sample names of the newly formed cluster           
                    c.analyze(verbose=True) ## Most of the work!
                    if c.enrichedIn:  ## always True if c.mostly = 0.0
                        hitinfo = max(c.enrichedIn)                        
                        hitinfo.addInfo(dunn.compactness(clusterdict[cidx+n], clusterdict.values(), dm), dunn.dunn(clusterdict.values(), dm)) ## dunn index slow! consider to just set to 0
                        hitinfo.generateBitvector(len(names), clusterdict[cidx + n])
                        if args.html:
                            recursiveDendrogramColoring(Y, cidx+n, lookupEnvoIdColor(hitinfo.envoID)) 
                        if hitinfos.has_key(hitinfo.envoID): ## take only the best hit (in terms of f-measure) for a certain envo-ID
                            hitinfos[hitinfo.envoID] = max(hitinfo, hitinfos[hitinfo.envoID])
                        else:
                            hitinfos[hitinfo.envoID] = hitinfo
                            bitvectors.append(hitinfo.envoID)
                        ## always add to classic hitinfo
                        hitinfosClassic.append(hitinfo)

                del clusterdict[cl1]
                del clusterdict[cl2]
                

        ## Generate Output spreadsheet
        outCsv = '%s/enrichedClusters_%s.csv'%(args.outdir, os.path.basename(args.samplenames))
        csv = open(outCsv, "w")
        print >> csv, "EnvoID\tcount\tCluster size\tEnvo Total\tSupported Studies\tDAG-Level Envo-Term\tPrecision\tRecall\tF\tlog P-value (Fisher's)\tHomogeneity\tSeparation\tDunn-Index"
        for envoID in bitvectors:
            hitinfos[envoID].pprint(csv)
        csv.close()
        print "output spread sheet successfully produced:", outCsv
        
        ## Visualization/Output 
        if args.html:
            envoCols = np.transpose(np.dstack([globalSampleDict[name].colors() for name in names]), (1,2,0))  ## 60x|samples|x3 matrix was |samples|x60x3 matrix
            ecoCols  = np.transpose(np.dstack([globalSampleDict[name].colors("ecocolor") for name in names]), (1,2,0))  ## |samples|x60x3 matrix
            studyCols = np.transpose(np.array([[hex2color(studyColorDict[globalSampleDict[names[i]].study])] for i in range(len(names))]), (1,0,2))
            sampleTooltips = [globalSampleDict[name].tooltip for name in names] 
            #annotatedDendrogram(Y, [globalSampleDict[names[i]].study for i in range(len(names))], [cols])
            plt.clf()
            plt.subplot(1,1,1)
            if len(bitvectors)>0:
                bitvectors = np.array([hitinfos[envoID].bitvector for envoID in bitvectors])
                maxvalue = float(bitvectors.max())
                normedBitvectors = 0.999*bitvectors/maxvalue
                order = annotatedDendrogramInfomatrices2(Y, [] , [studyCols, envoCols, ecoCols, normedBitvectors], [None, None, None, "Greys"], [0.05, 0.2, 0.2, 0.4], dendrogramColors)
            else:
                order = annotatedDendrogramInfomatrices2(Y, [] , [studyCols, envoCols, ecoCols], [None, None, None], [0.15, 0.35, 0.35], dendrogramColors)

            image = '%s/environmentDendro_%s.png'%(outdir, os.path.splitext(os.path.basename(args.samplenames))[0])
            imageSVG = os.path.splitext(image)[0]+".svg"
            print "Producing image", image

            dpi = len(names)/12
            plt.savefig(image, dpi=dpi, orientation='portrait')
            plt.savefig(imageSVG, orientation='portrait')
            plt.close()
            ## checking how large the image turned out to be. Hack, works only with imagemagick's identify. Not tested on any other platform than linux!
            xsize, ysize = map(int, os.popen("identify %s"%image).read().split()[2].split("x"))
            generateHTML(os.path.basename(image), np.array(sampleTooltips)[order], xsize, ysize, title="%s (%s)"%(environment, len(names)))
        #plt.show()
    finally:
        conn.close()

