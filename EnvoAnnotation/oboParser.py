"""
Simple parsing of obo files (in particular the OBO version of EnvO, obtained from obofoundry.org)

Author: Andreas Henschel
Licence: GPL v3.0

Usage:
python oboParser.py data/ontology/envo5.obo

Note that envo5 is a manually curated version of envo available at obofoundry.org.
According to Buttigieg et al, 2013, check also:
http://purl.obolibrary.org/obo/envo.obo

Saves 2 files: the parsed graph structure (networkx) -> envo5.pcl
and a dictionary with all synonyms -> envo5Terms.pcl (pickled file)
Actually quite fast, consider including this in envoTools, so no need to run it separately every time, envo is updated
"""
import networkx as nx
import cPickle
from collections import defaultdict
import sys, os
import argparse

parser = argparse.ArgumentParser(description = __doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("obofile", help="obo ontology file, as can be acquired from obofoundry.org")
args = parser.parse_args()
obofile = args.obofile

#if len(sys.argv)>1:
#    obofile = sys.argv[1]
basename, extension = os.path.splitext(obofile)
#synonymDictFile = basename + "_synonymDict.pcl"
#envoGraphFile   = basename + "_1.pcl"

if os.path.exists(obofile):
    synonymDictFile = basename + "Terms.pcl"
    envoGraphFile   = basename + ".pcl"
else:
    print "input file '%s' not found, exiting..." %obofile
    sys.exit(1)
    
print "Output will go to:", envoGraphFile
tags = "id name abstract def is_a relationship synonym".split() #: part_of
parsing=False
termDicts = {}
synonymDict = defaultdict(list)
for line in open(obofile):
    if line.startswith("["):
        if line.startswith("[Term]"):
            parsing=True
            termDict = defaultdict(list)
        else:
            parsing=False
    elif parsing and line.strip():
        line1 = line.split("!")[0].split(":")
        tag, value =  line1[0], ":".join(line1[1:]).strip()
        if tag in tags:
            if tag in ["synonym"] and '"' in value: ## lists for these attributes
                synonym = value.split('"')[1] #re.findall('\"\w*\"', value)[0][1:-1]
                termDict[tag].append(value)
                ## ignore synonyms provided for abstract terms
                if termDict['abstract'] != 'True':
                    synonymDict[synonym].append(termDict['id'])
                
            elif tag in ["relationship", 'is_a']:
                if tag == "relationship":
                    tag, value = value.split()
                termDict["relation"].append((tag, value))
            else:
                termDict[tag] = value
                if termDict['abstract'] != 'True' and tag == 'name':
                    synonymDict[value].append(termDict['id'])
                    
    if parsing and termDict and not line.strip():
        termDicts[termDict['id']] = termDict

##
f = open(synonymDictFile, "w")
cPickle.dump(synonymDict, f)
f.close()

G = nx.DiGraph()
## Add nodes
for id, termDict in termDicts.items():
    G.add_node(id, attr_dict=termDict)
## Add Edges
for node in G.nodes():
    if G.node[node].has_key('is_a'):
        superCategories = G.node[node]['is_a']
        for isa in superCategories:
            G.add_edge(node, isa, rel='is_a')
    if G.node[node].has_key('relation'):
        relations = G.node[node]['relation']
        for rel, rnode  in relations:
            G.add_edge(node, rnode, rel=rel)
nx.write_gpickle(G, envoGraphFile)

## Subgraph: sg = nx.DiGraph([ (u,v,d) for u,v,d in G.edges(data=True) if d['rel']=='is_a'])
