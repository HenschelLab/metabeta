"""
Script to evaluate habitat annotation with EnvO terms
Author: Andreas Henschel
Licence: GPL v3.0

For internal use to test how well mapping works.
Compares the best matches to manually assigned EnvO terms in Qiime-DB (considered 'Gold standard', obviously not perfect)
python habitatOntologyMappingEvaluation.py > habitatOntologyMappingEvaluation.csv

For each provided habitat description, the tool finds the best ontology term, based on weighted Jaccard index similarity.

Environment ontology, pre-parsed and saved as a networkx graph,
see oboParser.py (which produces envo?.pcl and envoTerms?.pcl)

After preprocessing, all habitat terms are matched against the
ontology terms (incl. their synonyms) using something Jaccard index
for bags of words: bestmatch(h) = argmax o_i Jaccard(h,o_i) where h is
the habitat term and o_i are the ontology terms.

The word phrases are filtered with a stop list (from NLTK) and
weighted based on their occurrence in the Brown corpus as well as
their appearance in habitat descriptions, thus calibrating phrase
similarity for use with IDF-Weighted Jaccard similarity (Weighted
Set-Based String Similarity, Hadjieleftheriou and Srivastava)

We evaluate the accuraccy of EnvO prediction in the following way: we predict the best EnvO term based on weighted Jaccard Index for word phrases (see Phrase Similarity) for non-redundant sample descriptions and their annotations from all.
We then calculate the minimal distance of the predicted and the annotated EnvO terms (given in Qiime-DB as Environmental matter, Environmental feature and Environmental Biome).
The distance is the shortest path in the undirected EnvO Graph, such that exact matches have graph distance 0, EnvO terms in direct subclass-superclass relation have distance 1, direct sibling nodes 2, cousins 4, etc.
The results are shown in Fig \ref{envoAnnotationEval}. It can be seen that from the 712 samples, most automatic annotations are in exact agreement (294) or in direct sub/superclass relation (117) with the manual annotation or
The probabilities for two random nodes to be 0, 1, 2 and 3 steps apart are 0.06\%, 0.20\%, 1.18\% and 3.42\% respectively. Thus 73\% of our annotations are below a 5\% significance level of a random predictor. % reformulate?
Note that manual EnvO annotations are a source of error as well and contribute to disagreement. We corrected for only very few blatant misannotations and consider therefore our accuracy estimates to be a conservative lower bound.
%It includes several semantically unlinked classification systems (including WWF, Udvardy and Baileys biomes) that can cause misleadingly graph distances despite strong semantic relations.
%Further, submitted annotations to Qiime-DB are also a source of error.

"""
import re, sys, pdb
import math
from nltk.corpus import brown
from nltk.probability import *
from itertools import chain
from dump import dump
from envoTools import Ontology
from habitatOntologyMapping import *
import MySQLdb
from MySQLdb.cursors import DictCursor

class VerifiableHabitatDescription(HabitatDescription):
    def __init__(self, eventID, habitat, study):
        self.goldStandardEnvos = set([])
        self.study = study
        HabitatDescription.__init__(self, eventID, habitat)
    def addEnvoTerm(self, envoTerm):
        self.goldStandardEnvos.add(envoTerm)
    def verifyEnvo(self, envo):
        (self.distance, self.pvalue), self.bestMatch = min([(envo.shortestPath(self.annotations[0][4], envoTerm), envoTerm) for envoTerm in self.goldStandardEnvos])
    def csv(self):
        goldStandardString = ["%s\t%s"%(envoTerm, envo.G.node[envoTerm]['name']) for envoTerm in self.goldStandardEnvos if envoTerm in envo.G.node.keys()]
        if not goldStandardString:
            return
        predictedTerm = envo.G.node[self.annotations[0][4]]['name'] if self.annotations[0][4] else ""
        st = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t"% (self.sampleEventID, self.habitatDescription,self.study,
                                                      self.annotations[0][4],
                                                      predictedTerm,
                                                      self.annotations[0][2],                                                
                                                      self.distance, self.pvalue, self.bestMatch) + "\t".join(goldStandardString)
                                                      
        #if st.count("\t") in [1,10,12]:
        #    pdb.set_trace()
        return st
    def __eq__(self, other):
        return self.habitatDescription == other.habitatDescription and self.goldStandardEnvos == other.goldStandardEnvos
    def __hash__(self):
        ## to be non-redundant, allow a set reduction that considers samples with equal textual description and equal annotated Envo Terms (Gold standard) as equivalent 
        return hash((self.habitatDescription, tuple(sorted(list(self.goldStandardEnvos)))))
    
if __name__ == "__main__":
    ## TODO: make this universal, so it can be run through the amazon cloud server
    conn = MySQLdb.connect(db="ServerMicroBiome", host="cis1-db", user="ahenschel", passwd="angi4rf")
    curs = conn.cursor(DictCursor)
    
    ## Create Ontology
    ontoDir = "data/ontology"
    envo = Ontology("%s/envo5Terms.pcl" % ontoDir, "%s/envo5.pcl" % ontoDir, "envo")

    envo.initPhraseSimilarity()
    envo.initShortestPath()

    #Read in habitat descriptions
    ## Select verifiable samples from QiimeDB
    #habitats = [HabitatDescription(*lineparse(line, args)) for line in open(args.mapfile) if not line.lstrip().startswith("#")]
    sampleQuery = """SELECT sample_event_ID, isolation_source, title, study, OntologyTerm, OntologyID
                     FROM samples_unified NATURAL JOIN `samples_EnvO_annotation_unified`
                     WHERE study LIKE 'QDB%'"""
    curs.execute(sampleQuery)
    verifiableHabitats = {}
    for rec in curs.fetchall():        
        sampleID = rec['sample_event_ID']
        if not sampleID in verifiableHabitats:
            habitatdescription = "%s, %s" % (rec['title'].strip(), rec['isolation_source'].strip())
            verifiableHabitats[sampleID] = VerifiableHabitatDescription(rec['sample_event_ID'], habitatdescription.replace("_", " "), rec['study'])
        verifiableHabitats[sampleID].addEnvoTerm(rec['OntologyID'])
    conn.close()
    ## Header for tab-sep spreadsheet output
    print "Sample ID	Sample Description (title, isolation source)	Study ID	Predicted EnvoID	Predicted Envo Term	Jaccard Score	Distance	P-value	Best Match	Annotated EnvoID1	Predicted Envo Term1	Annotated EnvoID2	Predicted Envo Term2	Annotated EnvoID3	Predicted Envo Term3"
    for habitat in set(verifiableHabitats.values()):
        envo.findBestMatch(habitat)
        habitat.verifyEnvo(envo)
        csv = habitat.csv()
        if csv: print csv
        
        
