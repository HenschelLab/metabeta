"""
Tool to annotate habitat descriptions with EnvO terms
Author: Andreas Henschel
Licence: GPL v3.0

Sample usage:
python habitatOntologyMapping.py data/sample_details.txt -c 2 > sample_details.html
python habitatOntologyMapping.py data/sample_details.txt -c 2 --tab > example/sample_details2.csv

For each provided habitat description, the tool finds the best ontology term, based on weighted Jaccard index similarity.
The output html page highlights the words that the match is based on.
For the best ontology maps, the page provides the ancestry of Ontology terms as mouse-over tooltip.
(make sure you have tooltips.js and style.css in the same directory as the produced html file)

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
"""
import re, sys
import math
import nltk
from nltk.corpus import brown
from nltk.probability import *
from itertools import chain
from dump import dump
from envoTools import Ontology
import argparse

wordPattern = re.compile("[a-zA-Z]{2}")     ## three letter chars required to be considered as a word
nonwordPattern = re.compile("[^a-zA-Z]{3}") ## three non-letter chars -> not a word
stoplist = set([word.strip() for word in open("data/english.stop")])
habitatIdCounter = 0

def filterWord(word):
    return not word in stoplist and wordPattern.search(word) and not nonwordPattern.search(word) 
    #return wordPattern.search(word) and not nonwordPattern.search(word)

def filterTokens(tokens):
    return [token for token in tokens if filterWord(token)]
def preprocess(line):
    return nltk.word_tokenize(line.split(";")[0].lower())
def idf(f, D):
    return math.log(1 + (D/f))

html1 = '''
<html>
  <head>
    <script type="text/javascript" language="javascript" src="tooltip.js"></script>
    <link rel="stylesheet" type="text/css" href="style.css" />
    <style type="text/css">
      .envo {color: #009000;}
      .mesh {color: #900000;}
      table {font-family: arial,sans-serif;font-size: x-small;} th {background: #336699; color: #FFFFFF; text-align: left;}
    </style>
  </head>
  <body>Below is a list of habitat descriptions and their best matches to Ontology terms. Mouse over ontology terms to see more categorization details (parent classes of the ontology term). Best match is based on text mining methods (weighted Jaccard index for bags of words). The score is between 0 and 1, 0 indicating no similarity and 1 indicating a perfect match.'''
html2 = "</body></html>"
tableHeader = '''<TABLE cellspacing="0" cellpadding="4" border="1" bordercolor="#224466" width="100%"><TR><TH>Habitat description</TH><TH>ENVO Ontology Terms</TH><TH>Jaccard</TH><TH>MESH Ontology Terms</TH><TH>Jaccard</TH></TR>'''
tableHeaderCSV = ['SampleEventID', 'Habitat description', 'EnvO Ontology Term', 'EnvO ID', 'Jaccard']
tableFooter = '</TABLE>'
headline = "<H2> %s </H2>"

class HabitatDescription:
    def __init__(self, eventID, habitat):
        global habitatIdCounter        
        #self.sampleEventID, self.habitatDescription = habitatData.split("\t")
        self.sampleEventID = eventID
        self.habitatDescription = habitat
        self.bagOfWords = set(filterTokens(preprocess(self.habitatDescription)))
        self.tokens = nltk.word_tokenize(self.habitatDescription)
        self.tokenAnnotations = [[] for i in self.tokens]
        self.annotations = []
        self.id = habitatIdCounter
        habitatIdCounter += 1        
    def annotate(self, ontoTerm, ontoClass): ## could be simplified! 
        self.annotations.append(ontoClass)
        if not ontoTerm: return
        for i in range(len(self.tokens)):
            word = self.tokens[i].lower()
            if word in ontoTerm.lower():
                self.tokenAnnotations[i].append(ontoClass)
    def htmlAnnotation(self):
        ## result could be simplified
        ## TODO: handle multiple annotations per word
        html = "<TR><TD>"
        for token, annotation in zip(self.tokens, self.tokenAnnotations):
            if annotation:
                html += ' <span class="%s">%s</span> ' % (annotation[0][0], token)                
            else:
                html += " %s" % token

        html += "</TD>" + "".join(['<TD><span class="info" onmouseover="tooltip.show(\'%s\');" onmouseout="tooltip.hide()"> %s/%s </span></TD><TD>%.3f</TD>' % (a[3], a[1], a[4], a[2]) for a in self.annotations])
        if len(self.annotations) == 1: html += "<TD>-</TD> <TD>n/a </TD>" # mocking up missing MeSH annotation, quick and dirty
        html += "</TR>"
        self.html = html
        self.bestScore = reduce(lambda x,y: x*y, [a[2] for a in self.annotations]) ## used to be max ...
        a = self.annotations[0]
        self.csvList = [self.sampleEventID, self.habitatDescription, a[1], a[4], str(a[2])]
    def __cmp__(self, hb):
        return cmp(self.bestScore, hb.bestScore)
    
def lineparse(line, options):
    fields = line.rstrip().split(options.delimiter)
    return fields[0], fields[options.column+1]
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mapfile", help="map file containing habitat descriptions, such as isolation source, experiment description or publication title/abstract. Can deal with MIMARKS/Qiime map file format http://qiime.org/documentation/file_formats.html (at least put sample id in column 1 and specify the description column)")
    parser.add_argument("-c", "--column", help="specify column number, default: last column", type=int, default=-1) ## read description from last column by default
    parser.add_argument("-d", "--delimiter", help="specify delimiter in text file, if it is not TAB", default="\t")
    parser.add_argument("-g", "--html", help="generate html", action="store_true", default=False)
    parser.add_argument("-t", "--tab", help="generate tab separated csv output", action="store_false", default=False)
    #parser.add_argument("-e", "--ecosystem", help="generate html for each ecosystem category", action="store_true", default=False)
    args = parser.parse_args()
    
## Count word freqs from brown but also from our environment articles (title + habitat descriptions)
    words = FreqDist()
    for sentence in brown.sents():
        for word in sentence:
            words.inc(word.lower())
    for line in open("data/sample_details.txt"):
        for token in nltk.word_tokenize(" ".join(line.split("\t")[2:])):
            words.inc(token.lower())
            
    D = float(sum(words.values()))
    maxIdf = idf(0.5, D)
    weights = dict([(w, idf(f, D)) for (w,f) in words.items()])


    #Read in habitat descriptions
    habitats = [HabitatDescription(*lineparse(line, args)) for line in open(args.mapfile) if not line.lstrip().startswith("#")]

    ## Create Ontology
    ontoDir = "data/ontology"
    envo = Ontology("%s/envo5Terms.pcl" % ontoDir, "%s/envo5.pcl" % ontoDir, "envo")
    envo.initPhraseSimilarity()
    for habitat in habitats:
        envo.findBestMatch(habitat)

    for habitat in habitats:
        habitat.htmlAnnotation()
    habitats.sort()

    if args.html and not args.tab:
        print html1
        #category = 'Gut' ## generalize!
        #print headline % category
        print tableHeader
        #envo.listHabitats(category.lower())
        for habitat in habitats[::-1]:
            print habitat.html
        print tableFooter
        print html2
    else:
        for habitat in habitats[::-1]:
            print args.delimiter.join(habitat.csvList)
    
## cluster habitats by ontology: get all pairwise distances -> distance matrix -> hier. cluster

