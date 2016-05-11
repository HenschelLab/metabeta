Script repository and Database dump for 
Microbial Community Meta-analysis

Authors: Andreas Henschel, Muhammad Zohaib Anwar, Vimitha Manohar
Licence: GNU GPL v3.0

Dependencies:

Non-stand Python modules, installable through pip/anaconda
PyCogent (UniFrac, PCoA plots)
NumPy, SciPy, 
Biom
NLTK (phrase similarity, text mining)
MySQLdb
networkx https://networkx.github.io


Most scripts are in heavy development state but come with some documentation. Invoke python scripts with -h to display the associated help page and usage instructions.

Database creation
=================
These scripts are not meant to work anywhere and require local adaptation! Instead we recommend to install the database dump directly:
http://ecophyl.info/html/DataDownload/

populateChaffronDataset.py -> Script is used for populating the database with Chaffron datasets. If needed you can populate your own database with these studies. see http://genome.cshlp.org/content/20/7/947.full
qiimedb2sql.py -> Script is used for populating the Qiime DB studies in to the database. you can use the dump itself or increase the number of studies with this script.
colorEnvo.py -> Script is used for Ecosystem defeinition of samples based on the Envo categories.

BIOM table generation
http://ecophyl.info/html/Biom/

Running scripts that need access: you have two options, in fact three: you can roll your own database, install the database from the SQL dump or access the remote server (default).
Create a file "cluster.rc" and put it in the directory of the python script that requires database access (e.g. posthocEnvoEnrichmentTest.py). 
A sample version that points to the remote server is provided, modify this file if you have a local copy of the database.
It simply contains the parameters database, host, user, passwd and follows the MySQL configuration file syntax: http://dev.mysql.com/doc/refman/5.1/en/option-files.html

EnvO annotation
===============
oboParser.py - simple text parser for obo files, dumps the networkx graph structure and a dictionary EnvoTerms(incl.synonyms) -> 
phraseSimilarity.py - script to calculate the similarity of two phrases as represented by bags of words, using weighted Jaccard index

envoTools.py - graph structure of Env. Ontology, provides functionality to calculate graph distances, parents, children etc.
based on networkx model 

habitatOntologyMapping.py - performs the actual mapping of samples to EnvO terms, input: csv file with sample descriptions, EnvO. See help (command line option -h)

Evaluation scripts for Envo annotation
======================================
We evaluate the accuraccy of EnvO prediction in the following way: we predict the best EnvO term based on weighted Jaccard Index for word phrases (see Phrase Similarity) for non-redundant sample descriptions and their annotations from all.
We then calculate the minimal distance of the predicted and the annotated EnvO terms (given in Qiime-DB as Environmental matter, Environmental feature and Environmental Biome).
The distance is the shortest path in the undirected EnvO Graph, such that exact matches have graph distance 0, EnvO terms in direct subclass-superclass relation have distance 1, direct sibling nodes 2, cousins 4, etc.

It can be seen from the results that from the 712 samples, most automatic annotations are in exact agreement (294) or in direct sub/superclass relation (117) with the manual annotation or
The probabilities for two random nodes to be 0, 1, 2 and 3 steps apart are 0.06\%, 0.20\%, 1.18\% and 3.42\% respectively. Thus 73\% of our annotations are below a 5\% significance level of a random predictor. % reformulate?
Note that manual EnvO annotations are a source of error as well and contribute to disagreement. We corrected for only very few blatant misannotations and consider therefore our accuracy estimates to be a conservative lower bound.

Scripts:

habitatOntologyMappingEvaluation.py - evaluates EnvO annotation using Qiime-DB annotations (subclasses Habitat class from habitatOntologyMapping.py), produces a spreadsheet, see habitatOntologyMappingEvaluation2.csv
habitatOntologyMappingEvaluationHistogram.py - short script that summarizes result table of habitatOntologyMappingEvaluation.py and produces a histogram of graph distances between predicted and manual EnvO terms, corrects for few mistakes in manual annotation.


Adaptive Rarefaction, Beta-diversity calculation and Hierarchical Clustering
============================================================================

adaptiveRarefaction.py, see script help (option -h)

Precalculated Beta-diversity distance matrices are provided in numpy matrix format (together with a pickle file containing the sample identifiers in the same order as the matrix) in the Data directory, 
as they are very computationally expensive (especially when using adaptive rarefaction), see adaptiveRarefaction.py


Post-hoc enrichment test for EnvO
=================================
upgma.py takes as input a beta-diversity distance matrix from adaptiveRarefaction.py and creates a dendrogram based on scipy's linkage function. Stores result in a (n-1) x 4 matrix, also numpy format.
posthocEnvoEnrichmentTest.py - This is the most central script in the set! It tests whether clusters are enriched in some EnvO term. Enrichment tests include parent categories. Input: clustering result from upgma.py, the distance matrix from adaptiveRarefaction.py, the list of sample names used (required to be in same order as in the distance matrix), EnvO etc. Produces a spreadsheet output and on demand an HTML visualization (requires ImageMagick and JavaScript/Style sheets as provided in the example HTML output). 

In order to produce Radio buttons for coloring of PCoA plots by EnvO, Ecosystem and study, we provide the modified qiime lib file (not the bin file with the same name!):
make_2d_plots.py

It needs to replace the corresponding file in QIIME. Note, we use the QIIME library version 1.8.0.

Enjoy!
Contact: ahenschel@masdar.ac.ae
