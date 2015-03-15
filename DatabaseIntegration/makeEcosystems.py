"""
Distributed under GNU V3 lisence

The script works with Envo file from environmentontology.org/downloads.
Envo file is used with the custom definition of Ecosystems based on EnvO terms. 

Requirements:
htmlcols.py from utils
envoTools.py from utils
Envo file from utils
envo5Terms.pcl file from utils

"""
import sys
import networkx as nx
import MySQLdb
from MySQLdb.cursors import DictCursor
from envoTools import Ontology
from collections import defaultdict
import itertools
from htmlcols import *
import numpy as np



def mapEcosystems2envos(allEnvos):
    """
    produces a dict:
    {'Biofilm': [sampleID1, sampleID2 ...], 'Soil': [sampleID1, ...]}
    """
    ecosystems2samples = defaultdict(list)
    for idx, ecosystem in enumerate(ecosystems):
        obsEnvosInEcosystem = set()
        for ecosystemCat in categories[ecosystem]:            
            if envo.synonymDict.has_key(ecosystemCat):
                envoClass = envo.synonymDict[ecosystemCat][0]
                subsumedEnvos = [en for en in allEnvos if nx.has_path(envo.G, en, envoClass)]
                obsEnvosInEcosystem = obsEnvosInEcosystem.union(subsumedEnvos)
            else:
                print >> sys.stderr, "Warning: ignoring toplevel ecosystem branch %s " % ecosystemCat        
        ecosystems2samples[ecosystem] = obsEnvosInEcosystem
    return ecosystems2samples



if __name__ == "__main__":
    categories = {'Biofilm': ['biofilm', 'microbial mat material', 'biofilm material', 'microbial mat'],
                  'Anthropogenic': ['anthropogenic feature', 'anthropogenic abiotic mesoscopic feature', 'anthropogenic environmental material', 'anthropogenic habitat', 'bioreactor', 'biofilter'],
                  'Marine': ['marsh', 'marine biome', 'marine sediment', 'saline water', 'coastal inlet', 'marine water body', 'saline hydrographic feature', 'coast', 'archipelago', 'marine channel', 'seashore', 'reef', 'undersea feature', 'black smoker', 'marine feature', 'marine snow', 'coastal wetland', 'saline water habitat'],
                  'Freshwater': ['freshwater wetland', 'glacial feature', 'reservoir', 'freshwater habitat', 'aquifer', 'fresh water', 'freshwater lake', 'freshwater biome'],
                  'Soil': ['soil', 'mountain', 'mountain range', 'karst', 'terrestrial biome', 'plantation', 'mud', 'depression', 'pebble sediment', 'clay', 'terrestrial habitat', 'sandy sediment', 'landslide', 'beach', 'desert', 'subterrestrial habitat', 'sediment'],
                  'Animal/Human': ['bodily fluid', 'animal food product', 'animal-associated habitat'],
                  'Plant': ['plantation', 'plant-associated habitat', 'plant food product', 'rhizosphere'],
                  'Geothermal': ['volcanic feature', 'geothermal power plant', 'volcanic hydrographic feature'],
                  'Hypersaline': ['haline habitat', 'hypersaline'],
                  }

    colors = {'Plant': 'DarkGreen', 'Geothermal': 'SaddleBrown', 'Soil': 'Gold', 'Biofilm': 'SlateGray', 'Animal/Human': 'DarkViolet', 'Freshwater': 'b', 'Marine': 'Cyan', 'Anthropogenic': 'DarkOrange', 'Hypersaline':'AliceBlue'}
    ecosystems = categories.keys()
    
    #Database details
    conn = MySQLdb.connect(db="#YourDatabase", host="", user="", passwd="")
    curs = conn.cursor(DictCursor)

    #EnvO File from environmentontology.org/downloads
    envo = Ontology("envoTerms.pcl", "envo.pcl", "envo")

    curs.execute("SELECT DISTINCT(OntologyID) FROM samples_EnvO_annotation_unified WHERE Ontology='envo'")
    allEnvos = [rec["OntologyID"] for rec in curs.fetchall()]
    ecosystems2envos = mapEcosystems2envos(allEnvos)
    allCatEnvos = set(itertools.chain.from_iterable(ecosystems2envos.values()))
    uncatEnvos = set(allEnvos).difference(allCatEnvos)
    try:
        for ecosystem, envos in ecosystems2envos.items():
            ecocol = colors[ecosystem]
            for envo in envos:
		length=len(htmlcols)
                envocol = htmlcols[np.random.randint(length)]
		print	"INSERT INTO envoColors_test (OntologyID, ecosystem, ecocolor, color) VALUES ('%s', '%s', '%s', '%s')" % (envo, ecosystem, ecocol, envocol)
                curs.execute("INSERT INTO envoColors_test (OntologyID, ecosystem, ecocolor, color) VALUES ('%s', '%s', '%s', '%s')" % (envo, ecosystem, ecocol, envocol))
    finally:
        conn.close()
    
            
