"""
Data Downloaded:
Qiime-db from microbio.me/emp, the 'official' link:
SRA from ncbi.nlm.nih.gov/sra, the 'official' link:

A bunch of Biom tables (OTU vs sample matrices with Meta data per sample (often empty) and per OTU (lineage))
Also comes with a mapping table, which contains EnvO annotation (make sure EnvO versions fit!)

Study names in the database are given with data source prefix. (QDB=QiimeDb, SRA=SRA) 

Parsing biom tables and EnvO information from mapping, populate SQL database:

sample data for each instance of iterTable -> Table samples_Qiime (create Identifier)

Functions can be used in main depending upon the data format
"""

from biom.parse import parse_biom_table
import MySQLdb
from MySQLdb.cursors import DictCursor
from envoTools import Ontology
import numpy as np
import glob
import pdb

manfix = {'egg':"hen's egg product", 'organic material feature':'organic material', 'Bay': 'bay', 'rainforest division (420)': 'Rainforest Division (420)', 'Cheese':'cheese product'}
def envoLookup(term):
    if not term.startswith("ENVO:"):
        print "Warning: No Envo term provided", term
        return None, None
    term = " ".join(term.split("ENVO:")[1].split())
    term = manfix.get(term, term)
    result = envo.synonymDict[term]
    if not result:
        #term = " ".join([word.capitalize() for word in term.split()])
        term = term.capitalize()
        result = envo.synonymDict[term]
        if not result:
            print "Warning, envo term not found:", term
            return (None, None)
    return (term, result[0])

def populate_samples(mappingfile):
    source="QDB" #SRA
    study = source+"_%04d" % int(mappingfile.split("/")[-1].split("_")[1])
    
    for line in open(mappingfile):
        if line.strip().startswith("#"):
            continue
        fields = line.strip().split("\t")
        sql = 'INSERT INTO samples_qdb (sample_event_ID, study) VALUES ("%s", "%s")' % (fields[0], study)
        curs.execute(sql)
        
def populate_sample_EnvO_annotation(mappingfile):
    for line in open(mappingfile):
        if line.strip().startswith("#"):
            cols = line.strip().split("\t")
            envoIdx = [idx for idx, col in enumerate(cols) if col.startswith("ENV_")]
            continue
        fields = line.strip().split("\t")
        envoFields = [envoLookup(fields[envoCol]) for envoCol in envoIdx]
        sampleID = fields[0] ## Should this be changed?
        #sqld ='DELETE FROM samples_EnvO_annotation_qdb WHERE sample_event_ID="%s"'% sampleID
        #nr = curs.execute(sqld)
        #print "Deleting %s entries from" % nr, sampleID
        try:
            for term, envoID in envoFields:
                if envoID is not None:
                    sql = 'INSERT INTO samples_EnvO_annotation_qdb VALUES ("%s", "%s", "%s", "envo", "", 1)' % (sampleID, term, envoID)
                    curs.execute(sql)
        except:
            pdb.set_trace()
            
def populate_OTUs_samples(otuTable):
    for sampleData, sampleID, m in otuTable.iterSamples():
        for idx, abundance in enumerate(sampleData):
            if abundance:
                curs.execute('INSERT INTO OTUs_samples_qdb VALUES ("%s", "%s", %s)' % ( otuTable.ObservationIds[idx], sampleID, abundance))
                
def populate_OTUs(otuTable):
    for otuData, otu, meta in otuTable.iterObservations():
        curs.execute('SELECT * FROM OTUS_qdb WHERE otu_id = "%s"' % otu)
        if not curs.fetchone():            
            curs.execute('INSERT INTO OTUS_qdb VALUES ("%s", "%s", "")' % ( otu, "; ".join(meta.get('taxonomy', []))))

def fix_samples_studies():
    curs.execute('SELECT title FROM GG13_samples GROUP BY title')
    for idx, rec in enumerate(curs.fetchall()): #UPDATE  `Microbes`.`samples` SET  `study` =  'CHA_0001' WHERE
        curs.execute('UPDATE GG13_samples SET study = "SRA_%04d" WHERE title="%s"'%(idx, rec["title"])) ### TO BE CONTINUED!!!

def fix_samples_studies2():
    curs.execute('SELECT * FROM samples NATURAL JOIN samples_EnvO_annotation WHERE isolation_source LIKE  "%hypersaline%"')
    for rec in curs.fetchall():
        curs.execute('INSERT INTO samples_EnvO_annotation VALUES ("%s", "ID:0000" )' % rec["samples_event_ID"]) #sample_event_ID OntologyTerm OntologyID Ontology OntologyLineage Score

if __name__ == "__main__":
    ## Provide your local database info here!
    conn = MySQLdb.connect(db="GlobalMicroBiome", host="", user="", passwd="")
    curs = conn.cursor(DictCursor)
    fix_samples_studies2()
    
    """datadir = "/data/EarthMicrobiomeProject/TGZ/TODO/*/"
    ontoDir = "/home/zain/Projects/OttoTextMining/OntologyData"
    envo = Ontology("%s/envoTerms3.pcl" % ontoDir, "%s/envo3.pcl" % ontoDir, "envo")
    #curs.execute("TRUNCATE TABLE samples_qdb")
    #curs.execute("TRUNCATE TABLE OTUs_samples_qdb")
    #curs.execute("TRUNCATE TABLE OTUS_qdb")
    for subdir in glob.glob(datadir)[28:]:
        if glob.glob("/data/EarthMicrobiomeProject/study_%s" % subdir.rstrip("/").split("/")[-1].split("_")[1]):
            print "subdir %s exists, skipping" % subdir
            continue
        print "##############", subdir, "####################"
        mapping = glob.glob("%s*mapping_file.txt" % subdir)
        biom = glob.glob("%s*.biom" % subdir)
        if not biom or not mapping: continue
        biomfile = biom[0]
        mappingfile = mapping[0]

        otuTable = parse_biom_table(open(biomfile, 'U'))    
        populate_OTUs_samples(otuTable)
        populate_OTUs(otuTable)
        populate_samples(mappingfile)
        populate_sample_EnvO_annotation(mappingfile)
        #break"""
    conn.close()
    

