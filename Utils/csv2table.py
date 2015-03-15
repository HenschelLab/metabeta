from matplotlib.mlab import csv2rec
import sys
import MySQLdb
from MySQLdb.cursors import DictCursor
import numpy as np
import nltk
import re

wordPattern = re.compile("[a-zA-Z]{2}")     ## three letter chars required to be considered as a word
nonwordPattern = re.compile("[^a-zA-Z]{3}") ## three non-letter chars -> not a word
stoplist = set([word.strip() for word in open("/home/zain/Projects/NLP/english.stop")])
ps = nltk.PorterStemmer()

def filterWord(word):
    return not word in stoplist and wordPattern.search(word) and not nonwordPattern.search(word) 
    #return wordPattern.search(word) and not nonwordPattern.search(word)

def filterTokens(tokens):
    return [ps.stem(token) for token in tokens if filterWord(token)]
def preprocess(line):
    return nltk.word_tokenize(line.split(";")[0].split(",")[0].lower())

def dtype2SQL(s):
    if s.startswith("<i"):
        return "INT"
    elif s.startswith("<f"):
        return "FLOAT"
    elif s.startswith("|S") and int(s[2:])>255:
        return "TEXT"
    elif s.startswith("|S"):
        return "VARCHAR(255)"
    elif s.startswith("|O4"):
        return "DATE"
    else:
        print "Warning, unknown datatype: %s" % s

def reform(s):
    if "isoformat" in dir(s):
        return s.isoformat()
    elif type(s) == np.float64 and np.isnan(s):
        return "NULL"
    elif not s:
        return "NULL"
    else:
        return s
    
def populate(csvFile, db, table, tablenames, delimiter="\t"):
    data = csv2rec(csvFile, delimiter=delimiter, names=tablenames)

    conn = MySQLdb.connect(db=db, host="localhost", user="zain", passwd="angi4rf")
    curs = conn.cursor(DictCursor)

    cols = ", ".join(["`%s` %s" %(data.dtype.names[i], dtype2SQL(data.dtype[i].str)) for i in range(len(data.dtype.names))])

    curs.execute("DROP TABLE IF EXISTS %s.%s " % (db, table))
    tableQ = "CREATE TABLE %s (id INT NOT NULL AUTO_INCREMENT, %s, PRIMARY KEY (id)) ENGINE = MYISAM" % (table, cols)
    print tableQ
    curs.execute(tableQ)
    for row in data:
        #row = tuple(row) + (" ".join(filterTokens(preprocess(row[-1]))),) ## stoplist/tokenize? -- delete!!!
        insertQ = "INSERT INTO %s VALUES %s" %(table, str(tuple([0]+[str(reform(el)).strip() for el in row])))
        curs.execute(insertQ)
    conn.close()

if __name__ == "__main__":    
    ##sys.argv.append("/home/zain/Projects/OttoTextMining/NcbiTaxonomy/names.dmp")
    csvFile = "/home/zain/Downloads/gg_sample_events_otus.tgz"
    db = "Microbes"
    table = "Mering_sample_events_otus"
    tablenames = "line_number, sample_event_ID, authors_list, title, isolation_source, sequences_count, OTU_TAX_ID".split(", ")
    db, table = sys.argv[2:4]
    populate(sys.argv[1], db, table, tablenames)

