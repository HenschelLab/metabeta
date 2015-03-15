import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from itertools import groupby

def process(line):
    fields = line.split("\t")
    sampleID = fields[0]
    studyID = fields[2]
    jaccard = float(fields[5])
    ontoDistance = int(fields[6])
    ## manual correction, see justification below!
    if studyID == "QDB_0317" and fields[4] == "hand" and fields[10] == "sebum":
        ontoDistance = 0
    return (ontoDistance, jaccard)

spreadsheet = "/home/zain/Projects/KnowYourEnv/GitRep/habitatOntologyMappingEvaluation2.csv"
distanceJaccard = [process(line) for line in open(spreadsheet) if not line.startswith("Sample ID")]
distanceJaccard.sort()

d,j = zip(*distanceJaccard)
dc = Counter(d)
#disthist, edges = np.histogram(d, bins=max(d))
positions = np.arange(max(d)+1) - 0.5

plt.close('all')
plt.bar(positions, [dc[i] for i in range(max(d)+1)])
plt.xticks(np.arange(max(d)+1))
plt.ylabel("Frequency")
plt.xlabel("EnvO Graph Distance between automated and manual annotation")
#f, axarr = plt.subplots(2, sharex=True)
#axarr[0].bar(positions, [dc[i] for i in range(max(d)+1)])

## Boxplots:
#data0 = dict([(distance, zip(*list(jscores))[1]) for distance, jscores in groupby(distanceJaccard, lambda x: x[0])])
#data = [data0.get(i,[]) for i in range(max(d)+1)]
#axarr[1].boxplot(data)    
#axarr[1].set_xticklabels(np.arange(max(d)+1))
plt.show()


## Differences in sebum production are postulated to be a potential sources as to why women and men differ in microbiomes hand palm. "Hand" (Envo ID:0000023)
## Fierer N, Hamady M, Lauber CL, Knight R. The influence of sex, handedness, and washing on the diversity of hand surface bacteria. Proceedings of the National Academy of Sciences of the United States of America. 2008;105(46):17994-17999. doi:10.1073/pnas.0807920105.
## we suspect further errors to originate from inaccurate author annotation in Qiime DB but we will for the sake of simplicity avoid manual correction and consider the provided accuracy of our automated EnvO annotation algorithm as a lower bound.
## We also omit results from study 2013, for which we can not retrieve the associated publication due to a now disfunct Earth Microbiome Database.
