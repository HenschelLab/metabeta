"""
Populating Chaffron dataset, available at
Requires that the data is downloaded from
http://mblnx-kallisto.uzh.ch:8888/microbial_coexistence/

please set the path (datadir) to the directory where the unpacked Chaffron data set resides!
"""

from csv2table import populate
import glob

datadir = "../Data/" ## modify this!
db = "ServerMicroBiome"
tablenames = "line_number, sample_event_ID, authors_list, title, isolation_source, sequences_count, OTU_TAX_ID".split(", ")

for csvFile in glob.glob("%s/gg_sample_details_otus_filtered_file.*"):
    basename, t1, distance, t2  = csvFile.split("/")[-1].split(".")
    if distance == "01": continue
    table = "_".join(basename.split("_")[1:4] + [distance])
    populate(csvFile, db, table, tablenames)
