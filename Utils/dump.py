import cPickle

def dump(data, filename="/tmp/data.pcl"):
    datafile = open(filename, "a")
    cPickle.dump(data, datafile)
    datafile.close()

def load(filename="/tmp/data.pcl"):
    datafile = open(filename)
    data = cPickle.load(datafile)
    datafile.close()
    return data
