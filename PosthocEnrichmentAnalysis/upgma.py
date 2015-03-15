""" 
Running hierarchical clustering and just saving the 4x(n-1) linkage clustering result
best to submit through ~/Python/qsub:
python qsub.py /home/handreas/Projects/KnowYourEnv/pairwiseUnifrac_heatmap2.py
as it might take longer than the ssh connection is capable of
"""
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import numpy as np
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", help= "Distance Matrix as produced with adaptiveRarefaction. Must be in NumPy (npy) format")
    parser.add_argument("-o", "--output", help= "Name the output file (which is in npy format as well)")
    args = parser.parse_args()
    
    dm = np.load(args.input)
    Y = sch.linkage(squareform(dm), method="average")
    np.save(args.output, Y)
