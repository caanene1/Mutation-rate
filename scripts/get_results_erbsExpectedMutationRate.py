#!/usr/bin/env python3 

# script to collect randomised mutations in a metafile
# written by Sabarinathan Radhakrishnan (Sabarinathan R et al. Nature 2016, 532:264-267)
# modified by Jiekun (Jackie) Yang

import os
import sys
import numpy as np
import pandas as pd
import glob
import scipy as sp
import scipy.stats
import gzip
import argparse
import logging

# paths
path=os.environ['path']

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-w', dest='WGSused', help='WGS data used (e.g brcaeu)')
parser.add_argument('-e', dest='ERBSused', help='ERBS data set used (e.g allERBS)')
parser.add_argument('-o', dest='outdir', help='output directory')
parser.add_argument('--cores', dest='cores', type=int, default=os.cpu_count(), help='Maximum number of CPU cores to use')
parser.add_argument('--debug', dest='debug', action='store_true')
args = parser.parse_args()

if args.WGSused == None or args.ERBSused == None or args.outdir == None:
    parser.print_usage()
    sys.exit(1)

# set output path
opath = "%s/results/%s" % (path, args.outdir) # output path

# default params
rand=1000 # number of random sampling performed 
flank=1000 # flanking from the ERBS mid-point 

### functions
def mean_confidence_interval(line):
    confidence=0.95
   # data = [line['rand_' + str(i)] for i in range(1, (rand+1))]
    data = line
    n = len(data)
    m, se, sd = np.mean(data), sp.stats.sem(data), np.std(data)
    h = se * sp.stats.t._ppf((1 + confidence) / 2, n-1) # confidence interval
    return m, m - h, m + h, sd

def extractDetails(sfile):
    # choose the observed mutation rate file to get the number of binding sites in each position for normalization
    pos={}
    er = pd.read_csv(sfile, sep="\t", header=0)
    # create the dictionary to save the number of binding sites for each position
    erBS = dict(zip(er['position'],er['bp']))
 
    # expectedMutRate file name
    mfilename="%s/results/%s/%s/expected/expectedMutRate_rand.bin.gz" % (path, args.WGSused, args.ERBSused)
    # unzip and read the expectedMutRate file
    a = np.fromstring(gzip.open(mfilename, "rb").read(), dtype="int32")

    finalDic={}
    #compute the average for each random sampling
    for sam in range(0, rand):
        # intialize array to store the mutation position counts for each sampling
        for i in range(-flank, flank+1):
            pos[i]=0

        # extract the output for one sampling
        start = int((len(a)/rand)*sam)
        count = np.unique(a[start:(start+int(len(a)/rand))], return_counts=True) # count for one random sampling
        for i in range(0, len(count[0])):
            pos[count[0][i]] += count[1][i]

        # normalize the counts with total binding sites
        for cnt, i in enumerate(range(-flank, flank+1)):
            if i in erBS and erBS[i] != 0:
                pos[i] = pos[i]/erBS[i]
            else:
                pos[i] = 0
            if not "rand_{}".format(sam + 1) in finalDic.keys():
                finalDic["rand_{}".format(sam + 1)] = {i:pos[i]}
            else:
                finalDic["rand_{}".format(sam + 1)].update({i:pos[i]})

    # save the final output for each of the 1000 random sampling
    finalOutput=pd.DataFrame.from_dict(finalDic,orient='columns')

    # compute mean, confidence interval and standard deviation
    t = finalOutput.apply(mean_confidence_interval, axis=1)

    return(t)

# collect outputs
expected={}
sfile="%s/results/%s/%s/observedMutRate.csv" % (path, args.WGSused, args.ERBSused)

# Configure the logging
level = logging.DEBUG if args.debug else logging.INFO
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=level)
logging.debug(args)

logging.debug("Start")
result = extractDetails(sfile)
expected["mean"] = [i[0] for i in result]
expected["CIL"] = [i[1] for i in result]
expected["CIU"] = [i[2] for i in result]
expected["sd"] = [i[3] for i in result]

logging.debug("Done")

# save the final output in a dataframe to print
toprint = pd.DataFrame.from_dict(expected, orient='columns')
position = [ val for val in range(-flank, flank+1) ]
toprint.insert(0, 'position', position)
outfile = "%s/expectedMutRate_meta.csv" % (opath)
toprint.to_csv(outfile,header=True,sep="\t",index=False)
os.system("gzip %s" % outfile)