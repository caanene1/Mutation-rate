# script to compute mutation rate enrichment
# written by Sabarinathan Radhakrishnan (Sabarinathan R et al. Nature 2016, 532:264-267)
# modified by Jackie Yang

# import modules
import os
import argparse
import pandas as pd
import numpy as np
import scipy as sp
import scipy.stats
from math import sqrt
from math import log10
from math import log2
from statsmodels.sandbox.stats.multicomp import multipletests as mlpt

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='projectPath', help='same as path in config.file')
parser.add_argument('-w', dest='WGSused', help='WGS data used (e.g brcaeu)')
parser.add_argument('-e', dest='ERBSused', help='ERBS data set used (e.g allERBS)')
parser.add_argument('-i', dest='avgpeaksize', help='avergae peak size in the ERBS data set used (e.g allERBS)')
args = parser.parse_args()

if args.WGSused == None or args.ERBSused == None:
    parser.print_usage()
    sys.exit(1)

# set output path
opath = "%s/results/%s/%s/" % (args.projectPath, args.WGSused, args.ERBSused) # output path

# functions for enrichment analysis
def enrichmentTest(imutCounts, omutCounts, intsCounts, ontsCounts):
    
    # chi-square test
    chi2, p, dof, ex = sp.stats.chi2_contingency([[imutCounts, omutCounts], [(intsCounts), (ontsCounts)]], correction=False)
    # fold change
    foldChange = imutCounts/ex[0][0]
    
    return foldChange, p

def fishersTest(imutCounts, omutCounts, intsCounts, ontsCounts):
    # fisher's test
    oddsration, fisherPvalue = sp.stats.fisher_exact([[imutCounts, omutCounts], [intsCounts, ontsCounts]], alternative="greater")
    
    return fisherPvalue

# ERBS list with motif counts and size

#inside=200
inside=int(args.avgpeaksize)

flag=0
count = 0
output = pd.DataFrame(columns=[ 'insideLen', 'insideRatio', 'Flank1000_Ratio', 'Flank1000_FC','Flank1000_chiPvalue' ])
filename = opath+"observedMutRate.csv"

data = pd.read_csv(filename, header=0, sep="\t", index_col=None )

insideM = data['cnt'][ (data['position'] >= -inside) & (data['position'] <= inside) ].sum()
insideN = data['bp'][ (data['position'] >= -inside) & (data['position'] <= inside) ].sum()
arr = []
for flank in [1000]:

    outsideM = data['cnt'][ (data['position'] >= -(inside + flank)) & (data['position'] < -inside) ].sum() 
    outsideM += data['cnt'][ (data['position'] > inside) & (data['position'] <= (inside+flank)) ].sum()
    outsideN = data['bp'][ (data['position'] >= -(inside + flank)) & (data['position'] < -inside) ].sum() 
    outsideN += data['bp'][ (data['position'] > inside) & (data['position'] <= (inside+flank)) ].sum()
    
    ous = [outsideM, outsideN]
    if insideM > 0 or outsideM > 0:
        fold, p = enrichmentTest(insideM, outsideM, insideN, outsideN)
    else:
        fold = 0; p = 1
    arr += [ str(ous), fold, p ]
    
ins = [insideM, insideN]    

output.loc[0] = [ inside, str(ins) ] +  arr
    
output['Flank1000_chiPvalueAdj']= mlpt(output['Flank1000_chiPvalue'], method='fdr_bh')[1]
output['Flank1000_FC_log2'] = output['Flank1000_FC'].apply(lambda x: log2(x) if x > 0 else log2(1))
minPvalue = output['Flank1000_chiPvalueAdj'][output['Flank1000_chiPvalueAdj']!=0].min()
output['Flank1000_chiPvalueAdj_neglog10'] = output['Flank1000_chiPvalueAdj'].apply(lambda x: log10(x)*-1 if x > 0 else log10(minPvalue)*-1)

output.to_csv(opath+"mutation_enrichment.csv", header=True, sep="\t", index=False)
