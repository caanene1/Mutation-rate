# script to identify sets of regions with similar DHS hypersensitivity
# this script is downloaded from Perera D et al. Nature 2016, 532:259-263

"""Echo command line arguments


"""                    
                    
import sys                    
import os                    
import getopt                    
from io import StringIO                    
import random
import subprocess
from subprocess import *                    
                    
class run_code:                    
                                    
    def __init__(self):                
        self.h_hi = None
        self.h_lo = 0
        self.selectNum = 0    
        self.inputFile1 = None    
        self.inputFile2 = None    
        self.outFile1 = None
        self.outFile2 = None
        self.values = []
        self.DHSc1 = 3
        self.DHSc2 = 3
                                
    def countMutations(self):
        # read list 2 into memory
        list2 = []
        f2 = open(self.inputFile2,'r')
        for line2 in f2:
            tmp = line2.split("\t")
            tmp[-1] = tmp[-1].strip()
            list2.append((line2,float(tmp[self.DHSc2])))
        f2.close()

        print("Beginning hypersensitivity matching\n")
        out1 = open(self.outFile1,'w')
        out2 = open(self.outFile2,'w')
        f1 = open(self.inputFile1,'r')
        for line in f1:
            tmp = line.split("\t")
            tmp[-1] = tmp[-1].strip()
            tmplist = []
            if float(tmp[self.DHSc1]) <= self.h_hi and float(tmp[self.DHSc1]) >= self.h_lo:
                for l2 in list2:           
                    if float(tmp[self.DHSc1]) >= l2[1]-5 and float(tmp[self.DHSc1]) <= l2[1]+5:
                        tmplist.append(l2[0])
                if len(tmplist) > 0:
                    out1.write(line)
                    out2.write(random.choice(tmplist))
        f1.close()
        out1.close()
        out2.close()
        print("Done\n")

                            
    def ParseCommandLine(self, Arguments):                
        try:            
            (Options, Args) = getopt.getopt(Arguments, "i:j:h:o:p:c:d:k:")
        except:            
            print("Unknown option entered")        
            print(UsageInfo)        
            sys.exit(1)        
        OptionsSeen = {}            
        for (Option, Value) in Options:            
            OptionsSeen[Option] = 1        
            if Option == "-i":    
                self.inputFile1 = Value                
            elif Option == "-j":    
                self.inputFile2 = Value    
            elif Option == "-h":        
                self.h_hi = float(Value)
            elif Option == "-k":
                self.h_lo = float(Value)
            elif Option == "-o":        
                self.outFile1 = Value            
            elif Option == "-p":        
                self.outFile2 = Value
            elif Option == "-c":
                self.DHSc1 = int(Value)
            elif Option == "-d":
                self.DHSc2 = int(Value)
            else:        
                print("** Unknown option:", Option, Value)    
        if not self.inputFile1:            
            print("\n* Error: -i missing")        
            print(UsageInfo)        
            sys.exit(1)        
                    
    "UsageInfo = """"                
    
    Parameters:                
    -i [FILENAME] input1 list
    -j [FILENAME] input2 list                
    -h [INT] maximum hypersensitivity
    -k [INT] minimum hypersensitivity
    -o [FILENAME] output1 file
    -p [FILENAME] output2 file
    -c [INT] DHS col file1
    -d [INT] DHS col file2

    Example:                
    hypersensitivityMatching.py -i1 dataset1.txt -j dataset2.txt -h 500 -o results1.txt -p results2.txt     
    """

# Runs the actual program one function at a time.
def Main(bootstrap = None):
    global MAX_RESULTS_FILES_TO_PARSE
    if not bootstrap:
        bootstrap = run_code()
        bootstrap.ParseCommandLine(sys.argv[1:])
        bootstrap.countMutations()


# This is optimisation stuff, this basically never needs to be edited.
if __name__ == "__main__":                
    try:                
        import psyco            
        psyco.full()            
    except:                
        print("\npsyco not found - running without optimization\n")            
    Main()                
