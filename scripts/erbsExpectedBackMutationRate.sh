#!/bin/bash

# script to compute the expected mutation rate for a set of ERBS
# written by Sabarinathan Radhakrishnan (Sabarinathan R et al. Nature 2016, 532:264-267)
# modified by Jiekun (Jackie) Yang

sourcedir=$(dirname $BASH_SOURCE) # find the path where scripts are available
source $sourcedir/config.file # read config file to get data path variables

show_help() {
cat << EOF
Usage: ${0##*/} -w <WGSused> -e <ERBSused> -n <number of cores>
        -h  display this help and exit
        -w  WGS data used (e.g., brcaeuc2t)
        -e  ERBS used (e.g., ERBS1Q)
        -n  number of cores (default:4)
Example: ${0##*/} -w brcaeu -e allERBS -n 6
EOF
}

if [[ $# == 0 ]];then show_help;exit 1;fi

background=0
while getopts "w:e:n:h" opt; do
    case "$opt" in
        h) show_help;exit 0;;
        w) WGSused=$OPTARG;;
        e) ERBSused=$OPTARG;;
        n) ncpus=$OPTARG;;
       '?')show_help >&2 exit 1 ;;
    esac
done

# check input variables
if [[ -z $WGSused || -z $ERBSused ]];then 
   echo -e "\n\n\t***input missing check the usage below**\n\n"; show_help; exit 1;
fi

# assign number of cores
if [[ -z $ncpus ]];then ncpus=4;fi

outfile="expectedMutRate"
mutprob="$mutpath/signature_probabilities_template.tsv"

# source path (ie., the observed mutation rate)
srcpath="$path/results/$WGSused/$ERBSused"

# check if the input source file exits
if [ ! -f "$srcpath/observedMutRate.bed.gz" ];then
   echo -e "input file not available here $srcpath/observedMutRate.bed.gz";exit 1;
fi
if [ ! -f $mutprob ];then
   echo -e "mutation probability file not available here $mutprob";exit 1;
fi

col="${WGSused}_${ERBSused}"


# make a tmp directory
if [ -d $path ]; then 
   tmpdir=$(mktemp -d -p $path)
 else
   tmpdir=$(mktemp -d)
fi

# output path, the path where to store the final output
outpath="$srcpath/expected"

# create output directory
if [ ! -d "$outpath" ];then mkdir "$outpath";fi

# Compute expected Mutation rate 
# ===============================
# copy the observed mutation rate file to tmp directory
cp $srcpath/observedMutRate.bed.gz $tmpdir/.
gzip -d $tmpdir/observedMutRate.bed.gz

# perform the background mutation rate computation, script from Jordi Deu-Pons
#python $path/scripts/samplingMutation_v01.py -c Prob_${col} -s $mutprob -f $tmpdir/ -o $tmpdir/ -p "*[0-9].bed" --cores=$ncpus
python $path/scripts/samplingMutation_v01.py -c Prob_${col} -s $mutprob -f $tmpdir/ -o $tmpdir/ -p "*.bed" --cores=$ncpus

# move the output files
gzip -9 $tmpdir/rand.bin
mv $tmpdir/rand.bin.gz $outpath/${outfile}_rand.bin.gz

# collect the results, script from Sabarinathan Radhakrishnan
export path
python $path/scripts/get_results_erbsExpectedMutationRate.py -w $WGSused -e $ERBSused -o $WGSused/$ERBSused/ --cores=$ncpus

# remove the tmp directory
rm -rf $tmpdir
