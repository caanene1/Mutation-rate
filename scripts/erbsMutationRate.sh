#!/bin/bash

# script to compute mutation rate at different sets of ERBS
# written by Sabarinathan Radhakrishnan (Sabarinathan R et al. Nature 2016, 532:264-267)
# modified by Jiekun (Jackie) Yang

# read configuration file
sourcedir=$(dirname $BASH_SOURCE) # find the path where scripts are available
source $sourcedir/config.file # read config file to get data path variables

show_help() {
cat << EOF
Usage: ${0##*/} -w <WGSused> -e <ERBSused>
        -h  display this help and exit
        -w  WGS data used (e.g., brcaeu)
        -e  ERBS used (e.g., allERBS)
Example: ${0##*/} -w brcaeu -e allERBS
EOF
}

if [[ $# == 0 ]];then show_help;exit 1;fi

while getopts "w:e:h" opt; do
    case "$opt" in
        h) show_help;exit 0;;
        w) WGSused=$OPTARG;;
        e) ERBSused=$OPTARG;;
       '?')show_help >&2 exit 1 ;;
    esac
done

if [[ -z $WGSused || -z $ERBSused ]];then 
   echo -e "\n\n\t***input missing check the usage below**\n\n"; show_help; exit 1;
fi

# make a tmp directory
if [ -d $path ]; then 
   tmpdir=$(mktemp -d -p $path);
 else
   tmpdir=$(mktemp -d);
fi

# path for different input data
ERBSfile="$erbsPath/${ERBSused}.bed";  

# check if the TFBS file is accessible
if [ ! -f $ERBSfile ];then echo "input ERBS file $ERBSfile not available";exit 1;fi

# create output directory
if [ ! -d "$path/results" ];then mkdir "$path/results";fi
if [ ! -d "$path/results/$WGSused" ];then mkdir "$path/results/$WGSused";fi
if [ ! -d "$path/results/$WGSused/$ERBSused" ];then mkdir "$path/results/$WGSused/$ERBSused";fi

outfile="observedMutRate"

cut -f 1-3 $ERBSfile | sort -u | bedtools sort >$tmpdir/${outfile}_pre.bed

if [ `cat $tmpdir/${outfile}_pre.bed | wc -l` -eq 0 ];then exit;fi

# Extend ERBS from its mid-point 
# ==============================
# take +/- 1000 nucleotides from ERBS mid-point; since the intervals are zero based added one to the start position before computing mid-point
cat $tmpdir/${outfile}_pre.bed | perl -sane '$F[1]++;$midpoint=int(($F[2]+$F[1])/2);print $F[0],"\t",$midpoint,"\t",$midpoint,"\t",$F[1],"\t",$F[2],"\n";' | bedtools slop -i stdin -g $hg19genome -b $flank >$tmpdir/${outfile}_flank.bed

# Get the tri-nucleotide info.
# ============================
# first extract the sequence
# Added the "chr" to solve the error issue, also removed the "lank.bed | sed 's/chr//g' | bedtools getfasta" and added 
perl -sane '$F[0]=~s/chr//g;$F[1]-=2;$F[2]+=1;print "chr",$F[0],"\t",$F[1],"\t",$F[2],"\t",join("|",@F),"\n";' $tmpdir/${outfile}_flank.bed | bedtools getfasta -fi $hg19Fasta -bed stdin -fo $tmpdir/seq.txt -name
# added a-zA
sed -i 's/::[a-zA-Z0-9]*:[0-9]*-[0-9]*$//g' $tmpdir/seq.txt

# convert trinucleotides into numeric ids and remove any overlapping motifs in the flank
perl $path/scripts/makeTriNucProfile.pl $tmpdir/seq.txt $flank $allTFBSmerged | sort -k1,1 -k2,2n >$tmpdir/tmp.bed


# Apply filters
# =============
# filter if any positions that overlap with CDS, blacklisted region and low mappability
# changed below
bedtools intersect -wa -a $tmpdir/tmp.bed -b $blackListV1 -f 1.0 -v | bedtools intersect -wa -a stdin -b $cds -f 1.0 -v | bedtools intersect -wa -a stdin -b $mappable36mer -f 1.0 -v >$tmpdir/tmp_core.bed

# Compute Mutation rate 
# =====================
# map and count the number of mutations that overlap each position
bedtools intersect -c -a $tmpdir/tmp_core.bed -b $mutpath/$WGSused.txt.gz -f 1.0 | sort -k1,1 -k2,2n >$tmpdir/${outfile}.bed
# compute per position mutation rate
echo -e "position\tbp\tcnt\tAvg" >$tmpdir/${outfile}.csv
cat $tmpdir/${outfile}.bed | perl -sane '$bp{$F[3]}++;$cnt{$F[3]}+=$F[5];END{foreach $d (keys(%bp)){print $d,"\t",$bp{$d},"\t",$cnt{$d},"\t",($cnt{$d}/$bp{$d}),"\n";}}' | sort -n -k1,1 | grep -v pos >>$tmpdir/${outfile}.csv

# move the output files
gzip -9 $tmpdir/${outfile}.bed
cp $tmpdir/${outfile}.bed.gz $tmpdir/${outfile}.csv $path/results/$WGSused/$ERBSused/.

#remove the tmp directory
#rm -rf $tmpdir