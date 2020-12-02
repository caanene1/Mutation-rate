#!/bin/bash

# script to compute the probability of occurrence of the 96 tri-nucleotide changes in each cancer type at the selected set of ERBS 
# written by Jiekun (Jackie) Yang

# read config file
sourcedir=$(dirname $BASH_SOURCE) # find the path where scripts are available
source $sourcedir/config.file # read config file to get data path variables

# usage
show_help() {
cat << EOF
Usage: ${0##*/} -w <WGSused> -e <ERBSused>
		-h  display this help and exit
    	-w  WGS data used (e.g., brcaeuc2t)
		-e  ERBS used (e.g., ERBS1Q)
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

# check input
if [[ -z $WGSused || -z $ERBSused ]];then 
    echo -e "\n\n\t***input missing check the usage below**\n\n"; show_help; exit 1;
fi

# result files
echo $WGSused
echo $ERBSused
respath="$path/results/$WGSused/$ERBSused"
zcat $respath/observedMutRate.bed.gz | cut -f 1,2,3,5,6 | sort -u -S 80% >$respath/observedMutRate.bed
bedtools intersect -wao -a $respath/observedMutRate.bed -b $mutpath/${WGSused}.txt.gz -f 1.0 | sort -k1,1 -k2,2n >$respath/observedMut.bed
for i in `seq 1 64`; 
do
	awk '{if($4=="'$i'"){print;}}' $respath/observedMutRate.bed | wc -l | perl -ne 'print "$_" x3'; 
done >$respath/ref_num.txt
for i in `seq 1 64 `; 
do
  case $i in
    1|2|3|4|17|18|19|20|33|34|35|36|49|50|51|52 )
        nucs=(T G C)
        ;;
    5|6|7|8|21|22|23|24|37|38|39|40|53|54|55|56 )
        nucs=(A G C)
        ;;
    9|10|11|12|25|26|27|28|41|42|43|44|57|58|59|60 )
        nucs=(A T C)
        ;;
    13|14|15|16|29|30|31|32|45|46|47|48|61|62|63|64 )
        nucs=(A T G)
        ;;
  esac
  for nuc in "${nucs[@]}"; 
  do 
	awk '{if($4=="'$i'" && $10=="'$nuc'"){print;}}' $respath/observedMut.bed | wc -l;
  done; 
done >>$respath/alt_num.txt

paste -d "\t" $respath/ref_num.txt $respath/alt_num.txt | awk 'BEGIN{FS=OFS="\t"}{print $1, $2, $2/$1}' - > $respath/mutSig.tsv
sed -i "1 i\Ref_${WGSused}_${ERBSused}\tAlt_${WGSused}_${ERBSused}\tProb_${WGSused}_${ERBSused}" $respath/mutSig.tsv
paste -d "\t" $mutpath/signature_probabilities_template.tsv $respath/mutSig.tsv > $mutpath/temp.tsv
mv $mutpath/temp.tsv $mutpath/signature_probabilities_template.tsv

rm -f $mutpath/temp.tsv
rm -f $respath/observedMut*.bed
