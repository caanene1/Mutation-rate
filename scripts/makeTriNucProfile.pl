#!/usr/bin/perl -w

use strict;

# scritp to generate the tri-nucleotide profile for a given chromosomal range
# writted by Sabarinathan Radhakrishnan

# prepare the trinucleotide list
my @nuc=qw(A T G C);
my %triList;my $cnt=1;
foreach my $f (@nuc)
{
  foreach my $s (@nuc)
  {
   foreach my $t (@nuc)
   {
     my $tri=$f.$s.$t;
     #print $cnt,"\t",$tri,"\n";
     $triList{$tri}=$cnt++;
   }
  }
}

# read the all motif files
open(IN,"zcat $ARGV[2] |") or die "error in opening the TFBS co-ordinates"; 

my %TFregions;
while(<IN>)
{
 my @s=split(/\s+/, $_);
 $s[0]=~s/chr//g; 
 for(my $i=$s[1]+1;$i<=$s[2];$i++) # add plus one to the start position, if the motifs are zero based
 {
  $TFregions{"$s[0]-$i"}=0; # save the chromosmal position of TFbound regions
 }
}
close IN;

# input sequence file
open(IN,$ARGV[0]) or die "error in opening the input sequence file";
# flank size
my $flank=$ARGV[1]; #flank size from mid-point

while(<IN>)
{
 chomp; my $header;my $seq;
 if($_=~/^>/){$header=$_;$seq=<IN>;chomp $seq;$header=~s/>//g;}
 my @s=split(/\|/, $header);
 my $chr=$s[0];
 my $start=$s[1];my $end=$s[2];# start and end of the flanking regions
 my $mstart=$s[3];my $mend=$s[4];# start and end of the motif regions
 my %motifregions;
 for(my $i=$mstart;$i<=$mend;$i++)
 { 
   $motifregions{"$chr-$i"}++;
 }
 next, if ($seq !~/^[A-Z]/);
 my $pos;
 $pos=-$flank;

 my $j=0;
 for(my $i=$start+2;$i<=($end-1);$i++)
 {
   my $sub=substr($seq,$j,3);   
   if(exists $triList{$sub}){
   # check if the base is in the flank side don't overlap with any tfbound region
   if(!exists $motifregions{"$chr-$i"} && !exists $TFregions{"$chr-$i"}) 
   {print "chr$chr","\t",$i,"\t",$i,"\t",$pos,"\t",$triList{$sub},"\n";}
   # if the base is present in motif region, don't apply filter
   elsif(exists $motifregions{"$chr-$i"})
   {print "chr$chr","\t",$i,"\t",$i,"\t",$pos,"\t",$triList{$sub},"\n";}
   }      
   $j++;
   $pos++;
 } 
}
close IN;