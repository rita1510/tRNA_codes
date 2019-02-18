#!/usr/bin/perl

use strict;
use List::Util;
use List::MoreUtils qw(uniq);
use Array::Utils qw(:all);

my @trna;
my %trna_hash;
my %trna_type;
my $tcount=0;

my $stretchto=100;
open(trna,"<", "tRNA_ucsc_hg38.txt")or die("unable to open trna");

<trna>;
while(my $line=<trna>)
{
chomp($line);
my @array=split(/\t/,$line);

	    push @trna, [$array[1], $array[2], $array[3], $array[6],$array[7],$array[9]];
        my $size=$array[3]-$array[2];
        print $size." ".$array[9]."\n";
        
}
close(trna);
my @files=("Pol2_ctrl_hg38_F4_s_rmdup","Pol2_chip_hg38_F4_s_rmdup", "Pol3_Oler_hek293_hg38_F4_s_rmdup","SRR1230840_1_trimmedx6_aligned_all_hg38_F4_s_rmdup");
for (my $j=0; $j<=$#files; $j++)
{
my $file=$files[$j];


open(metagene, "<", "metagene_".$file.".txt") or die ("unable to open Polmetagene");
open(metagene_binned, ">", "metagene_stretched_".$stretchto."_".$file.".txt") or die ("unable to open Polmetagene");

while(my $line=<metagene>)
{
chomp($line);
my @array=split(/\t/,$line);
for(my $v=0; $v<=199;$v++)
{
    print metagene_binned $array[$v]."\t";
}
my $length=$#array+1-400;

my @uvalue;
for(my $v=200; $v<=200+$length;$v++)
{
    $uvalue[$v-200]=$array[$v];
}
my $binsize=$stretchto/$length;
print $binsize."\n";
my @value;
my $w=1;
    while($w<=$stretchto)
    {
        for(my $j=0; $j<$length; $j++)
        {
            if($w>=$j*$binsize&&$w<($j+1)*$binsize)
            {
               if($w+1>=($j+1)*$binsize)
               {
                $value[$w]=(($j+1)*$binsize-$w)*$uvalue[$j]+($w+1-($j+1)*$binsize)*$uvalue[$j+1]
               }
               elsif($w+1<($j+1)*$binsize)
               {
                $value[$w]=$uvalue[$j]
               }
            }
        }
        $w+=1;
    }
    $value[$stretchto]=$uvalue[$length];
for(my $v=1; $v<=$stretchto;$v++)
{
    print metagene_binned $value[$v]."\t";
}
for(my $v=200+$length+1; $v<=$#array;$v++)
{
    print metagene_binned $array[$v]."\t";
}
print metagene_binned "\n";
}
close(metagene);
close(metagene_binned);
}