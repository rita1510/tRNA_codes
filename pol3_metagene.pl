#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use Array::Utils qw(:all);
my @trna;
my %trna_hash;
my %trna_type;
my $tcount=0;


open(trna,"<", "tRNA_ucsc_hg38.txt")or die("unable to open trna");

<trna>;
while(my $line=<trna>)
{
chomp($line);
my @array=split(/\t/,$line);
my @array1=split(/pos /,$array[9]);
my @array2=split(/<BR>/,$array1[1]);
my @array3=split(/-/,$array2[0]);

print $array3[0]." ".$array3[1]."\n";
	    push @trna, [$array[1], $array[2], $array[3], $array[6],$array[7],$array3[0],$array3[1]];
        
}
close(trna);

my @files=("Pol2_ctrl_hg38_F4_s_rmdup","Pol2_chip_hg38_F4_s_rmdup", "Pol3_Oler_hek293_hg38_F4_s_rmdup","SRR1230840_1_trimmedx6_aligned_all_hg38_F4_s_rmdup");
for (my $j=0; $j<=$#files; $j++)
{
my $file=$files[$j];




open(Polvals, "<", $file.".bedgraph") or die ("unable to open Polvals");
my $q=0;
print $q."\n";
my %polvals;
my $norm=qx(wc -l $file.bedgraph);
while (my $line=<Polvals>)
{
    $q+=(1/$norm);
    print $q." Pol3\n";
    chomp($line);
        my @array=split(/\t/,$line);


    for (my $v=$array[1]; $v<$array[2]; $v++)
    {
    $polvals{$array[0]}[$v]=$array[3];
    }


}
close(Polvals);

my @trna_pol3;
open(Polmetagene, ">", "metagene_$file.txt") or die ("unable to open Polmetagene");

for (my $u=0; $u<=$#trna; $u++)
{
    my $length=0;

    if ($trna[$u][3]eq"+")
    {
    my $start=$trna[$u][1];
    my $stop=$trna[$u][2];

    my $a=$start-200;
    my $b=$stop+200;
      for(my $v=$a; $v<=$b; $v++)
      {
        if($trna[$u][5]>0)
        {
        if($v<$start+$trna[$u][5]||$v>$start+$trna[$u][6])
        {
        my $val=$polvals{$trna[$u][0]}[$v]+0;
        print Polmetagene $val."\t";
        $length++;
        }
        }
        else
        {
        my $val=$polvals{$trna[$u][0]}[$v]+0;
        print Polmetagene $val."\t";
        $length++;
        }
      }
      
    }
    elsif ($trna[$u][3]eq"-")
    {
        my $start=$trna[$u][2];
    my $stop=$trna[$u][1];

    my $b=$start+200;
    my $a=$stop-200;

    
      for(my $v=$b; $v>=$a; $v--)
      {
        if($trna[$u][5]>0)
        {
        if($v<$stop+$trna[$u][5]||$v>$stop+$trna[$u][6])
        {
        my $val=$polvals{$trna[$u][0]}[$v]+0;
        print Polmetagene $val."\t";
                $length++;

        }
        }
        else
        {
        my $val=$polvals{$trna[$u][0]}[$v]+0;
        print Polmetagene $val."\t";
                $length++;

        }
      }
    }
          print $length."\n";

print Polmetagene "\n";
}
close(Polmetagene);
}
