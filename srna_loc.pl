#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use List::Util qw( min max sum);
use Array::Utils qw(:all);
    use POSIX;

open(seqs,"<", "AGO_Dicer_trna_7nt_allDicers_0_Ago25_323_41_19_unique.fasta")or die("unable to open seqs");
my @seqs;
while(my $line=<seqs>)
{
    chomp($line);
my $line1=<seqs>;
chomp($line1);
push @seqs, $line1;

}


open(trnaseq,"<", "ucsc_trna_hg38_7nt.fa")or die("unable to open trnaseq");
my %trnas;
my $trna;
while(my $line=<trnaseq>)
{
chomp($line);
if($line=~">")
{
my @array=split(/ /, $line);
$trna=substr($array[0],1);
#print $trna."\n";
}
else
{
    $trnas{$trna}.=$line;
}
}


close(trnaseq);
open(tab,">", "Rreduce.tab")or die("unable to open seqs");

my $key;
my %pos;
my %posseq;
foreach $key (keys %trnas)
{
    #print length($trnas{$key})."\n";
    for(my $u=0; $u<=$#seqs;$u++)
    {
        if(index($trnas{$key},$seqs[$u])>=0)
        {
    print tab $key." ".$seqs[$u]." ".index($trnas{$key},$seqs[$u])." ".(length($seqs[$u])+index($trnas{$key},$seqs[$u])-1)." ".length($trnas{$key})."\n";
    push @{$pos{$key}},[index($trnas{$key},$seqs[$u]),(length($seqs[$u])+index($trnas{$key},$seqs[$u])-1),length($trnas{$key}),$seqs[$u]];
        @{$pos{$key}}=sort{$a->[0]<=>$b->[0]||$a->[1]<=>$b->[1]} @{$pos{$key}};
        }
    }
    if(index($trnas{$key},"CGAGTCTCGGTGGAACCTCCA")>=0)
        {
        print $key."\n";
         }
    
}
die;
my %posreduced;
my %allseqcomb;
my %group;
my %rgroup;

open(posreduced,">", "posreduced.tab")or die("unable to open posreduced");

foreach $key (keys %pos)
{
    my @combine;
    for(my $u=0; $u<=$#{$pos{$key}};$u++)
    {
     for(my $v=$u; $v<=$#{$pos{$key}};$v++)
        {
            if($pos{$key}[$u][1]>=$pos{$key}[$v][0]+10&&$pos{$key}[$u][0]+10<=$pos{$key}[$v][1])
            {
            push @{$combine[$u]},$v;
         #   print $pos{$key}[$u][0]." ".$pos{$key}[$u][1]." ".$pos{$key}[$v][0]." ".$pos{$key}[$v][1]." ".$u." ".$v."\n";
            }
        }
    }
 my $u=0;
     my $i=0;

 while($u<=$#combine)
 {
    my @starts;
    my @stops;
    my @allseqs;
    for(my $v=0; $v<=$#{$combine[$u]}; $v++)
    {

      push @starts, $pos{$key}[$combine[$u][$v]][0];
      push @stops, $pos{$key}[$combine[$u][$v]][1];
      push @allseqs, $pos{$key}[$combine[$u][$v]][3]
        
     # print $combine[$u][$v]." ";
    }
    #print"\n";
    my $start=min(@starts);
    my $stop=max(@stops);
    
    $posreduced{$key}[$i]=[$start,$stop,($stop-$start+1)];
    $allseqcomb{$key}[$i]=@allseqs;
    for(my $g=0; $g<=$#allseqs; $g++)
    {
        push @{$group{$allseqs[$g]}},$key."_".$i;
        push  @{$rgroup{$key."_".$i}},$allseqs[$g];
        push @{$posseq{$allseqs[$g]}},$start."_".$stop."_".$pos{$key}[0][2];
        
    }
    $i++;
    $u=max(@{$combine[$u]})+1;
    #print $i." ".$u."\n";
 }
    for(my $v=0; $v<=$#{$posreduced{$key}}; $v++)
    {
        print posreduced $key."\t".$posreduced{$key}[$v][0]."\t".$posreduced{$key}[$v][1]."\t".$posreduced{$key}[$v][2]."\t".substr($trnas{$key},$posreduced{$key}[$v][0],$posreduced{$key}[$v][2])."\n";
    }
}
close(posreduced);

my %groupinfo;
my %posinfo;

my @GR;
foreach $key (keys %rgroup)
{
    

for(my $u=0; $u<=$#seqs;$u++)
{
    for(my $v=0; $v<=$#{$group{$seqs[$u]}};$v++)
    {
        if($key eq  $group{$seqs[$u]}[$v])
        {
            push @{$groupinfo{$key}},@{$group{$seqs[$u]}};
            @{$groupinfo{$key}}=uniq(@{$groupinfo{$key}});
                        push @{$posinfo{$key}},@{$posseq{$seqs[$u]}};
            @{$posinfo{$key}}=uniq(@{$posinfo{$key}});
                      #  print $key." ".$#{$groupinfo{$key}}."\n";


        }
 }
}

push @GR, join(" ",@{$groupinfo{$key}});
#print join(" ",@{$groupinfo{$key}})."\n";
}

@GR = uniq(@GR);
my @grseq;
for(my $u=0; $u<=$#GR;$u++)
{
my @array=split(" ",$GR[$u]);
    for(my $v=0; $v<=$#seqs; $v++)
    {
       # print $#array."\n";
        for(my $w=0; $w<=$#array; $w++)
        {
            for(my $i=0; $i<=$#{$group{$seqs[$v]}};$i++)
            {
            if($group{$seqs[$v]}[$i] eq $array[$w])
            {
            push @{$grseq[$u][1]},$seqs[$v];
            push @{$grseq[$u][2]},$v;
            for(my $k=0; $k<=$#{$posseq{$seqs[$v]}};$k++)
            {
            push @{$grseq[$u][3]},$posseq{$seqs[$v]}[$k];
            }
            @{$grseq[$u][3]}=uniq(@{$grseq[$u][3]});
                        @{$grseq[$u][1]}=uniq(@{$grseq[$u][1]});

            }
            }
        }

}
}

open(out,">", "sRNA_positions_out.tab")or die("unable to open posreduced");
open(stats,">", "sRNA_positions_stats.tab")or die("unable to open posreduced");
open(out2,">", "sRNA_positions_stats_start_stop.tab")or die("unable to open posreduced");
my @positions;
my @positions_ungr;

for(my $u=0; $u<=$#grseq;$u++)
{
    print stats $u."\t";
    for(my $v=0; $v<=$#{$grseq[$u][1]};$v++)
    {
        print out $grseq[$u][1][$v]." ";
         print stats $grseq[$u][1][$v]." ";
         print out2 $grseq[$u][1][$v]." ";

    }
    print out "\t";
    my @stats;
    my @stats2;

        for(my $v=0; $v<=$#{$grseq[$u][3]};$v++)
    {
       # print out $grseq[$u][3][$v]." ";
        my @array=split("_",$grseq[$u][3][$v]);
        print $array[1]." ".$array[2]."\n";
        my $percentage = ($array[1]-7)/($array[2]-14)*100;
                my $percentage2 = ($array[0]-7)/($array[2]-14)*100;

        print out "start: ".($array[0]-7)." stop: ".($array[1]-7)." tRNA length: ".($array[2]-14)." relative postion: ".int($percentage)." %\t";

        push @stats, $percentage;
        push @stats2, $percentage2;
        push @positions_ungr,[$percentage2,$percentage];
    }
    print out "\n";

    print stats "\t".@stats."\t".min(@stats)."\t".max(@stats)."\t".mean(@stats)."\n";
        print out2 "\t".@stats."\t".mean(@stats2)."\t".mean(@stats)."\n";
    push @positions,[mean(@stats2),mean(@stats)];
}
open(vect,">", "sRNA_positions_vector.tab")or die("unable to open vect");
open(vect_ungr,">", "sRNA_positions_vector_ungrouped.tab")or die("unable to open vect");

my @vector;
my @min;
for(my $u=0; $u<=$#positions; $u++)
{
push @min, $positions[$u][0];

}
@min=sort{$a <=> $b}@min;
print $min[0]."\n";



for(my $u=0; $u<=$#positions; $u++)
{
for(my $v=ceil($positions[$u][0])-ceil($min[0]); $v<=ceil($positions[$u][1])-ceil($min[0]); $v++)
{
        $vector[$v][0]=$v+ceil($min[0]);
    $vector[$v][1]=0;
}
}

for(my $u=0; $u<=$#positions; $u++)
{
for(my $v=ceil($positions[$u][0])-ceil($min[0]); $v<=ceil($positions[$u][1])-ceil($min[0]); $v++)
{
    $vector[$v][1]++;
}
}
for(my $u=0; $u<=$#vector; $u++)
{
print vect $vector[$u][0]."\t".$vector[$u][1]."\n";
}
my @vector_ungr;
my @min_ungr;
for(my $u=0; $u<=$#positions_ungr; $u++)
{
push @min_ungr, $positions_ungr[$u][0];

}
@min_ungr=sort{$a <=> $b}@min_ungr;
print $min_ungr[0]."\n";



for(my $u=0; $u<=$#positions_ungr; $u++)
{
for(my $v=ceil($positions_ungr[$u][0])-ceil($min_ungr[0]); $v<=ceil($positions_ungr[$u][1])-ceil($min_ungr[0]); $v++)
{
        $vector_ungr[$v][0]=$v+ceil($min_ungr[0]);
    $vector_ungr[$v][1]=0;
}
}

for(my $u=0; $u<=$#positions_ungr; $u++)
{
for(my $v=ceil($positions_ungr[$u][0])-ceil($min_ungr[0]); $v<=ceil($positions_ungr[$u][1])-ceil($min_ungr[0]); $v++)
{
    $vector_ungr[$v][1]++;
}
}
for(my $u=0; $u<=$#vector_ungr; $u++)
{
print vect_ungr $vector_ungr[$u][0]."\t".$vector_ungr[$u][1]."\n";
}

sub mean {
    return sum(@_)/@_;
}