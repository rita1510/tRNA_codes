#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use Array::Utils qw(:all);


my @trna;
my %trna_hash;
my %trna_type;

open(trna,"<", "tRNA_ucsc_hg38.txt")or die("unable to open trna");

<trna>;
while(my $line=<trna>)
{
chomp($line);
my @array=split(/\t/,$line);

	    push @trna, [$array[1], $array[2], $array[3], $array[6],$array[7],$array[9]];
        push @{$trna_hash{$array[1]}},[$array[2], $array[3],$array[7]];
        my $size=$array[3]-$array[2];
     #   print $size." ".$array[7]."\n";
        
}
close(trna);

open(dicer, "<","MACS_Dicer_vs_Ctrl_q05_broad/NA_peaks.broadPeak") or die("unable to open dicer peaks");
my %dicer;
while (my $line=<dicer>)
{
chomp($line);
my @array=split(/\t/,$line);
push @{$dicer{$array[0]}}, [$array[1]+1,$array[2]];
}
close(dicer);

open(pol3, "<","MACS_Pol3_vs_Ctrl_q05_broad/NA_peaks.broadPeak") or die("unable to open pol3 peaks");
my %pol3;
while (my $line=<pol3>)
{
chomp($line);
my @array=split(/\t/,$line);
push @{$pol3{$array[0]}}, [$array[1]+1,$array[2]];
}
close(pol3);

my $key;

foreach $key (keys %trna_hash)
{
    for(my $u=0; $u<=$#{$trna_hash{$key}}; $u++)
    {
        for(my $v=0; $v<=$#{$dicer{$key}}; $v++)
        {
            if($trna_hash{$key}[$u][1]>=$dicer{$key}[$v][0]-100 &&$trna_hash{$key}[$u][0]<=$dicer{$key}[$v][1]+100)
            {
               $trna_hash{$key}[$u][3]="dicer" 
            }
        }
        for(my $v=0; $v<=$#{$pol3{$key}}; $v++)
        {
            if($trna_hash{$key}[$u][1]>=$pol3{$key}[$v][0] &&$trna_hash{$key}[$u][0]<=$pol3{$key}[$v][1])
            {
               $trna_hash{$key}[$u][4]="pol3" 
            }
        }
    }
}
    my @number_trnas;#[0]-> all trnas. [1]->with dicer, [2]->with pol3; [3]->with dicer and pol3

foreach $key (keys %trna_hash)
{

    for(my $u=0; $u<=$#{$trna_hash{$key}}; $u++)
    {
        $number_trnas[0]++;
        $trna_type{$trna_hash{$key}[$u][2]}[0]++;
        if($trna_hash{$key}[$u][3]eq"dicer")
        {
           $number_trnas[1]++;
            $trna_type{$trna_hash{$key}[$u][2]}[1]++;

        }
        if($trna_hash{$key}[$u][4]eq"pol3")
        {
           $number_trnas[2]++;
                       $trna_type{$trna_hash{$key}[$u][2]}[2]++;

        }
        if($trna_hash{$key}[$u][3]eq"dicer"&&$trna_hash{$key}[$u][4]eq"pol3")
        {
           $number_trnas[3]++;
                       $trna_type{$trna_hash{$key}[$u][2]}[3]++;

        }
        if($trna_hash{$key}[$u][3]eq"dicer"&&$trna_hash{$key}[$u][4]ne"pol3")
        {
           $number_trnas[4]++; 
        }
        if($trna_hash{$key}[$u][3]ne"dicer"&&$trna_hash{$key}[$u][4]ne"pol3")
        {
           $number_trnas[5]++; 
        }
        if($trna_hash{$key}[$u][3]ne"dicer"&&$trna_hash{$key}[$u][4]eq"pol3")
        {
           $number_trnas[6]++; 
        }
    }
}

print "q05 ".$number_trnas[0]." ".$number_trnas[1]." ".$number_trnas[2]." ".$number_trnas[3]." ".$number_trnas[4]." ".$number_trnas[5]." ".$number_trnas[6]."\n";

foreach $key (keys %trna_type)
{
    print $key."\t".$trna_type{$key}[0]."\t".$trna_type{$key}[1]."\t".$trna_type{$key}[2]."\t".$trna_type{$key}[3]."\n";
}
#subroutines=============================================================================
#subroutine to get the maximum of an array=======================================
sub maxval
{

my @array=@_;
my $max_value = $array[0];

for(my $i = 0; $i<=$#array;$i++)
{
  if($array[$i] > $max_value)
        {
        $max_value = $array[$i];
        }
}
return $max_value;
}
#=============maxindex==========================
sub maxindex
{

my @array=@_;
my $max_index;
my $max_value=$array[0];

for(my $i = 0; $i<=$#array;$i++)
{
  if($array[$i] > $max_value)
        {
	$max_value = $array[$i];
	$max_index = $i;
        }
}
return $max_index;
}

#==============minval==========================
sub minval
{

my @array=@_;
my $min_value = $array[0];

for(my $i = 0; $i<=$#array;$i++)
{
  if($array[$i] < $min_value)
        {
        $min_value = $array[$i];
        }
}
return $min_value;
}
#=============maxindex==========================
sub minindex
{

my @array=@_;
my $min_index;
my $min_value=$array[0];

for(my $i = 0; $i<=$#array;$i++)
{
  if($array[$i] < $min_value)
        {
	$min_value = $array[$i];
	$min_index = $i;
        }
}
return $min_index;
}