#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);


open(ago1,"<", "AGO1_seqs.fasta")or die("unable to open ago1");

my %agos;
while(my $line=<ago1>)
{
chomp($line);
    if ($line=~">")
    {
    my @array=split(/ /, $line);
    my $line1=<ago1>;
    chomp($line1);
    push @{$agos{$array[1]."_".$array[4]}},[$line1,"AGO1"];
    }
}
close(ago1);
open(ago2,"<", "AGO2_seqs.fasta")or die("unable to open ago2");

while(my $line=<ago2>)
{
chomp($line);
    if ($line=~">")
    {
    my @array=split(/ /, $line);
    my $line1=<ago2>;
    chomp($line1);
    push @{$agos{$array[1]."_".$array[4]}},[$line1,"AGO2"];
    }
}
close(ago2);
open(ago3,"<", "AGO3_seqs.fasta")or die("unable to open ago3");

while(my $line=<ago3>)
{
chomp($line);
    if ($line=~">")
    {
    my @array=split(/ /, $line);
    my $line1=<ago3>;
    chomp($line1);
    push @{$agos{$array[1]."_".$array[4]}},[$line1,"AGO3"];
    }
}
close(ago3);
open(ago4,"<", "AGO4_seqs.fasta")or die("unable to open ago4");

while(my $line=<ago4>)
{
chomp($line);
    if ($line=~">")
    {
    my @array=split(/ /, $line);
    my $line1=<ago4>;
    chomp($line1);
    push @{$agos{$array[1]."_".$array[4]}},[$line1,"AGO4"];
    }
}
close(ago4);

open(dicer1,"<", "Dicer1_seqs.fasta")or die("unable to open dicer1");
my %dicers;
while(my $line=<dicer1>)
{
chomp($line);
    if ($line=~">")
    {
    my @array=split(/ /, $line);
    my $line1=<dicer1>;
    chomp($line1);
    push @{$dicers{$array[1]."_".$array[4]}},[$line1,"Dicer1"];

    }
}
close(dicer1);

open(dicer2,"<", "Dicer2_seqs.fasta")or die("unable to open dicer2");
while(my $line=<dicer2>)
{
chomp($line);
    if ($line=~">")
    {
    my @array=split(/ /, $line);
    my $line1=<dicer2>;
    chomp($line1);
    push @{$dicers{$array[1]."_".$array[4]}},[$line1,"Dicer2"];
   # print $array[1]."_".$array[4]." ".$line1."\n";
   }
}
close(dicer2);
open(dicer3,"<", "Dicer3_seqs.fasta")or die("unable to open dicer3");
while(my $line=<dicer3>)
{
chomp($line);
    if ($line=~">")
    {
    my @array=split(/ /, $line);
    my $line1=<dicer3>;
    chomp($line1);
    push @{$dicers{$array[1]."_".$array[4]}},[$line1,"Dicer3"];
      #  print $array[1]."_".$array[4]." ".$line1."\n";

    }
}
close(dicer3);

my $key;
my $a=0;
open(seqs,">", "Dicer_Ago.fasta")or die("unable to open seqs");

foreach $key (keys %dicers)
{
                     #   print $key." ".$#{$dicers{$key}}."\n";

        if ($#{$dicers{$key}}>0)
        {
                if ($#{$agos{$key}}>-1)
                {
                $a++;

                print seqs ">".$key."\n";
                print seqs $agos{$key}[0][0]."\n";
                }
        
        }
        
}

print $a."\n";