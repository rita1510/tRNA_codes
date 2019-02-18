#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);

open(file1,"<", "AGO_Dicer_trna_7nt_Dicer1.fasta")or die("unable to open file1");
my %dicer1;
my %rdicer1;
my %dicers;
my %rdicers;
my %distr1;
my %distr2;
my %distr3;
my %agodistr;
while(my $line=<file1>)
{
        chomp($line);

    my @array=split(/ /, $line);
    my $line1=<file1>;
    chomp($line1);
if($array[3]>=1 && $array[4]>=1)#400
{
    push @{$dicer1{$array[0]."_".$array[1]."_".$array[2]}}, $line1;
    push @{$rdicer1{$line1}}, $array[0]."_".$array[1]."_".$array[2];
    print $line1." ".$array[0]."_".$array[1]."_".$array[2]."\n";
    push @{$dicers{$array[0]."_".$array[1]."_".$array[2]}}, "dicer1";
    push @{$rdicers{$line1}}, "dicer1";
 @{$rdicers{$line1}}=uniq(@{$rdicers{$line1}});
}
            push @{$distr1{$line1}}, $array[4];
            push @{$agodistr{$line1}}, $array[3];
            

}
close (file1);

open(file2,"<", "AGO_Dicer_trna_7nt_Dicer2.fasta")or die("unable to open file2");
#my %dicer2;
while(my $line=<file2>)
{
        chomp($line);

    my @array=split(/ /, $line);
    my $line1=<file2>;
    chomp($line1);
if($array[3]>=1 && $array[4]>=1)#40
{
    push @{$dicer1{$array[0]."_".$array[1]."_".$array[2]}}, $line1;
        push @{$rdicer1{$line1}}, $array[0]."_".$array[1]."_".$array[2];

        push @{$dicers{$array[0]."_".$array[1]."_".$array[2]}}, "dicer2";
                push @{$rdicers{$line1}}, "dicer2";
                 @{$rdicers{$line1}}=uniq(@{$rdicers{$line1}});

}
            push @{$distr2{$line1}}, $array[4];
            push @{$agodistr{$line1}}, $array[3];

}
close (file2);

open(file3,"<", "AGO_Dicer_trna_7nt_Dicer3.fasta")or die("unable to open file3");
#my %dicer3;
while(my $line=<file3>)
{
        chomp($line);

    my @array=split(/ /, $line);
    my $line1=<file3>;
    chomp($line1);
if($array[3]>=1 && $array[4]>=1)#15
{
    push @{$dicer1{$array[0]."_".$array[1]."_".$array[2]}}, $line1;
        push @{$rdicer1{$line1}}, $array[0]."_".$array[1]."_".$array[2];
        push @{$dicers{$array[0]."_".$array[1]."_".$array[2]}}, "dicer3";
        push @{$rdicers{$line1}}, "dicer3";
         @{$rdicers{$line1}}=uniq(@{$rdicers{$line1}});

}
                    push @{$distr3{$line1}}, $array[4];
            push @{$agodistr{$line1}}, $array[3];

}
close (file3);
#open(file4,">", "AGO_Dicer_misno_allDicers_0_Ago80_800_100_45.fasta")or die("unable to open file4");
#open(file4,">", "AGO_Dicer_misno_allDicers_0_Ago25_323_41_19.fasta")or die("unable to open file4");
open(file4,">", "AGO_Dicer_trna_7nt_allDicers_0_Ago1_1_1_1.fasta")or die("unable to open file4");

my $key;
foreach $key (keys %dicers)
{
    if($#{$dicers{$key}}>0)
    {
    print file4 $key."\n";
    print file4 $dicer1{$key}[0]."\n";
    }
}


close(file4);



#open(file5,">", "rev_AGO_Dicer_misno_allDicers_restricted_0_Ago30.fasta")or die("unable to open file5");
#
#my $key;
#foreach $key (keys %rdicers)
#{
#    if($#{$rdicers{$key}}>0)
#    {
#    print file5 $rdicer1{$key}[0]."\n";
#    print file5 $key."\n";
#    }
#}
#close(file5);
#
#open(distr1,">", "Distr_Dicer1.txt")or die("unable to open distr1");
#open(distr2,">", "Distr_Dicer2.txt")or die("unable to open distr2");
#open(distr3,">", "Distr_Dicer3.txt")or die("unable to open distr3");
#
#open(agodist,">", "Distr_AGO.txt")or die("unable to open agodist");
#foreach $key (keys %distr1)
#{
#    print distr1 $distr1{$key}[0]."\n";
#}
#foreach $key (keys %distr2)
#{
#    print distr2 $distr2{$key}[0]."\n";
#}
#foreach $key (keys %distr3)
#{
#    print distr3 $distr3{$key}[0]."\n";
#}
#foreach $key (keys %agodistr)
#{
#    print agodist $agodistr{$key}[0]."\n";
#}
#close(distr1);
#close(distr2);
#close(distr3);