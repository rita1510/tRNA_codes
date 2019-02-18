#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);

my $S1_norm=5518809+9744029;
my $S2_norm=6534397+11500309;
my $S3_norm=6174733+10794058;
my $D1_norm=7067107+10605728;
my $D2_norm=5360819+7886349;
my $D3_norm=6013141+8418682; 



open(misnoseq,"<", "hg38_mirna_sno.fa")or die("unable to open misnoseq");
my %seq;
my $misno;
while(my $line=<misnoseq>)
{
chomp($line);
if($line=~">")
{
my @array=split(/ /, $line);
$misno=substr($array[0],1);
print $misno."\n";
}
else
{
    $seq{$misno}.=$line;
}
}
close (misnoseq);


open(file,"<", "AGO_Dicer_misno_allDicers_0_Ago1_1_1_1_unique.fasta")or die("unable to open file");
my %AD;
my %AD_seq;

while(my $line=<file>)
{
chomp($line);
if($line=~">")
{
$misno=$line;
}
else
{
    $AD{$misno}=$line;
    $AD_seq{$line}=$misno;

}
}
close(file);


open(sam_S1,"<", "S1_mappings_misno_sort.sam")or die("unable to open sam_S1");
my @S1;
my %S1_hash;
while(my $line=<sam_S1>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {

            push @S1, \@array;
            my $pr=substr($seq{$array[2]},$array[3]-1,length($array[9]));
            push @{$S1_hash{$pr}}, $array[2];

    }
}
close(sam_S1);
print $#S1."\n";

open(sam_S2,"<", "S2_mappings_misno_sort.sam")or die("unable to open sam_S2");
my @S2;
my %S2_hash;
while(my $line=<sam_S2>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {

            push @S2, \@array;
            my $pr=substr($seq{$array[2]},$array[3]-1,length($array[9]));
            push @{$S2_hash{$pr}}, $array[2];
        
    }
}
close(sam_S2);
print $#S2."\n";

open(sam_S3,"<", "S3_mappings_misno_sort.sam")or die("unable to open sam_S3");
my @S3;
my %S3_hash;
while(my $line=<sam_S3>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {

            push @S3, \@array;
            my $pr=substr($seq{$array[2]},$array[3]-1,length($array[9]));
            push @{$S3_hash{$pr}}, $array[2];
        
    }
}
close(sam_S3);
print $#S3."\n";

open(sam_D1,"<", "D1_mappings_misno_sort.sam")or die("unable to open sam_D1");
my @D1;
my %D1_hash;
while(my $line=<sam_D1>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {

            push @D1, \@array;
            my $pr=substr($seq{$array[2]},$array[3]-1,length($array[9]));
            push @{$D1_hash{$pr}}, $array[2];
        
    }
}
close(sam_D1);
print $#D1."\n";

open(sam_D2,"<", "D2_mappings_misno_sort.sam")or die("unable to open sam_D2");
my @D2;
my %D2_hash;
while(my $line=<sam_D2>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {

            push @D2, \@array;
            my $pr=substr($seq{$array[2]},$array[3]-1,length($array[9]));
            push @{$D2_hash{$pr}}, $array[2];
        
    }
}
close(sam_D2);
print $#D2."\n";

open(sam_D3,"<", "D3_mappings_misno_sort.sam")or die("unable to open sam_D3");
my @D3;
my %D3_hash;
while(my $line=<sam_D3>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {

            push @D3, \@array;
            my $pr=substr($seq{$array[2]},$array[3]-1,length($array[9]));
            push @{$D3_hash{$pr}}, $array[2];
        
    }
}
close(sam_D3);
print $#D3."\n";

open(scores,">", "Score_S1_S2_S3_D1_D2_D3_misno_0_Ago1_1_1_1")or die("unable to open scores");

print scores "Sequence\ttRNA\tS1\tS2\tS3\tD1\tD2\tD3\n";

my $key;
foreach $key(keys %AD_seq)
{
    my $S1_score=($#{$S1_hash{$key}}+1)/$S1_norm;
    my $S2_score=($#{$S2_hash{$key}}+1)/$S2_norm;
    my $S3_score=($#{$S3_hash{$key}}+1)/$S3_norm;
    my $D1_score=($#{$D1_hash{$key}}+1)/$D1_norm;
    my $D2_score=($#{$D2_hash{$key}}+1)/$D2_norm;
    my $D3_score=($#{$D3_hash{$key}}+1)/$D3_norm;
    if($S1_score>0 || $S2_score>0 || $S3_score>0)
    {
    print scores $key."\t".$AD_seq{$key}."\t".$S1_score."\t".$S2_score."\t".$S3_score."\t".$D1_score."\t".$D2_score."\t".$D3_score."\n";
    }
}