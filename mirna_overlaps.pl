#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);

open(misnobed,"<", "hg38_mirna_sno_allfields.bed")or die("unable to open misnobed");
my %introns;

<misnobed>;

while(my $line=<misnobed>)
{
    my @array=split(/\t/, $line);
    $introns{"hg38_wgRna_".$array[4]}=[$array[2],$array[3]];

}
close (misnobed);

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
close(misnoseq);


open(sam_AGO1,"<", "AGO1_mappings_misno_sort.sam")or die("unable to open sam_AGO1");
my @AGO1;
my %AGO_hash;
while(my $line=<sam_AGO1>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {
        if($array[3]>=$introns{$array[2]}[0]&&$array[3]+$length<=$introns{$array[2]}[1])
        {
            #print $introns{$array[2]}[0]." ".$introns{$array[2]}[1]." ".$array[2]." ".$array[3]." INTRONMAPPER\n";
        }
        else
        {
            push @AGO1, \@array;
            push @{$AGO_hash{$array[2]}}, \@array;

        }
    }
}
close(sam_AGO1);
print $#AGO1."\n";

open(sam_AGO2,"<", "AGO2_mappings_misno_sort.sam")or die("unable to open sam_AGO2");
my @AGO2;
while(my $line=<sam_AGO2>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {
        if($array[3]>=$introns{$array[2]}[0]&&$array[3]+$length<=$introns{$array[2]}[1])
        {
            #print $introns{$array[2]}[0]." ".$introns{$array[2]}[1]." ".$array[2]." ".$array[3]." INTRONMAPPER\n";
        }
        else
        {
            push @AGO2, \@array;
            push @{$AGO_hash{$array[2]}}, \@array;

        }
    }
}
close(sam_AGO2);
print $#AGO2."\n";

open(sam_AGO3,"<", "AGO3_mappings_misno_sort.sam")or die("unable to open sam_AGO3");
my @AGO3;
while(my $line=<sam_AGO3>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {
        if($array[3]>=$introns{$array[2]}[0]&&$array[3]+$length<=$introns{$array[2]}[1])
        {
            #print $introns{$array[2]}[0]." ".$introns{$array[2]}[1]." ".$array[2]." ".$array[3]." INTRONMAPPER\n";
        }
        else
        {
            push @AGO3, \@array;
            push @{$AGO_hash{$array[2]}}, \@array;

        }
    }
}
close(sam_AGO3);
print $#AGO3."\n";

open(sam_AGO4,"<", "AGO4_mappings_misno_sort.sam")or die("unable to open sam_AGO4");
my @AGO4;
while(my $line=<sam_AGO4>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {
        if($array[3]>=$introns{$array[2]}[0]&&$array[3]+$length<=$introns{$array[2]}[1])
        {
            #print $introns{$array[2]}[0]." ".$introns{$array[2]}[1]." ".$array[2]." ".$array[3]." INTRONMAPPER\n";
        }
        else
        {
            push @AGO4, \@array;
            push @{$AGO_hash{$array[2]}}, \@array;
            
        }
    }
}
close(sam_AGO4);
print $#AGO4."\n";

open(sam_Dicer1,"<", "Dicer_rep1_mappings_misno_sort.sam")or die("unable to open sam_Dicer1");
my @Dicer1;
my %Dicer_hash;
my %dicers;
while(my $line=<sam_Dicer1>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {
        if($array[3]>=$introns{$array[2]}[0]&&$array[3]+$length<=$introns{$array[2]}[1])
        {
            #print $introns{$array[2]}[0]." ".$introns{$array[2]}[1]." ".$array[2]." ".$array[3]." INTRONMAPPER\n";
        }
        else
        {
            push @Dicer1, \@array;
            push @{$Dicer_hash{$array[2]}}, \@array;
           push @{$dicers{$array[2]."=".$array[3]."=".length($array[9])}},"dicer1";
            @{$dicers{$array[2]."=".$array[3]."=".length($array[9])}}=uniq(@{$dicers{$array[2]."=".$array[3]."=".length($array[9])}});      

        }
    }
}
close(sam_Dicer1);
print $#Dicer1."\n";

open(sam_Dicer2,"<", "Dicer_rep2_mappings_misno_sort.sam")or die("unable to open sam_Dicer2");
my @Dicer2;
while(my $line=<sam_Dicer2>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {
        if($array[3]>=$introns{$array[2]}[0]&&$array[3]+$length<=$introns{$array[2]}[1])
        {
            #print $introns{$array[2]}[0]." ".$introns{$array[2]}[1]." ".$array[2]." ".$array[3]." INTRONMAPPER\n";
        }
        else
        {
            push @Dicer2, \@array;
            push @{$Dicer_hash{$array[2]}}, \@array;
           push @{$dicers{$array[2]."=".$array[3]."=".length($array[9])}},"dicer2";
            @{$dicers{$array[2]."=".$array[3]."=".length($array[9])}}=uniq(@{$dicers{$array[2]."=".$array[3]."=".length($array[9])}});      

        }
    }
}
close(sam_Dicer2);
print $#Dicer2."\n";

open(sam_Dicer3,"<", "Dicer_rep3_mappings_misno_sort.sam")or die("unable to open sam_Dicer3");
my @Dicer3;
while(my $line=<sam_Dicer3>)
{
    my @array=split(/\t/, $line);
    my $length=length($array[9]);
    if($array[1]==0 && $length>=19 && $length<=22)
    {
        if($array[3]>=$introns{$array[2]}[0]&&$array[3]+$length<=$introns{$array[2]}[1])
        {
            #print $introns{$array[2]}[0]." ".$introns{$array[2]}[1]." ".$array[2]." ".$array[3]." INTRONMAPPER\n";
        }
        else
        {
            push @Dicer3, \@array;
            push @{$Dicer_hash{$array[2]}}, \@array;
           push @{$dicers{$array[2]."=".$array[3]."=".length($array[9])}},"dicer3";
            @{$dicers{$array[2]."=".$array[3]."=".length($array[9])}}=uniq(@{$dicers{$array[2]."=".$array[3]."=".length($array[9])}});      

        }
    }
}
close(sam_Dicer3);
print $#Dicer3."\n";
my $key;
my %ago_num;
foreach $key (keys %AGO_hash)
{
    for(my $u=0; $u<=$#{$AGO_hash{$key}}; $u++)
    {
        $ago_num{$key."=".$AGO_hash{$key}[$u][3]."=".length($AGO_hash{$key}[$u][9])}++;
    }
}
print "AGO_num\n";
my $key;
my %dicer_num;
foreach $key (keys %Dicer_hash)
{
    for(my $u=0; $u<=$#{$Dicer_hash{$key}}; $u++)
    {
        $dicer_num{$key."=".$Dicer_hash{$key}[$u][3]."=".length($Dicer_hash{$key}[$u][9])}++;
        print $#{$dicers{$key."=".$Dicer_hash{$key}[$u][3]."=".length($Dicer_hash{$key}[$u][9])}}."\n";

    }
}
print "Dicer_num\n";

my $key;
#my %coords_ov;
#foreach $key (keys %AGO_hash)
#{
#    for(my $u=0; $u<=$#{$AGO_hash{$key}}; $u++)
#    {
#            for(my $v=0; $v<=$#{$Dicer_hash{$key}}; $v++)
#            {
#                    if($AGO_hash{$key}[$u][3]==$Dicer_hash{$key}[$v][3]&&length($AGO_hash{$key}[$u][9])==length($Dicer_hash{$key}[$v][9]))
#                    {
#                        push @{$coords_ov{$key}}, $AGO_hash{$key}[$u][3]."_".length($AGO_hash{$key}[$u][9]);
#                    }
#            }
#    }
#    @{$coords_ov{$key}}=uniq(@{$coords_ov{$key}});
#}
open(file,">", "AGO_Dicer_misno.fasta")or die("unable to open file");

foreach $key (keys %ago_num)
{
    if(exists($dicer_num{$key}) && $#{$dicers{$key}}>0)
    {
    print $key."\n";
    my @array=split(/=/,$key);
    print file ">".$array[0]." ".$array[1]." ".$array[2]." ".$ago_num{$key}." ".$dicer_num{$key}."\n";
    my $pr=substr($seq{$array[0]},$array[1]-1,$array[2]);
    print file $pr."\n";
    print  ">".$array[0]." ".$array[1]." ".$array[2]." ".$ago_num{$key}." ".$dicer_num{$key}."\n";
    print  substr($seq{$array[0]},$array[1]-1,$array[2])."\n";
    
    }
}


