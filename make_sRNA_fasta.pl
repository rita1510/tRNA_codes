#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);


open(sam,"<", "Rybak_Wolf_25416952_Dicer_PAR_CLIP_3.bowtie.sort.reheader_chroms.sam")or die("unable to open sam");
my %sequence;
while(my $line=<sam>)
{
    chomp($line);
 my @array=split(/\t/, $line);
$sequence{$array[0]}=$array[9];
}
close(sam);

open(sRNA,"<", "Dicer_3_exons_intersect_names_restricted")or die("unable to open sRNA");
my %names;
my @coords;
my @name;
my %coord;
while(my $line=<sRNA>)
{
    chomp($line);
 my @array=split(/\t/, $line);
@{$names{$array[5]."_".$array[6]."_".$array[7]}}=$array[8];
push @coords, $array[5]."\t".$array[6]."\t".$array[7]."\tn\t255\t".$array[4];
push @name, $array[8];
$coord{$array[8]}=$array[5]."_".$array[6]."_".$array[7];
}
close(sRNA);

@coords=uniq(@coords);
@coords=sort{$a cmp $b}(@coords);

@name=uniq(@name);

open(fasta,">", "Dicer3_test.fasta")or die("unable to open sRNA");
open(cor,">", "Dicer3_intersect_coords.bed")or die("unable to open sRNA");

for(my $u=0; $u<=$#name; $u++)
{
    print fasta ">".$name[$u]." ".$coord{$name[$u]}."\n";
    print fasta $sequence{$name[$u]}."\n";
}
for(my $u=0; $u<=$#coords; $u++)
{
print cor $coords[$u]."\n";
}
close(cor);