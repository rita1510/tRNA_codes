#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);


open(misnobed,"<", "hg38_mirna_sno_allfields.bed")or die("unable to open misnobed");
my %exons;

<misnobed>;

while(my $line=<misnobed>)
{
    chomp($line);
    my @array=split(/\t/, $line);
    $exons{">hg38_wgRna_".$array[4]}=[$array[2],$array[3],$array[9]];
print $array[4]."\n";
}
close (misnobed);

open(scores,"<", "Score_S1_S2_S3_D1_D2_D3_misno_0_Ago1_1_1_1")or die("unable to open scores");
my %scores;
<scores>;
while(my $line=<scores>)
{
        chomp($line);

    my @array=split(/\t/, $line);
        my @array1=split(/_/, $array[1]);
    $scores{$array[0]}=[$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],$array[7],">hg38_wgRna_".$array1[2]];
    #print $scores{$array[0]}[8]."\n";
}
close(scores);

my $key;
open(scores_mirna,">", "Score_S1_S2_S3_D1_D2_D3_0_Ago1_1_1_1_miRNA")or die("unable to open scores");
open(scores_sno,">", "Score_S1_S2_S3_D1_D2_D3_0_Ago1_1_1_1_sno")or die("unable to open scores");
print scores_mirna "Sequence\ttRNA\tS1\tS2\tS3\tD1\tD2\tD3\n";
print scores_sno "Sequence\ttRNA\tS1\tS2\tS3\tD1\tD2\tD3\n";

foreach $key (keys %scores)
{
    print $exons{$scores{$key}[8]}[2]."\n";
    if($exons{$scores{$key}[8]}[2]eq"miRNA")
    {
        print scores_mirna $scores{$key}[0]."\t".$scores{$key}[1]."\t".$scores{$key}[2]."\t".$scores{$key}[3]."\t".$scores{$key}[4]."\t".$scores{$key}[5]."\t".$scores{$key}[6]."\t".$scores{$key}[7]."\n";
    }
    else
    {
        print scores_sno $scores{$key}[0]."\t".$scores{$key}[1]."\t".$scores{$key}[2]."\t".$scores{$key}[3]."\t".$scores{$key}[4]."\t".$scores{$key}[5]."\t".$scores{$key}[6]."\t".$scores{$key}[7]."\n";
    }
}
close(scores_mirna);
close(scores_sno);