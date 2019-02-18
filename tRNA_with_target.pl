#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);


open(miranda,"<", "miranda-AGO_Dicer_trna_7nt_allDicers-HEK_dicer_results4_dicer005_ago005_drosha05_all.txt")or die("unable to open miranda");
my %miranda;
while(my $line=<miranda>)
{
    my @array=split(/'/, $line);

        if($line=~"Query")
        {
          #  print $array[1]."\n";
            my $string = uc reverse $array[1];
            $string =~ tr/5\-\ //d;
        $miranda{$string}="miranda-target";
        }
}
close (miranda);

open(scores,"<", "../trna_mappings/Score_S1_S2_S3_D1_D2_D3_AD_7nt_allDicers_0_Ago25_323_41_19")or die("unable to open scores");
my %scores;
<scores>;
while(my $line=<scores>)
{
        chomp($line);

    my @array=split(/\t/, $line);
        my @array1=split(/ /, $array[1]);

    $scores{$array[0]}=[$array[0],$array[1],$array[2],$array[3],$array[4],$array[5],$array[6],$array[7],$array1[0]];
}
close(scores);

my $key;
open(scores_target,">", "Score_S1_S2_S3_D1_D2_D3_AD_7nt_allDicers_0_Ago25_323_41_19_miranda-target")or die("unable to open scores_target");
print scores_target "Sequence\ttRNA\tS1\tS2\tS3\tD1\tD2\tD3\n";

foreach $key (keys %miranda)
{
    if (exists $scores{$key})
    {
        
    
    print scores_target $scores{$key}[0]."\t".$scores{$key}[1]."\t".$scores{$key}[2]."\t".$scores{$key}[3]."\t".$scores{$key}[4]."\t".$scores{$key}[5]."\t".$scores{$key}[6]."\t".$scores{$key}[7]."\n";
    }
    else
    {
        print $key."\n";
    }
}
close(scores_target);

