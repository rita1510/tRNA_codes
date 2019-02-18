#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );
open(fch,"<", "HEK_dicer_results1_dicer005_ago005_drosha05_upregulated.txt")or die("unable to open fch");
my %fch;
 while(my $line=<fch>)
{
chomp($line);
my @array=split(/\t/, $line);
$fch{$array[0]}=[$array[2],$array[5],$array[6],$array[7]];
}

open(miranda,"<", "miranda-AGO_Dicer_trna_7nt_allDicers-HEK_dicer_results4_dicer005_ago005_drosha05_all.txt")or die("unable to open miranda");
my %miranda;
my $string="";
my %targets;
my %sequence;
my %score;
my %Energy;
my @trnas;
while(my $line=<miranda>)
{
chomp($line);
        if($line=~"Query")
        {
        my @array=split(/'/, $line);

          #  print $array[1]."\n";
            $string = uc reverse $array[1];
            $string =~ tr/5\-\ //d;
        }
        if($line=~">>")
        {
        my @array=split(/\t/, $line);
        my $gene=$array[1];
        my $trna=$array[0];
        my @array1=split(" ",$array[9]);
        my $targets=$array[9];
            for(my $u=0; $u<=$#array1; $u++)
            {
           push @{$targets{$gene}},$array1[$u];
            }
            @{$targets{$gene}}=uniq(@{$targets{$gene}});
                        $trna =~ tr/>//d;
push @trnas, $trna;
            push @{$miranda{$gene}},$trna;
            @{$miranda{$gene}}=uniq(@{$miranda{$gene}});
                
            push @{$sequence{$gene}},$string;
            @{$sequence{$gene}}=uniq(@{$sequence{$gene}});
        push @{$score{$gene}}, $array[2];
                push @{$Energy{$gene}}, $array[3];

        }
        if($line=~"Complete")
        {
        $string="";
        }
        
}
close (miranda);

@trnas=uniq(@trnas);
print $#trnas."\n";
die;

my $key;
open(geneswitht,">", "genes_with_mirandatargets_0_Ago25_323_41_19-HEK_dicer_results4_dicer005_ago005_drosha05_all")or die("unable to open geneswitht");

foreach $key (keys %miranda)
{
    my $min = min @{$Energy{$key}};
    my $max = max @{$score{$key}};

    print geneswitht $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$min."\t".$max."\t";

    for(my $u=0; $u<$#{$targets{$key}}; $u++)
    {
        print geneswitht $targets{$key}[$u]." ";
    }
    print geneswitht $targets{$key}[$#{$targets{$key}}]."\t";
    for(my $u=0; $u<$#{$miranda{$key}}; $u++)
    {
        print geneswitht $miranda{$key}[$u]." ";
    }
    print geneswitht $miranda{$key}[$#{$miranda{$key}}]."\t";


    for(my $u=0; $u<=$#{$sequence{$key}}; $u++)
    {
        print geneswitht $sequence{$key}[$u]." ";
    }
    print geneswitht $sequence{$key}[$#{$sequence{$key}}]."\n";

}

