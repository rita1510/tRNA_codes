#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );
open(mir,"<", "../trna_mappings/mirtar.tab")or die("unable to open mir");
my %mir;
while(my $line=<mir>)
{
chomp($line);
my @array=split(/\t/, $line);

#if($array[7]!~"Weak")
#{
    $mir{$array[3]}="mir";
#}

}
close(mir);

print keys(%mir)."\n";

open(geneswitht1,"<", "table_allfields_genes_with_mirandatargets_0_Ago25_323_41_19_HEK_dicer_results4_dicer005_ago005_drosha05_all")or die("unable to open geneswitht");
my @genes1;
my @genes;

while(my $line=<geneswitht1>)
{
chomp($line);
my @array=split(/\t/, $line);
if($array[2]>0)
{
    push @genes1,$array[0]."_".$array[1]."_".$array[2];
}
}
close(geneswitht1);
    @genes1=uniq(@genes1);

for(my $u=0; $u<=$#genes1; $u++)
{
    my @array=split(/_/,$genes1[$u]);
        print $array[0]."\n";

    $genes[$u]=[$array[0],$array[1],$array[2]];
}

open(genesnomir1,">", "targetgenes_nomir_ordered.txt")or die("unable to open genesnomir1");
my @genesnomir;
for(my $u=0; $u<=$#genes; $u++)
{
            print $genes[$u][1]." ".$mir{$genes[$u][1]}."\n";

    if($mir{$genes[$u][1]}eq"mir")
    {
    }
    elsif($genes[$u][1]ne"")
    {
        print $genes[$u][1]."\n";
        push @genesnomir, [$genes[$u][0],$genes[$u][1],$genes[$u][2]];
    }
}
@genesnomir=sort{$b->[2]cmp$a->[2]}@genesnomir;

for(my $u=0; $u<=$#genesnomir; $u++)
{
    print $genesnomir[$u][0]." ".$genesnomir[$u][1]." ".$genesnomir[$u][2]."\n";
    
    print genesnomir1 $genesnomir[$u][0]."\t".$genesnomir[$u][1]."\t".$genesnomir[$u][2]."\n";
}

close(genesnomir1);