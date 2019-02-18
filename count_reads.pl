#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use Statistics::Basic qw(:all);

my @reads;

open(reads,"<","Dicer_rep1_mappings_trna_sort.sam") or die ("unable to open reads"); #ENSEMBL
my $q=0;
while (my $line=<reads>)
{
   # $q+=(1/13313248);
    #print $q."\n";
    
    my @array=split(/\t/,$line);
    push @reads,$array[0];
}
close(reads);

@reads=uniq(@reads);
print @reads."\n";
