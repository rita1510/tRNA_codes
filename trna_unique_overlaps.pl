#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);

open(fasta,"<", "AGO_Dicer_trna_7nt_allDicers_0_Ago1_1_1_1.fasta")or die("unable to open trnabed");

my %seq;
my $trna;
while(my $line=<fasta>)
{
chomp($line);
if($line=~">")
{
$trna=$line;
print $trna."\n";
}
#else
#{
my $line1=<fasta>;
chomp($line1);
    push @{$seq{$line1}},$trna;
#}
}
close (trnaseq);

open(file,">", "AGO_Dicer_trna_7nt_allDicers_0_Ago1_1_1_1_unique.fasta")or die("unable to open file");
my $key;
foreach $key (keys %seq)
{
print file $seq{$key}[0]."\n";
print file $key."\n";
}
close(file);

