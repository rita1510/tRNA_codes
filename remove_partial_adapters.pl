#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);

my $adapter="GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG";
open(fastq,"<", "CAWCVACXX_D1_17s003026-1-1_Davis_lane517s003026_sequence_trimmed.fq")or die("unable to open fastq");
open(fastq_removed,">", "CAWCVACXX_D1_17s003026-1-1_Davis_lane517s003026_sequence_trimmed_removed.fq")or die("unable to open fastq_removed");

while(my $line1=<fastq>)
{
   $line2=<fastq>;
   $line3=<fastq>;
   $line4=<fastq>;
 chomp($line1);
  chomp($line2);
 chomp($line3);
 chomp($line1);
 print $line1."\n";
print $line2."\n";  
 print $line3."\n";  
 print $line4."\n";  

}