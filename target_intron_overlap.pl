#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );
open(targets,"<", "targetpos_Ago25_323_41_19-HEK_dicer_results4_dicer005_ago005_drosha05_all.bed")or die("unable to open targets");

#open(targets,"<", "randomtargetpos_upto513.bed")or die("unable to open targets");
my %targets;
my %gtargets;
 while(my $line=<targets>)
{
chomp($line);
my @array=split(/\t/, $line);

push @{$targets{$array[0]}},[$array[1],$array[2],$array[3],$array[5],$array[6]];
push @{$gtargets{$array[3]}},[$array[0],$array[1],$array[2],$array[5],$array[6]];

}
 
 close(targets);
 
 
open(introns,"<", "introns1.bed")or die("unable to open introns");
my %introns;
my %gintrons;
 while(my $line=<introns>)
{
chomp($line);
my @array=split(/\t/, $line);

push @{$introns{$array[0]}},[$array[1],$array[2],$array[3],$array[5]];
push @{$gintrons{$array[3]}},[$array[0],$array[1],$array[2],$array[5]];

}
 
 close(introns);
 
 open(strict,"<", "targetgenes_nomir_ordered_strict.txt")or die("unable to open strict");
my %strict;
 while(my $line=<strict>)
{
chomp($line);
my @array=split(/\t/, $line);

$strict{$array[0]}="strict"

}
 
 close(strict);
 
  open(nstrict,"<", "targetgenes_nomir_ordered.txt")or die("unable to open nstrict");
my %nstrict;
 while(my $line=<nstrict>)
{
chomp($line);
my @array=split(/\t/, $line);

$nstrict{$array[0]}="nstrict"

}
 
 close(nstrict);
 
open(genes,"<", "HEK_dicer_results1_dicer005_ago005_drosha05_upregulated.txt")or die("unable to open genes");
my %genes;
<genes>;
 while(my $line=<genes>)
{
chomp($line);
my @array=split(/\t/, $line);

$genes{$array[0]}=[$array[2],$array[7]];
}
 
 close(genes); 
 
 
 my $key;
my @genesgenes;
foreach $key (keys %gtargets)
{
  for(my $u=0; $u<=$#{$gtargets{$key}}; $u++)
  {
      for(my $v=0; $v<=$#{$gintrons{$key}}; $v++)
      {

          if($gtargets{$key}[$u][1]>=$gintrons{$key}[$v][1]&&$gtargets{$key}[$u][2]<=$gintrons{$key}[$v][2])
          {
            $gtargets{$key}[$u][5]="INTRON";
           # print "INTRON\n";
            push @genesgenes, $key;
          }
      }
  }
}
 @genesgenes=uniq(@genesgenes);
 print $#genesgenes."\n";
   open(upint,">", "upint.txt")or die("unable to open upint");
   open(downint,">", "downint.txt")or die("unable to open downint");
  my $upinint=0;
  my $downinint=0;
    my $upnotinint=0;
  my $downnotinint=0;
  
  
  my $intstrict=0;
  my $intnonstrict=0;
  my $nonintstrict=0;
  my $nonintnonstrict=0;
  my $allup=0;
  my $alldown=0;

foreach $key (keys %genes)
{

  my $inintron=0;
  for(my $u=0; $u<=$#{$gtargets{$key}}; $u++)
  {
    if($gtargets{$key}[$u][5]eq"INTRON")
    {
     $inintron++;
    }
  }
  my $strict="miRNA evidence";
  if($strict{$key}eq"strict")
  {
    $strict="No miRNA evidence"
  }
  elsif($nstrict{$key}eq"nstrict")
  {
    $strict="WEAK miRNA evidence"
  }

  if($genes{$key}[0]>0)
  {
    $allup++;
  print upint $key."\t".$genes{$key}[1]."\t".$strict."\t".$inintron."\t".($#{$targets{$key}}+1-$inintron)."\n";
    if($inintron>0)
    {
      $upinint++;
      if($strict eq "miRNA evidence")
      {
        $intstrict++;
      }
      elsif($strict eq "WEAK miRNA evidence")
      {
        $intnonstrict++;
      }
    }
    else
    {
            $upnotinint++;

      if($strict eq "miRNA evidence")
      {
        $nonintstrict++;
      }
      elsif($strict eq "WEAK miRNA evidence")
      {
        $nonintnonstrict++;
      }
    }
  }
    if($genes{$key}[0]<0)
  {
    $alldown++;
  print downint $key."\t".$genes{$key}[1]."\t".$strict."\t".$inintron."\t".($#{$targets{$key}}+1-$inintron)."\n";
if($inintron>0)
    {
      $downinint++;
      
    }
    else
    {
            $downnotinint++;

    }
  }
}


print "UP: ".$allup." ".$upinint." ".$intstrict." ".$intnonstrict." ".$upnotinint." ".$nonintstrict." ".$nonintnonstrict."\n";
print "DOWN: ".$alldown." ".$downinint." ".$downnotinint."\n";
 my $trnainintron_all=0;
 my $trnainexon_all=0;
 my %trnatarget;
 my @genetest1;
 my @genetest2;
  my @genetest3;

foreach $key (keys %gtargets)
{
   my $trnainintron=0;
 my $trnainexon=0;
     my $r=0;

  for(my $u=0; $u<=$#{$gtargets{$key}}; $u++)
  {
    if($gtargets{$key}[$u][5]eq"INTRON")
    {
      $trnainintron++;
      $trnatarget{$gtargets{$key}[$u][4]}[0]++;
      if (exists$genes{$key})
      {
      push @genetest1,$key;
      $r=1;
      }
    }
     else
     {
      $trnatarget{$gtargets{$key}[$u][4]}[1]++;
            if (exists$genes{$key})
      {
 
            push @genetest2,$key;
      }

     }
  }
if($r==0)
{
             push @genetest3,$key;
 
}
}
my $size = keys %trnatarget;
print $size."\n";
my $overlap=0;
my $notarget=0;
my $onlyintron=0;
my $onlyexon=0;


foreach $key (keys %trnatarget)
{
  if($trnatarget{$key}[0]>0)
  {
    $trnainintron_all++;
  }

    if($trnatarget{$key}[1]>0)
  {
    $trnainexon_all++;
  }
  if($trnatarget{$key}[0]>0&&$trnatarget{$key}[1]>0)
  {
    $overlap++;
  }
  if($trnatarget{$key}[0]<=0&&$trnatarget{$key}[1]<=0)
  {
    $notarget++
  }
    if($trnatarget{$key}[0]>0&&$trnatarget{$key}[1]<=0)
  {
    $onlyintron++;
  }
  if($trnatarget{$key}[0]<=0&&$trnatarget{$key}[1]>0)
  {
    $onlyexon++;
  }
}
@genetest1=uniq(@genetest1);
@genetest2=uniq(@genetest2);
@genetest3=uniq(@genetest3);

print "TRNA ".$trnainintron_all." ".$trnainexon_all." ".$overlap." ".$notarget." ".$onlyintron." ".$onlyexon."\n";
print $#genetest1." ".$#genetest2." ".$#genetest3."\n";
close(upint);
close(downint); 
 