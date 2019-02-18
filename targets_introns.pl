#!/usr/bin/perl

use strict;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

open(coord,"<", "HEK_dicer_results1_dicer005_ago005_drosha05_coord.txt")or die("unable to open coord");
my %coords;
 while(my $line=<coord>)
{
chomp($line);
my @array=split(/\t/, $line);
if(length($array[2])<4&&$array[2]ne"MT")
{
$coords{$array[6]}=["chr".$array[2],$array[4],$array[5],$array[3]];
}
elsif($array[2]eq"MT")
{
        $coords{$array[6]}=["chrM",$array[4],$array[5],$array[3]];

}
else
{
        #$coords{$array[6]}=[$array[2],$array[4],$array[5],$array[3]];

}
}

open(fch,"<", "HEK_dicer_results1_dicer005_ago005_drosha05_upregulated.txt")or die("unable to open fch");
my %fch;
 while(my $line=<fch>)
{
chomp($line);
my @array=split(/\t/, $line);
$fch{$array[0]}=[$array[2],$array[5],$array[6],$array[7]];
}
 
open(file1,"<", "../trna_mappings/AGO_Dicer_trna_7nt_Dicer1.fasta")or die("unable to open file1");
my %dicer1;
my %rdicer1;

while(my $line=<file1>)
{
        chomp($line);

    my @array=split(/ /, $line);
        my @array1=split(/>hg38_tRNAs_/, $array[0]);
#print $array1[1]."\n";
    my $line1=<file1>;
    chomp($line1);

    push @{$dicer1{$line1}}, [$array[3],$array[4],$array1[1]];#ago-count, dicer counts,tRNA name
   # print $#{$dicer1{$line1}}."\n";
    push @{$rdicer1{$array1[1]}}, $line1;
}
close (file1);

open(file2,"<", "../trna_mappings/AGO_Dicer_trna_7nt_Dicer2.fasta")or die("unable to open file2");
my %dicer2;
my %rdicer2;

while(my $line=<file2>)
{
        chomp($line);

    my @array=split(/ /, $line);
        my @array1=split(/>hg38_tRNAs_/, $array[0]);

    my $line1=<file2>;
    chomp($line1);

    push @{$dicer2{$line1}}, [$array[3],$array[4],$array1[1]];#ago-count, dicer counts
    push @{$rdicer2{$array1[1]}}, $line1;

}
close (file2);

open(file3,"<", "../trna_mappings/AGO_Dicer_trna_7nt_Dicer3.fasta")or die("unable to open file3");
my %dicer3;
my %rdicer3;

while(my $line=<file3>)
{
        chomp($line);

    my @array=split(/ /, $line);
            my @array1=split(/>hg38_tRNAs_/, $array[0]);

    my $line1=<file3>;
    chomp($line1);

    push @{$dicer3{$line1}}, [$array[3],$array[4],$array1[1]];#ago-count, dicer counts
    push @{$rdicer3{$array1[1]}}, $line1;

}
close (file3);
open(score,"<", "../trna_mappings/Score_S1_S2_S3_D1_D2_D3_AD_7nt_allDicers_0_Ago1_1_1_1")or die("unable to open score");
my %score123;
while(my $line=<score>)
{
        chomp($line);

    my @array=split(/\t/, $line);
    my $s1;
    my $s2;
    my $s3;
if($array[5]>0&&$array[2]>0)
{
        $s1=log($array[5]/$array[2])/log(2);
}
else
{
      $s1="NA"  
}
if($array[6]>0&&$array[3]>0)
{
        $s2=log($array[6]/$array[3])/log(2);
}
else
{
      $s2="NA"  
}
if($array[7]>0&&$array[4]>0)
{
        $s3=log($array[7]/$array[4])/log(2);
}
else
{
      $s3="NA"  
}
        $score123{$array[0]}=[$s1,$s2,$s3];
}
close (file3);

open(miranda,"<", "miranda-AGO_Dicer_trna_7nt_allDicers-HEK_dicer_results4_dicer005_ago005_drosha05_all.txt")or die("unable to open miranda");
my %miranda;
my $string="";
my $string1="";

my %targets;
my %sequence;
my %sequence1;

my %score;
my %Energy;
while(my $line=<miranda>)
{
chomp($line);
        if($line=~"Query")
        {
        my @array=split(/'/, $line);

          #  print $array[1]."\n";
            $string = uc reverse $array[1];
            $string =~ tr/35\-\ //d;
        }
         if($line=~"Ref:")
        {
        my @array=split(/'/, $line);

          #  print $array[1]."\n";
            $string1 = uc reverse $array[1];
            $string1 =~ tr/35\-\ //d;
        }
        if($line=~">>")
        {
        my @array=split(/\t/, $line);
        my $gene=$array[1];
        my @array2=split(/_/, $array[0]);
        
        my $trna=$array2[2];

        my @array1=split(" ",$array[9]);
            
            @{$targets{$gene."_".$trna}}=uniq(@{$targets{$gene."_".$trna}});

            push @{$miranda{$gene}},$trna;
            @{$miranda{$gene}}=uniq(@{$miranda{$gene}});
                
            push @{$sequence{$gene."_".$trna}},$string;
            push @{$sequence1{$gene."_".$trna}},$string1;

            #print $string1." 2\n";
            #print length($string1)."\n";

            @{$sequence{$gene."_".$trna}}=uniq(@{$sequence{$gene."_".$trna}});
            for(my $u=0; $u<=$#array1; $u++)
            {
                if($coords{$gene}[3]==1)
                {
                push @{$targets{$string."_".$gene}},($coords{$gene}[1]+$array1[$u]);
                }
                elsif($coords{$gene}[3]==-1)
                {
                push @{$targets{$string."_".$gene}},($coords{$gene}[2]-$array1[$u]);
                }
            }
            
            @{$targets{$string."_".$gene}}=uniq(@{$targets{$string."_".$gene}});
        $score{$gene."_".$trna."_".$string}=$array[2];
                $Energy{$gene."_".$trna."_".$string}=$array[3];

        }
        if($line=~"Complete")
        {
        $string="";
        }
        
}
close (miranda);

my $key;
open(geneswitht1,">", "table_allfields_genes_with_mirandatargets_0_Ago25_323_41_19-HEK_dicer_results4_dicer005_ago005_drosha05_all")or die("unable to open geneswitht");
open(targetpos,">", "targetpos_Ago25_323_41_19-HEK_dicer_results4_dicer005_ago005_drosha05_all.bed")or die("unable to open geneswitht");

foreach $key (keys %miranda)
{
        

        

    my $min0 = min @{$Energy{$key."_".$miranda{$key}[0]}};
    my $max0 = max @{$score{$key."_".$miranda{$key}[0]}};

    if(exists $score123{$sequence{$key."_".$miranda{$key}[0]}[0]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[0]}[0]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[0]}[0]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[0]}[0]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[0]}[0]}[2];
        
        
    print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[0]."\t".$sequence{$key."_".$miranda{$key}[0]}[0]."\t".$Energy{$key."_".$miranda{$key}[0]."_".$sequence{$key."_".$miranda{$key}[0]}[0]}."\t".$score{$key."_".$miranda{$key}[0]."_".$sequence{$key."_".$miranda{$key}[0]}[0]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[0]}[0]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[0]}[0]}[0][1]."\t".$dicer2{$sequence{$key."_".$miranda{$key}[0]}[0]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[0]}[0]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";
    my $tk=$sequence{$key."_".$miranda{$key}[0]}[0]."_".$key;
    print $tk." 1\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[0]}[0])-1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[0]}[0])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[0]}[0]."\n";
        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[0]}[0])+1)."\t".($targets{$tk}[$k]-1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[0]}[0])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[0]}[0]."\n";

        }
    }
        print geneswitht1 "\n";
    for(my $u=1; $u<$#{$sequence{$key."_".$miranda{$key}[0]}}; $u++)
    {
        if(exists $score123{$sequence{$key."_".$miranda{$key}[0]}[$u]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[0]}[$u]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[0]}[$u]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[0]}[$u]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[0]}[$u]}[2];

        print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[0]."\t".$sequence{$key."_".$miranda{$key}[0]}[$u]."\t".$Energy{$key."_".$miranda{$key}[0]."_".$sequence{$key."_".$miranda{$key}[0]}[$u]}."\t".$score{$key."_".$miranda{$key}[0]."_".$sequence{$key."_".$miranda{$key}[0]}[$u]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[0]}[$u]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[0]}[$u]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[0]}[$u]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[0]}[$u]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";
    my $tk=$sequence{$key."_".$miranda{$key}[0]}[$u]."_".$key;
    print $tk." 2\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[0]}[$u])-1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[0]}[$u])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[0]}[$u]."\n";

        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[0]}[$u]))."\t".($targets{$tk}[$k]+1)."\t";
                print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[0]}[$u])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[0]}[$u]."\n";

        }
    }
        print geneswitht1 "\n";
    }
    if($#{$sequence{$key."_".$miranda{$key}[0]}}>0)
    {
        
        if(exists $score123{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}[2];
        
       print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[0]."\t".$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]."\t".$Energy{$key."_".$miranda{$key}[0]."_".$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}."\t".$score{$key."_".$miranda{$key}[0]."_".$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}[0][1]."\t".$dicer2{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";
    my $tk=$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]."_".$key;
    print $tk." 3\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}])-1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]."\n";

        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence{$key."_".$miranda{$key}[0]}[$#{$sequence1{$key."_".$miranda{$key}[0]}}])+1)."\t".($targets{$tk}[$k]+1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence{$key."_".$miranda{$key}[0]}[$#{$sequence1{$key."_".$miranda{$key}[0]}}])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[0]}[$#{$sequence{$key."_".$miranda{$key}[0]}}]."\n";

        }
    }
        print geneswitht1 "\n";
    }
  for(my $u=1; $u<$#{$miranda{$key}}; $u++)
     {
    my $min = min @{$Energy{$key."_".$miranda{$key}[$u]}};
    my $max = max @{$score{$key."_".$miranda{$key}[$u]}};
    if($#{$sequence{$key."_".$miranda{$key}[$u]}}>0)
    {
        if(exists $score123{$sequence{$key."_".$miranda{$key}[$u]}[0]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[$u]}[0]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[$u]}[0]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[$u]}[0]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[$u]}[0]}[2];
        print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[$u]."\t".$sequence{$key."_".$miranda{$key}[$u]}[0]."\t".$Energy{$key."_".$miranda{$key}[$u]."_".$sequence{$key."_".$miranda{$key}[$u]}[0]}."\t".$score{$key."_".$miranda{$key}[$u]."_".$sequence{$key."_".$miranda{$key}[$u]}[0]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$u]}[0]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$u]}[0]}[0][1]."\t".$dicer2{$sequence{$key."_".$miranda{$key}[$u]}[0]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[$u]}[0]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";
    my $tk=$sequence{$key."_".$miranda{$key}[$u]}[0]."_".$key;
    print $tk." 4\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$u]}[0])-1)."\t";
                print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$u]}[0])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[$u]}[0]."\n";

        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$u]}[0])+1)."\t".($targets{$tk}[$k]+1)."\t";

        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$u]}[0])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[$u]}[0]."\n";
        }
    }
        print geneswitht1 "\n";
    }
     
    for(my $v=1; $v<$#{$sequence{$key."_".$miranda{$key}[$u]}}; $v++)
    {
        if(exists $score123{$sequence{$key."_".$miranda{$key}[$u]}[$v]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[$u]}[$v]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[$u]}[$v]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[$u]}[$v]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[$u]}[$v]}[2];
        
        print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[$u]."\t".$sequence{$key."_".$miranda{$key}[$u]}[$v]."\t".$Energy{$key."_".$miranda{$key}[$u]."_".$sequence{$key."_".$miranda{$key}[$u]}[$v]}."\t".$score{$key."_".$miranda{$key}[$u]."_".$sequence{$key."_".$miranda{$key}[$u]}[$v]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$u]}[$v]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$u]}[$v]}[0][1]."\t".$dicer2{$sequence{$key."_".$miranda{$key}[$u]}[$v]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[$u]}[$v]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";
    my $tk=$sequence{$key."_".$miranda{$key}[$u]}[$v]."_".$key;
    print $tk." 5\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$u]}[$v])-1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$u]}[$v])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[$u]}[$v]."\n";

        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$u]}[$v])+1)."\t".($targets{$tk}[$k]+1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$u]}[$v])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[$u]}[$v]."\n";

        }
    }
        print geneswitht1 "\n";
    }
    
    if(exists $score123{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}[2];
     print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[$u]."\t".$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]."\t".$Energy{$key."_".$miranda{$key}[$u]."_".$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}."\t".$score{$key."_".$miranda{$key}[$u]."_".$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}[0][1]."\t".$dicer2{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";
    my $tk=$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]."_".$key;
    print $tk." 6\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}])-1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]."\n";

        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}])+1)."\t".($targets{$tk}[$k]+1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[$u]}[$#{$sequence{$key."_".$miranda{$key}[$u]}}]."\n";

        }
    }
        print geneswitht1 "\n";
    }
 if($#{$miranda{$key}}>0)
{
  
        my $minu = min @{$Energy{$key."_".$miranda{$key}[$#{$miranda{$key}}]}};
    my $maxu = max @{$score{$key."_".$miranda{$key}[$#{$miranda{$key}}]}};
    
    if($#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}>0)
    {
        if(exists $score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}[2];
        
    print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[$#{$miranda{$key}}]."\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]."\t".$Energy{$key."_".$miranda{$key}[$#{$miranda{$key}}]."_".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}."\t".$score{$key."_".$miranda{$key}[$#{$miranda{$key}}]."_".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}[0][1]."\t".$dicer2{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";
my $tk=$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]."_".$key;
print $tk." 7\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0])-1)."\t";
         print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]."\n";

        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0])+1)."\t".($targets{$tk}[$k]+1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[0]."\n";

        }
    }
        print geneswitht1 "\n";
    }
    for(my $v=1; $v<$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}; $v++)
    {
         if($#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}>0)
        {
        if(exists $score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}[2];
        
 
        print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[$#{$miranda{$key}}]."\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]."\t".$Energy{$key."_".$miranda{$key}[$#{$miranda{$key}}]."_".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}."\t".$score{$key."_".$miranda{$key}[$#{$miranda{$key}}]."_".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}[0][1]."\t".$dicer2{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";
my $tk=$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]."_".$key;
print $tk." 8\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v])-1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]."\n";

        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v])+1)."\t".($targets{$tk}[$k]+1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$v]."\n";

        }
    }
        print geneswitht1 "\n";
        }
    }
            if(exists $score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]})
        {
        }
        else
        {
                $score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}=["no score","no score","no score"];
        }
        my $s1=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}[0];
        my $s2=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}[1];
        my $s3=$score123{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}[2];
        
    print geneswitht1 $key."\t".$fch{$key}[3]."\t".$fch{$key}[0]."\t".$fch{$key}[1]."\t".$fch{$key}[2]."\t".$miranda{$key}[$#{$miranda{$key}}]."\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]."\t".$Energy{$key."_".$miranda{$key}[$#{$miranda{$key}}]."_".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}."\t".$score{$key."_".$miranda{$key}[$#{$miranda{$key}}]."_".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}[0][0]."\t".$dicer1{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}[0][1]."\t".$dicer2{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}[0][1]."\t".$dicer3{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]}[0][1]."\t".$s1."\t".$s2."\t".$s3."\t".$coords{$key}[3]."\t";    
    my $tk=$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]."_".$key;
    print $tk." 9\n";
    print $#{$miranda{$key}}."\n";
    print $#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}."\n";
    for(my $k=0; $k<=$#{$targets{$tk}}; $k++)
    {
        if($coords{$key}[3]==1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-1)."\t".($targets{$tk}[$k]+length($sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}])-1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-2)."\t".($targets{$tk}[$k]+length($sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}])-2)."\t".$key."\t0\t+\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]."\n";

        }
 
        elsif($coords{$key}[3]==-1)
        {
        print geneswitht1 $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}])+1)."\t".($targets{$tk}[$k]+1)."\t";
        print targetpos $coords{$key}[0]."\t".($targets{$tk}[$k]-length($sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence1{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}])+1)."\t".($targets{$tk}[$k]+1)."\t".$key."\t0\t-\t".$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}[$#{$sequence{$key."_".$miranda{$key}[$#{$miranda{$key}}]}}]."\n";

        }
    }
 
        print geneswitht1 "\n";
 }
# for(my $u=0; $u<$#{$targets{$key}}; $u++)
    #{
   #     print geneswitht $targets{$key}[$u]." ";
    #}
    #print geneswitht $targets{$key}[$#{$targets{$key}}]."\t";
    #for(my $u=0; $u<=$#{$sequence{$key}}; $u++)
    #{
    #    print geneswitht $sequence{$key}[$u]." ";
    #}
    #print geneswitht $sequence{$key}[$#{$sequence{$key}}]."\n";
        
}
close(targetpos);
close(geneswitht1)