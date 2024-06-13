#!/usr/bin/env perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
use FindBin qw($Bin);
my $cwd = getcwd();

### Initialization of software and database ###
my $nogDir = "/data/database/micro/eggNOG";
my $nognames = "$nogDir/NOG.annotations.xls";
my $gi2nog = "$nogDir/NOG.members.tsv";
my $Rscript="/home/zjming/anaconda3/bin/Rscript";


### usage ###
sub Usage {
  print STDERR << "USAGE"; 
  Description: This script is aimed at identifying the COG class corresponding to query genes or proteins
  based on BLAST suite output result 
  Author: longsheng.xing\@allwegene.com
  Version: 1.0
=====================================================================================
Options
-i   --input   input sequence file in FASTA-format
-t   --thread  number of CPUs to use (default: 8)
-o   --out     output directory
-h   --help    print this help info
-k   --project name
=====================================================================================

USAGE

  exit();

}

my $blastOut;
my $outdir;
my $pro;
my $thread = 0;
my $help;
GetOptions(
    "input|i=s" => \$blastOut,
    "out|o=s" => \$outdir,
    "thread|t=i" => \$thread,
    "help|h" => \$help,
    "k=s" => \$pro,
	);
           
if( !defined($blastOut) || !defined($outdir) || 
    $help ) {
  &Usage();	
}

$thread = $thread == 0 ? 8 : $thread;

if( !-e $outdir ) {
  mkdir $outdir;
}
$outdir =~ s/\/$//;


#C(K)OG ID to alpgabet and annotation
open NAMES,"<$nognames" || die "Cannot open this file!$!";
<NAMES>;
my %nog2alphabet = ();
my %nog2name = ();
while(<NAMES>) {
  chomp;
  my @line = split /\t/,$_;
  if($line[0]!~"#")
  {
    $nog2alphabet{$line[0]} = $line[1];	
    $nog2name{$line[0]} = $line[2];
  }
}

close NAMES;

### reflection of eggNOG ID to C(K)OG ID ###
open NOG,"<$gi2nog" || die "Cannot open this file!$!";
my %GI2NOG;
while(<NOG>)
{
  chomp;
  my @line = split /\t/,$_;
  my @tp = split (/,/,$line[5]);
  for(my $i=0;$i<=$#tp;$i++)
  {
    $GI2NOG{$tp[$i]} = $line[1];
  }
}
close NOG;

#### core processing stage ###
open IN,"<$blastOut" || die "Cannot open this file!$!";

my %transcript2NOG = ();
my %transcript2FUN = ();
my %transcript2DESC = ();
while(<IN>) {
  chomp;
  my @line = split /\t/,$_;
  if(exists $GI2NOG{$line[1]})
  {
    $transcript2FUN{$line[0]} = $nog2alphabet{$GI2NOG{$line[1]}};
    $transcript2NOG{$line[0]} = $GI2NOG{$line[1]};
    $transcript2DESC{$line[0]} = $nog2name{$GI2NOG{$line[1]}};
  }
}

close IN;

### output COG annotation result for each transcript ###
my %NOG2count;
open OUT,">$outdir/$pro.anno.xls" || die "Cannot write to this file!$!";
print OUT "Cluster_id\tAbbreviation\teggNOG_ID\teggNOG_Description\n";
foreach my $item(keys %transcript2FUN)
{
  my $temp = $transcript2FUN{$item};
  if(length($temp) == 1)
  {
    $NOG2count{$temp}++;
  }
  else
  {
    my $n = length($temp);
    foreach (0 .. $n - 1)
    {
      my $substr = substr($temp,$_,1);
      $NOG2count{$substr}++;	
    }	
  }
  if($temp)
  {
    print OUT "$item\t$temp\t$transcript2NOG{$item}\t$transcript2DESC{$item}\n";
  }
}

close OUT;

#######################################
### statistics of all the NOG items ###
#######################################
my $NOG_summary = "$outdir/$pro.summary.xls";
open SUMMARY,">$NOG_summary" || die "Cannot write to this file!$!";

print SUMMARY "eggNOG\tnumber\n";
my @items = ('A','B','C','D','E','F','G','H','I','J','K','L','M',
             'N','O','P','Q','R','S','T','U','V','W','X','Y','Z');

foreach my $i(0 .. $#items) {
  my $cnt = 0;
  if (defined $NOG2count{$items[$i]}) {
    $cnt  = $NOG2count{$items[$i]};
  }
  print SUMMARY "$items[$i]\t$cnt\n";
}

close SUMMARY;


system("paste $NOG_summary $Bin/NOG_type > $outdir/$pro\.eggNOGsummary.xls");
system "mv $outdir/$pro\.eggNOGsummary.xls $NOG_summary";
system("export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH && Rscript $Bin/NOG.bar.R $NOG_summary $outdir $pro");
