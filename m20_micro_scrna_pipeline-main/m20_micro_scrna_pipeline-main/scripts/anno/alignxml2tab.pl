#!/usr/bin/perl
#xiaoyue.wang
=head1 Name

	alignxml2tab.pl -- parse the BLAST/diamond result and convert to tab format.

=head1 Description

	Parses XML format into the TAB format.
	TAB format header:
	1:Query_id  2:Query_length  3:Query_start  4:Query_end  5:Subject_id  6:Subject_length  7:Subject_start  
	8:Subject_end  9:Identity  10:Gap  11:Align_length  12:Score  13:E_value  14:Description

=head1 Version

	Author : wangxiaoyue   email:1557182368\@qq.com
	Version: 1.0	       Date : 2017-10-30
	
=head1 Usage

  	perl alignxml2tab.pl [options] input_file
	-nohead               If there is this parameter, will not print the header.
	-tophit     <num>     integer(default:1), to set how many subjects for a query to be displayed. 
	-topmatch   <num>     integer(default:1), to set suits(results of one subject match one query) to be displayed. 
	-eval       <num>     float or exponent(default:1e-5), to filter the alignments which worse than the E-value cutoff.
	-nofilter             Not filter, will print all align record.
	-m8       <outfile>   [optional], If there is this parameter,program will output blast M8 format
	-help                 output help information to screen.

=head1 Exmple
	
	1. perl alignxml2tab.pl test.nr.blast.xml > test.nr.blast.xls
	2. perl alignxml2tab.pl -tophit 2 -topmatch 1 -nohead test.nr.blast.xml > test.nr.blast.xls

=cut

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Data::Dumper;

my ($Nohead,$Tophit,$Topmatch,$Eval,$nofilter,$m8);
my ($Help);
GetOptions(
	"nohead"=>\$Nohead,
	"tophit:i"=>\$Tophit,
	"topmatch:i"=>\$Topmatch,
	"eval:f"=>\$Eval,
	"nofilter"=>\$nofilter,
	"m8:s"=>\$m8,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV==0 || $Help);
my $blast_file = shift;
$blast_file=File::Spec->rel2abs($blast_file);
$Tophit ||=1;
$Topmatch ||=1;
$Eval ||= 1e-5;

##convert blast raw result to tabular format
if(defined $m8){
	$m8=File::Spec->rel2abs($m8);
	&parse_blast($blast_file,$Tophit,$Topmatch,$Eval,$Nohead,$m8);
}else{
	&parse_blast($blast_file,$Tophit,$Topmatch,$Eval,$Nohead);
}

####################################################
################### Sub Routines ###################
####################################################


##parse the BLAST files, and output in tabular formats 
####################################################
sub parse_blast{
	my ($file,$tophit,$topmatch,$eval,$nohead,$m8file) = @_;
	$/="</Iteration>\n";
	open (BLAST,"<$file") || die ("Error: Could not open the file:$file!\n");
	if(defined $m8file && $m8file ne ""){
		open OUT,">$m8file" || die $!;
		print OUT "#Query_id\tSubject_id\t%identity\tAlign_length\tMismatches\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE_value\tScore\n" unless(defined $nohead);
	}
	print "#Query_id\tQuery_length\tQuery_start\tQuery_end\tSubject_id\tSubject_length\tSubject_start\tSubject_end\tIdentity\tGap\tAlign_length\tScore\tE_value\tDescription\n" unless(defined $nohead);
	while(<BLAST>){
		next if(/^\s*$/);
		next if($_!~/<Iteration_query-def>/);
		my ($query,$query_len);
		my @b=split /<\/?Iteration_hits>\n/,$_;
		$query=$1 if($b[0]=~/<Iteration_query-def>(.*)<\/Iteration_query-def>/);
		$query_len=$1 if($b[0]=~/<Iteration_query-len>(.*)<\/Iteration_query-len>/);
		if($b[1]=~/^\s*$/){
			print STDERR "$query : No hits found in DataBase!\n";
			next;
		}
		#print "################################\n$b[1]#####################################\n";
		my @hit=split /<\/Hit>\n/,$b[1];
		for my $h(@hit){
			my @hsp=split /<Hsp>/,$h;
			my ($hit_num,$hitID,$hit_len,$annotation);
			$hit_num=$1 if($hsp[0]=~/<Hit_num>(.*)<\/Hit_num>/);
			if($hsp[0]=~/<Hit_id>(.*)<\/Hit_id>/){
				$hitID = $1;
				if($hitID=~/^sp\|/){
					$hitID=~s/sp\|//g;
					$hitID=~s/\|$//g;
				}
			}
			if($hsp[0]=~/<Hit_def>(.*)<\/Hit_def>/){
				my $raw_anno=$1;
				$raw_anno=$hitID if($raw_anno eq "");
				if($raw_anno=~/\&gt;/){
					$annotation=(split /&gt;/, $raw_anno)[0];
                }elsif($raw_anno=~/\x01/ && $raw_anno=~/\#k__/){
                    $annotation=(split /\x01/,$raw_anno)[0]."#k__".(split /\#k__/,$raw_anno)[1];
				}elsif($raw_anno=~/\x01/){
					$annotation=(split /\x01/,$raw_anno)[0];
				}elsif($raw_anno=~/gi\|/){
					$annotation=(split /gi\|/,$raw_anno)[0];
				}else{
					$annotation=$raw_anno;
				}
			}
			if($hitID=~/^gnl\|BL_ORD_ID\|/){
				$hitID=$annotation;
			}
			$hit_len=$1 if($hsp[0]=~/<Hit_len>(.*)<\/Hit_len>/);
			for(my $i=1;$i<@hsp;$i++){
				my $match_num=$1 if($hsp[$i]=~/<Hsp_num>(.*)<\/Hsp_num>/);
				my $bits=int($1) if($hsp[$i]=~/<Hsp_bit-score>(.*)<\/Hsp_bit-score>/);
				my $e_value=$1 if($hsp[$i]=~/<Hsp_evalue>(.*)<\/Hsp_evalue/);
				my $qFrom=$1 if($hsp[$i]=~/<Hsp_query-from>(.*)<\/Hsp_query-from/);
				my $qTo=$1 if($hsp[$i]=~/<Hsp_query-to>(.*)<\/Hsp_query-to/);
				my $hFrom=$1 if($hsp[$i]=~/<Hsp_hit-from>(.*)<\/Hsp_hit-from/);
				my $hTo=$1 if($hsp[$i]=~/<Hsp_hit-to>(.*)<\/Hsp_hit-to/);
				my $identity=$1 if($hsp[$i]=~/<Hsp_identity>(.*)<\/Hsp_identity/);
				my $gap=$1 if($hsp[$i]=~/<Hsp_gaps>(.*)<\/Hsp_gaps/);
				my $length=$1 if($hsp[$i]=~/<Hsp_align-len>(.*)<\/Hsp_align-len/);
				if(defined $query && defined $query_len && defined $qFrom && defined $qTo && defined $hitID && defined $hit_len && defined $hFrom && defined $hTo && defined $identity && defined $gap && defined $length && defined $bits && defined $e_value && defined $annotation){
					my $percent1 = sprintf("%.2f", $identity*100/$length);
					my $percent = "$identity\/$length\($percent1\)";
					my $mismach=$length-$identity;
					if(defined $nofilter){
						print "$query\t$query_len\t$qFrom\t$qTo\t$hitID\t$hit_len\t$hFrom\t$hTo\t$percent\t$gap\t$length\t$bits\t$e_value\t$annotation\n";
						print OUT "$query\t$hitID\t$percent1\t$length\t$mismach\t$gap\t$qFrom\t$qTo\t$hFrom\t$hTo\t$e_value\t$bits\n" if(defined $m8file && $m8file ne "");
					}elsif($hit_num<=$tophit && $match_num<=$topmatch && $e_value<=$eval){
						print "$query\t$query_len\t$qFrom\t$qTo\t$hitID\t$hit_len\t$hFrom\t$hTo\t$percent\t$gap\t$length\t$bits\t$e_value\t$annotation\n";
						print OUT "$query\t$hitID\t$percent1\t$length\t$mismach\t$gap\t$qFrom\t$qTo\t$hFrom\t$hTo\t$e_value\t$bits\n" if(defined $m8file && $m8file ne "");
					}
				}else{
					die "<BLAST> $. line exists null value: $!\n";
				}
			}
		}
	}
	close BLAST;
	close OUT;
}
