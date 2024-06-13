#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin);
use File::Spec;


my($pep,$outdir,$card,$type,$threads,$evalue,$sample_name);
GetOptions(
    "p:s"=>\$pep,
    "o:s"=>\$outdir,
    "s:s"=>\$sample_name,
    "d:s"=>\$card,
    "type:s"=>\$type,
    "t:i"=>\$threads,
    "e:s"=>\$evalue,
);

if (!$pep || !$type || !$sample_name){
    &usage();
}

$card||="$Bin/protein_fasta_protein_homolog_model.fasta";
my $blast="/home/sunkun/Miniconda3/bin/";
$threads||=10;
$evalue||="1e-10";
$pep=File::Spec->rel2abs($pep);
$outdir ||= getcwd;
$outdir=File::Spec->rel2abs($outdir);
$card=File::Spec->rel2abs($card);


sub usage{
	print qq{
This script was used to annotated card.
usage:
perl $0 -p protein.fasta -o /path/xx/xx -d /xx/xx/prot
  -p          the query gene fasta [required]
  -s          sample_name / outfile_prefix. [required]
  -o          the output directory. [default: ./]
  -d          the database of card.
              /public/flash/works/Database/CARD/protein.(prot)
              /public/flash/works/Database/CARD/nucleotide.(nucl)
              [default: prot]
  -type       the file type you input.(prot or nucl?) [required]
  -t          the threads num set. [default: 10]
  -e          the evalue of blast. [default: 1e-10]

};
    exit;
}

my $out_name=(split /\//,$pep)[-1];
(-d $outdir) || system "mkdir -p $outdir/";
my $card_type=(split /\//,$card)[-1];
if($type eq 'prot'){
	if($card=~/nucl/){$blast.="tblastn";}else{$blast.="blastp";}
}elsif($type eq 'nucl'){
	if($card=~/nucl/){$blast.="blastn";}else{$blast.="blastx";}
}else{
	print "please input the correct type of your input file.\n";
	&usage(); 
}

print "begining blast,please wait .......\n";
system "$blast -query $pep -db $card -num_threads $threads -max_target_seqs 1 -evalue $evalue -outfmt \"6 qseqid sseqid salltitles qlen qstart qend slen sstart send evalue bitscore pident\" -out $outdir/$sample_name\.blast.xls";

print "blast is over,result at $outdir ,the card_result files\n\n";
print "clear up the final result now,please wait.......\n";

my %hash=();
open IN,"/home/sunkun/database/scMicro_refdata/CARD/aro_index.tsv" or die "please check out the aro_index.tsv\n";
while(<IN>){
	chomp;
	my ($ARO_Accession,$ARO_Name,$AMR_Gene_Family,$Drug_Class,$Resistance_Mechanism)=((split /\t/,$_)[0,5,8,9,10]);
    my @Name=split /\s+/, $ARO_Name;
    my $ARON='';
    if(@Name>3){
         $ARON=$Name[0]."_".$Name[1]."_".$Name[2];
     }else{
         $ARON=join("_",@Name);
     }
	$hash{$ARO_Accession}="$ARON\t$AMR_Gene_Family\t$Drug_Class\t$Resistance_Mechanism";
}
close IN;

open IN,"$outdir/$sample_name\.blast.xls" or die;
open OUT,">$outdir/$sample_name\.anno.xls" or die;
print OUT "Cluster_id\tARO_Accession\tARO_Name\tAMR_Gene_Family\tDrug_Class\tResistance_Mechanism\tTax\n";
while(<IN>){
	chomp;
	my @cut=split /\t/,$_;
	my $length=$cut[5]-$cut[4]+1;
	#next if($cut[-1]<30 or $length<20);
    next if($cut[-1]<20 or $length<20);
	my $seqid=$cut[0];
	my ($gi,$aro,$resistance_type)=(split /\|/,$cut[1])[1,-2,-1];
	my $tax||="none";
	$tax=(split /\]/,(split /\[/,$cut[2])[1])[0];
	if(exists $hash{$aro}){
		print OUT "$seqid\t$aro\t$hash{$aro}\t$tax\n";
	}
	else{next;}
}
close IN;
