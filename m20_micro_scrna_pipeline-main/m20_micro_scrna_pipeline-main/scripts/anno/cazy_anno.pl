#!/usr/bin/perl -w
#xiaoyue.wang\@oebiotech.com      2017.09.24
use strict;
use Cwd;
use Getopt::Long;
use File::Spec;

my %opts;
GetOptions (\%opts,"aa=s","pre=s","odir=s","e1=s","e2=s","cov=f","cpu=i","help!");
my $usage = <<"USAGE";
Program: $0 v1.0
Discription: CAZy database annotation...
Usage: perl $0 options
    -aa     <file>       input file (protein fasta), required
    -pre    <str>        output prefix, required
    -odir   <dir>        output file dir, default: ./
    -e1     <str>        use evalue < this value when alignment length > 80aa, default:1e-5
    -e2     <str>        use evalue < this value when alignment length <= 80aa, default:1e-3
    -cov    <float_num>  use covered fraction of CAZy Family > this value, default:0.3
    -cpu    <num>        number of processors to use, default:3
    -help                print this help
Example: perl $0 -aa sample.faa -pre sample -odir cazy/

USAGE

die $usage if ( ! $opts{aa} || ! $opts{pre} || $opts{help});
$opts{aa}=File::Spec->rel2abs($opts{aa});
$opts{odir} ||= getcwd;
$opts{odir}=File::Spec->rel2abs($opts{odir});
mkdir $opts{odir} unless( -d $opts{odir});
$opts{e1} ||= "1e-5";
$opts{e2} ||= "1e-3";
$opts{cov} ||= 0.3;
$opts{cpu} ||= 3;
my $out="$opts{odir}/$opts{pre}";

#################################################################################
my $hmmscan="/data/software/hmmer/3.2.1/bin/hmmscan";
my $hmmout = $out . ".dbCAN.hmmscan.out";
my $db_dir = "/data/database/micro/CAZy";
my $db = "$db_dir/dbCAN-HMMdb-V8.txt";
my $R="/home/zjming/anaconda3/bin/R";
#################################################################################

print "\nHmmscan $opts{aa} against dbCAN ...\n\n";
system("$hmmscan --cpu $opts{cpu} -o $hmmout $db $opts{aa}");

print "Parsing annotation result ...\n\n";

open IN,"$db_dir/class_definition.txt" or die "read $db_dir/class_definition.txt: $!\n";
<IN>;
my (@clas,%def);
while(<IN>){
    chomp;
    my @list=split /\t/;
    $def{$list[0]}=$list[1];
    push @clas,$list[0];
}
close IN;

open IN,"$db_dir/all.hmm.ps.len.ps" or die "read $db_dir/all.hmm.ps.len.ps: $!\n";
my %len;
while(<IN>){
    chomp;
    my @list=split /\t/;
    $len{$list[0]}=$list[1];
}
close IN;

open O1,">$out.out" or die "write $out.CAZy.out: $!\n";
open O2,">$out.anno.xls" or die "write $out.CAZy.anno.xls: $!\n";
print O1 "Cluster_id\tFamily\tFamily_len\tFamily_start\tFamily_end\tGene_start\tGene_end\tEvalue\tCovered_fraction\n";
print O2 "Cluster_id\tFamily\tEvalue\tCovered_Fraction\tClass\tClass_Definition\n";
open IN,$hmmout or die "read $hmmout: $!\n";
my (@blo,%evalue,%li1,%li2);
while(<IN>){
    if(/^\/\//){
	my $x=join("",@blo);
	my ($q)=($x=~/^Query:\s+(\S+)/m);
	while($x=~/^>> (\S+.*?\n\n)/msg){
	    my $a=$&;
	    my @c=split(/\n+/,$a);
	    $c[0]=~s/>> //;
	    $c[0]=~s/\s+.*//;
	    for(my $i=3;$i<=$#c;$i++){
		my @d=split(/\s+/,$c[$i]);
		if(!exists $len{$c[0]}){
		    print "Error: not exists the length of $c[0] in file $db_dir/all.hmm.ps.len\n";
		    exit;
		}
		my $covf=($d[8]-$d[7])/$len{$c[0]};
		my $alen=$d[11]-$d[10];
		my $fam=$c[0];
		$fam=~s/\.hmm$//;
		if($fam ne "dockerin" and $fam ne "cohesin" and $fam ne "SLH" and (($alen > 80 and $d[6] < $opts{e1} and $covf > $opts{cov}) or ($alen <= 80 and $d[6] < $opts{e2} and $covf > $opts{cov}))){
		    my $class=(split /\d/,$fam)[0];
		    if((!exists $evalue{$q}) or $d[6] < $evalue{$q}){
			$evalue{$q}=$d[6];
			$li1{$q}="$q\t$fam\t$len{$c[0]}\t$d[7]\t$d[8]\t$d[10]\t$d[11]\t$d[6]\t$covf\n";
			if(!exists $def{$class}){
			    print "Error: not exists the definition of $class in file $db_dir/class_definition.txt\n";
			    exit;
			}
			$li2{$q}="$q\t$fam\t$d[6]\t$covf\t$class\t$def{$class}\n";
		    }
		}else{
		    next;
		}
	    }
	}
	@blo=();
    }else{
	push(@blo,$_);
    }
}
foreach(sort keys %li1){
    print O1 $li1{$_};
}
my (%fc,%fl,%cc,%cl);
foreach(sort keys %li2){
    print O2 $li2{$_};
    chomp($li2{$_});
    my @list=split /\t/,$li2{$_};
    $fc{$list[1]}++;
    $fl{$list[1]}.="$list[0],";
    $cc{$list[4]}++;
    $cl{$list[4]}.="$list[0],";
}
close O1;
close O2;
close IN;

open Fam,"$db_dir/family_information.txt" or die "$!\n";
open O2,"$out.anno.xls" or die "$!\n";
open O3,">$out.anno.1.xls" or die "write $out.CAZy.anno.1.xls: $!\n";
my %hash;
while (my $input=<Fam>){
	chomp $input;
	my @a=split(/\t+/,$input);
	$hash{$a[0]}=$a[-1];
}

while (my $input2=<O2>){
	chomp $input2;
	my @a=split (/\t+/,$input2);
        my $info=$a[1];
        $info=~s/_.*//;
	print O3"$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$hash{$info}\n";
}
close Fam;
close O2;
close O3;
system("rm $out.anno.xls");
system("mv $out.anno.1.xls $out.anno.xls");

open O1,">$out.family.xls" or die "write $out.CAZy.family.xls: $!\n";
open O2,">$out.class.xls" or die "write $out.CAZy.class.xls: $!\n";
print O1 "Family\tGenes_Count\tGenes_List\tClass\tClass_Definition\n";
print O2 "Class\tGenes_Count\tGenes_List\tClass_Definition\n";
foreach(sort keys %fc){
    my $class=(split /\d/,$_)[0];
    $fl{$_}=~s/,$//;
    print O1 "$_\t$fc{$_}\t$fl{$_}\t$class\t$def{$class}\n";
}
foreach(@clas){
    if(exists $cc{$_}){
    	$cl{$_}=~s/,$//;
    	print O2 "$_\t$cc{$_}\t$cl{$_}\t$def{$_}\n";
    }
}
close O1;
close O2;

system("cut -f1,2,4 $out.class.xls > $out.class.txt");
open R,">$out.cmd.r" || die $!;
print R"library(ggplot2)\ncount <- read.table(file=\"$out.class.txt\", header=T, quote=\"\", sep=\"\\t\", check.names=F)
p=ggplot(data=count, aes(x=Class_Definition, y=Genes_Count,fill=Class_Definition))+geom_bar(stat=\"identity\", width = 0.5, position = position_dodge(0.7))+
guides(fill=FALSE)+
theme(panel.background=element_rect(fill='transparent',color =\"gray\"),  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,color = \"black\",size=9))+
theme(panel.grid =element_blank())+
geom_text(mapping = aes(label = count\$Genes_Count),size=3.5,vjust=-0.2)+
xlab(\"Class definition\") + ylab(\"Genes number\")
ggsave(\"$out.class.pdf\",width=8,height=7,plot=p)
ggsave(\"$out.class.png\",type=\"cairo-png\",width=8,height=7,plot=p)\n";
close R;
#system("export LD_LIBRARY_PATH=/home/fanyucai/lib/usr/lib64/:\$LD_LIBRARY_PATH");
system("$R --restore --no-save < $out.cmd.r");
system("rm -rf $out.cmd.r $out.class.txt");

print "Finished\n\n";
