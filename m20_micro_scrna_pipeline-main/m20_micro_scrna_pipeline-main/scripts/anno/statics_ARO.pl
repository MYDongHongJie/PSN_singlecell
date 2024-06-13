#!/usr/bin/perl -w
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use FindBin qw($Bin);
my $R="/home/sunkun/Miniconda3/envs/anvio-7.1/bin/R";
my($pep,$outdir,$result);
GetOptions(
    "i:s"=>\$pep,
    "o:s"=>\$outdir,
    "r:s"=>\$result,
);

sub usage{
        print qq{
This script was used to annotated card.
notice:please put the ncbi++ in your enviroment.
usage:
perl $0 -p protein.fasta -i /xx/final.xls -o /xx/result.txt
-i          input the final result of the blast card.
-o          output the result.
-r	    output picture name.

};
    exit;
}
if(!$pep||!$outdir||!$result){
    &usage();
}

open IN,"$pep" or die;
my %uniq=();
my %hash=();
<IN>;
while(<IN>){
	chomp;
	my @cut=split /\t/,$_;
	if(exists $uniq{$cut[0]}){
		if($uniq{$cut[0]} eq $_){next;}
		else{
			if(exists $hash{$cut[2]}){$hash{$cut[2]}++;}
			else{$hash{$cut[2]}=1;}
		}
	}
	else{
		if(exists $hash{$cut[2]}){$hash{$cut[2]}++;}
                else{$hash{$cut[2]}=1;}
		$uniq{$cut[0]}=$_;
	}
}
close IN;

my $i=1;
open OUT,">$outdir/$result\.statics.txt" or die;
print OUT "AROName\tNum\tpercent\n";
my $temp=0;my $other;
foreach my $keys(sort {$hash{$b} <=> $hash{$a} }keys %hash){
	if($i>10){
		$other+=$hash{$keys};
	}
	$temp=$temp+$hash{$keys};
	$i++;
}

$i=1;
foreach my $keys(sort {$hash{$b} <=> $hash{$a} }keys %hash){
	if($i<=10){
		my $percent=$hash{$keys}/$temp*100;
        my @Name=split /\s+/, $keys;
        my $AROName='';
        if(@Name>3){
            $AROName=$Name[0]."_".$Name[1]."_".$Name[2];
        }else{
            $AROName=join("_",@Name);
        }
		print OUT "$AROName\t$hash{$keys}\t";
		print OUT sprintf("%.2f", $percent)."%\n";
	}else{
		last;
	}
	$i++;
}

$other=0 if(!defined $other || $other eq "");
my $other_per=sprintf("%.2f",$other*100/$temp);
print OUT "other\t$other\t$other_per%\n";
close OUT;

open CMD,">$outdir/$result\_cmd.r" or die;
print CMD "
library(RColorBrewer)
pdf(\"$outdir/$result.pdf\",height=6,width=10)
par(mfrow=c(1,2))
col=c(\"lightblue\",\"wheat\",\"darkgreen\",\"limegreen\",\"lightseagreen\",\"lightsteelblue\",\"rosybrown\",\"tomato\",\"lightpink\",\"deeppink\",\"gray\",\"gold\",\"slateblue\",\"mediumvioletred\",\"deepskyblue\",\"darkorange\",\"midnightblue\",\"forestgreen\",\"violet\",\"darkgoldenrod\")
a<-read.table(\"$outdir/$result\.statics.txt\",header=T,sep=\"\t\")
legend<-paste(a\$AROName,a\$Num,sep=\": \")
labels<-a\$percent
y=a\$Num
pie(y,labels=labels,col=col,cex=0.6,edges = 300,radius = 1,border =\"white\",main=\"Functional Categories\")
plot.new()
legend(\"topleft\",legend=legend,col=col,fill=col,box.lwd=0,cex=1,border =\"white\")
dev.off()
\n";

system("$R --restore --no-save < $outdir/$result\_cmd.r");
system "rm $outdir/$result\_cmd.r";

close CMD;
