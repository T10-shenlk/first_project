#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################
#
#######################################################################################
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
# GetOptions
my ($fIn,$dOut,$m);

GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$dOut,
				"i:s"=>\$fIn,
				"m:s"=>\$m,
				) or &USAGE;
&USAGE unless ($fIn and $dOut);
# ------------------------------------------------------------------
# 
# ------------------------------------------------------------------
$m||=3;
mkdir $dOut unless(-d $dOut);
open (IN,$fIn) or die $!;
my $out_prefix=basename($fIn);
$out_prefix=~s/.xml$//;
$/="<Iteration>\n";

my @in=<IN>;
close IN;

my $file=join("",@in);
my ($header,$blastrs,$ender)=$file=~/^(.+?)(<Iteration>.+<\/Iteration>)(.+)/sm;

my @blast=$blastrs=~/(<Iteration>.+?<\/Iteration>)/smg;


#my $header=shift @blast;
#my $ender=pop @blast;

my $num=int(scalar(@blast) / $m ) + 1;
print $num,"\n";
for (my $i=0;$i < $m  ;$i++) {
	open OUT,">$dOut/$out_prefix.$i.xml";
	my $end=($i+1)*$num-1;
	$end > $#blast ? $end=$#blast : $end=$end;
	print OUT join("",$header,@blast[$i*$num..$end],$ender);
}


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
################################################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

################################################################################################################


sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

################################################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

################################################################################################################
sub USAGE {
	my $usage=<<"USAGE";
 ProgramName:
     Version:	$version
Program Date:	2012.07.02
 Description:	
       Usage:
		Options:
		-i <file>	input file,xxx format,forced

		-o <file>	output file,optional

		-h		help

USAGE
	print $usage;
	exit;
}
