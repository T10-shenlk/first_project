#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use ScriptManager;
use PerlPackage::gtf_parser::gtf_parser;
use FindBin qw($Bin $Script);
use Data::Dumper;
############################################################
#
#根据ENTRY NCBI ID
#
#
#############################################################
my ($kegg, $eggnog, $gtf, $outdir);
GetOptions(
	"help|h" =>\&USAGE,
	"o:s"=>\$outdir,
	"k:s"=>\$kegg,
	"e:s"=>\$eggnog,
	"g:s"=>\$gtf,
);
&USAGE unless (@ARGV && $kegg && $gtf);

$eggnog ||= $PATH_IN_ENV->{'ensembl_eggnog'};
$outdir ||=".";
mkdir $outdir unless -d $outdir;

my %result;

#read Ensembl Genebank *.dat
#target: ensembl id, Name, Description, GO
foreach my $ensembl (@ARGV){
	my $seqin;
	if($ensembl=~/.gz$/){
		open ENSEMBL, "gunzip -c $ensembl |" or die $!;
		$seqin = Bio::SeqIO->new( 
			-fh => \*ENSEMBL,
			-format => "Genbank"
		);
	}else{
		$seqin = Bio::SeqIO->new( 
			-file => $ensembl,
			-format => "Genbank"
		);
	}

	while(my $seq_obj = $seqin->next_seq()){
		#my $chr = $seq_obj->display_id();
		for my $feat_obj ($seq_obj->get_SeqFeatures) {
			if($feat_obj->primary_tag eq "gene"){
				my $gene_id=($feat_obj->get_tag_values("gene"))[0];
				$gene_id=~s/\.[0-9]+//g;
				$result{$gene_id}{"Gene ID"}=$gene_id;
				#$result{$gene_id}{"Chromosome"}=$chr;
				#$result{$gene_id}{"Start Site"}=$feat_obj->location->start;
				#$result{$gene_id}{"End Site"}=$feat_obj->location->end;
				#$result{$gene_id}{"Direction"}= $feat_obj->location->strand ==1 ? "+" : "-";
				$result{$gene_id}{"Name"}=($feat_obj->get_tag_values("locus_tag"))[0] if $feat_obj->has_tag('locus_tag');
				$result{$gene_id}{"Description"}=($feat_obj->get_tag_values("note"))[0] if $feat_obj->has_tag('note');
			}else{
				if($feat_obj->has_tag('gene')){
					my $gene_id=($feat_obj->get_tag_values("gene"))[0];
					$gene_id=~s/\.[0-9]+//g;
					if($result{$gene_id} && $feat_obj->has_tag('db_xref')){
						for my $val ($feat_obj->get_tag_values('db_xref')) {
							$result{$gene_id}{"GO"} .= $val.";" if $val=~/GO:\d{7}/;
						}					
					}
				}
			}
		}
	}
=pod
	my ($chr, $feature);
	while(<ENS>){
		chomp;
		($chr) = /^LOCUS\s+(\S+)\s+/ if /^LOCUS/;
		if(/^FEATURES/){ #change line breaker if meet FEATURES
			$feature=1;
			$/="\n     gene";
			next;
		}
		if(/\n[A-Z]+\s+|\/\// && $feature){ #change line breaker back after FEATURES
			#print "-------------------------\n".$_."\n";
			$feature=0;
			$/="\n";
			next;
		}

		if($feature){
			#print "-------------------------\n".$_."\n";
			my @lines=split /\n\s{5}\w+\s+/;
			my ($gene_id)= $lines[0] =~ /\s+\/gene=([^"]+)\s*\n/;
			next unless $gene_id;
			$gene_id=~s/\.[0-9]+//g;
			#print $gene_id."\n";
	
			$result{$gene_id}{"Gene ID"}=$gene_id;
			$result{$gene_id}{"Chromosome"}=$chr;
	
			# match complement(14404..29570) or 11869..14409
			($result{$gene_id}{"Direction"}, $result{$gene_id}{"Start Site"}, $result{$gene_id}{"End Site"})= $lines[0]=~/\s*(complement)?\(?(\d+)\.\.(\d+)\)?\s*\n/;
			$result{$gene_id}{"Direction"}= $result{$gene_id}{"Direction"} ? "-" : "+";

			#match /locus_tag="AC215217.1"
			($result{$gene_id}{"Name"})= $lines[0]=~/\/locus_tag="([^"]+)"/;

			#match /note="tubulin beta 8 class VIII [Source:HGNC
			# Symbol;Acc:HGNC:20773]"
			($result{$gene_id}{"Description"})= $lines[0]=~/\/note="([^"]+)"/;
			$result{$gene_id}{"Description"}=~s/\n\s{2,}/ /g if $result{$gene_id}{"Description"};

			shift @lines;
			#match /db_xref="GO:0008203"
			foreach my $line (@lines){
				my @GO;
						$result{$gene_id}{"GO"} .= join(';',@GO).";" if @GO= $line=~/db_xref=".*(GO:\d{7})"/g;
			}
			#print join("\t", map { $result{$gene_id}{$_} || "-" } ("Gene ID", "Chromosome", "Direction", "Start Site", "End Site", "Name", "Description", "GO") ),"\n";
		}else{
			#print "-------------------------\n".$_."\n";
		}
	}
	close ENS;
=cut
}

#read KEGG download file
#target: KEGG, EC, Entrez_geneID, UniProt
local $/="///";
open KEGG, $kegg or die $!;
while(<KEGG>){
	chomp;
	next if /^$/;

	my @lines=split /\n/;
	my ($ens_id) = $_=~ /ENTRY\s+(\w+)/;
	if ($ens_id && $result{$ens_id}) {
		($result{$ens_id}{"KEGG"}) = $_=~/ORTHOLOGY\s+(K\d+)\s+/;
		($result{$ens_id}{"EC"}) = $_=~/ORTHOLOGY.+\[(EC:[0-9\-\.]+)\]/;
		($result{$ens_id}{"Entrez_geneID"}) = $_=~/NCBI-GeneID: (\w+)/;
		($result{$ens_id}{"UniProtAC"}) = $_=~/UniProt: ([0-9A-Za-z ]+)/;
		$result{$ens_id}{"UniProtAC"}=~s/ /;/g if $result{$ens_id}{"UniProtAC"};
	}
}
close KEGG;
local $/="\n";


# read eggNOG - ensembl IDMAP
#target: eggNOG, eggNOG class
my %eggnog;
open EGG, $eggnog or die $!;
while(<EGG>){
	chomp;
	my @c = split /\t/;
	$eggnog{$c[0]}{"eggNOG"} = $c[1];
	$eggnog{$c[0]}{"eggNOG Class"} = $c[2];
}
close EGG;

#read gtf
#target: chromosome, start, end, direction, exons total length
#exons total length will be used to calculate rpkm
open OUT,">$outdir/Annotation.xls";
my @head=("Gene ID", "Chromosome", "Start Site","End Site", "Direction", "Length", "Name", "GO", "KEGG", "EC", "eggNOG Class", "eggNOG", "Description", "Entrez_geneID", "UniProtAC");
print OUT join("\t", @head)."\n";

my %sum;
my $parser = PerlPackage::gtf_parser::gtf_parser->new(file=>$gtf);
while(my $gene=$parser->next_gene()){
	if($gene->filter_feature(attrs =>{gene_biotype=>"protein_coding"})){
		if($result{$gene->{"gene_id"}}){
			my $gene_id=$gene->{"gene_id"};
			$result{$gene_id}{"Chromosome"} = $gene->{"seqname"};
			$result{$gene_id}{"Start Site"} = $gene->{"start"};
			$result{$gene_id}{"End Site"} = $gene->{"end"};
			$result{$gene_id}{"Direction"} = $gene->{"strand"};
			$result{$gene_id}{"Length"} = calculate_exon_cov($gene);

			#eggNOG
			my %protein; 
			my %eggClass;
			foreach my $cds ($gene->descendant(type=>"CDS", hasAttr=> ["protein_id"])){
				if($eggnog{$cds->{"attrs"}{"protein_id"}} && !$protein{$cds->{"attrs"}{"protein_id"}}){
					$result{$gene_id}{"eggNOG"} .= $eggnog{$cds->{"attrs"}{"protein_id"}}{"eggNOG"}.";";
					map { $eggClass{$_} = 1 } (split /;|,/, $eggnog{$cds->{"attrs"}{"protein_id"}}{"eggNOG Class"});
					$protein{$cds->{"attrs"}{"protein_id"}} = 1;
				}
			}
			$result{$gene_id}{"eggNOG"} =~s/;$// if $result{$gene_id}{"eggNOG"};
			$result{$gene_id}{"eggNOG Class"} = join(";", keys %eggClass);

			# print all annotation
			print OUT join("\t", map { $result{$gene_id}{$_} || "-" } @head ),"\n";
			map { $sum{$_}++ if $result{$gene_id}{$_} && $result{$gene_id}{$_} ne "-" } ("Gene ID", "GO", "KEGG", "EC", "eggNOG", "UniProtAC", "Entrez_geneID");
		}
	}
}
close OUT;

#print Annotation summary
open SUM,">$outdir/Annotation.stat";
print SUM "Database\tAnnotated\tPercent\n";
$sum{"Ensembl"}=$sum{"Gene ID"};
map { print SUM join("\t", $_, $sum{$_} || 0, 100*substr($sum{$_}/$sum{"Gene ID"},0,6) ),"\n" } 
	("Ensembl", "GO", "KEGG", "EC", "eggNOG", "UniProtAC", "Entrez_geneID");
close SUM;

sub USAGE {
		my $usage=<<"USAGE";
Usage: perl $0 *.dat.gz -k kegg -g gtf [-o out_dir -e Ensembl_EggNOG.v4.5]
Description: integrate annotation from ensembl and kegg
Options:
	-k <file>	kegg, kegg annotation file, forced
	-g <file>	gtf download from Ensembl, forced
	-e <file>	eggNOG id and ensembl protein id file, optional
	-o <dir>	output directory,optional,default "./"
	-h              help

USAGE
		print $usage;
		exit;
}
