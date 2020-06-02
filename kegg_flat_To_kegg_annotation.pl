#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use ScriptManager;
use PerlPackage::gtf_parser::gtf_parser;
use FindBin qw($Bin $Script);
use Data::Dumper;

my ($kegg,$outdir);
GetOptions(
        "help|h" =>\&USAGE,
        "o:s"=>\$outdir,
        "k:s"=>\$kegg,
);

#&USAGE unless (@ARGV && $kegg );
$outdir ||=".";
mkdir $outdir unless -d $outdir;

my %result;
local $/="///";
open KEGG, $kegg or die $!;
while(<KEGG>){
        chomp;
        next if /^$/;
	
        my @lines=split /\n/;
        my ($ens_id) = $_=~ /ENTRY\s+(\w+)/;
#        if ($ens_id || $result{$ens_id}) {
		($result{$ens_id}{"Gene ID"}) = $ens_id;
                ($result{$ens_id}{"KEGG"}) = $_=~/ORTHOLOGY\s+(K\d+)\s+/;
                ($result{$ens_id}{"EC"}) = $_=~/ORTHOLOGY.+\[(EC:[0-9\-\.]+)\]/;
                ($result{$ens_id}{"Entrez_geneID"}) = $_=~/NCBI-GeneID: (\w+)/;
                ($result{$ens_id}{"UniProtAC"}) = $_=~/UniProt: ([0-9A-Za-z ]+)/;
                $result{$ens_id}{"UniProtAC"}=~s/ /;/g if $result{$ens_id}{"UniProtAC"};
#        }
}
close KEGG;
local $/="\n";
open KEGG_Result, ">$outdir/kegg.txt" or die $!;
my @head=("Gene ID", "KEGG","EC","Entrez_geneID","UniProtAC","UniProtAC");
print KEGG_Result join ("\t",@head) ."\n";
foreach my $ens_id (keys %result){
	print KEGG_Result join("\t", map { $result{$ens_id}{$_} || "-" } @head ),"\n";
}
close KEGG_Result;

############################FUNCTION#########################################
sub USAGE {
                my $usage=<<"USAGE";
Usage: perl $0 -k kegg  [-o out_dir]
Description: integrate annotation from ensembl and kegg
Options:
        -k <file>       kegg, kegg annotation file, forced
        -o <dir>        output directory,optional,default "./"
        -h              help

USAGE
                print $usage;
                exit;
}
