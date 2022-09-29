#!/usr/bin/perl -w

use strict;
no warnings 'uninitialized';

##########################################################
## Written by KS Geist 									##
## 12 July 2021											##
##														##
## Counts the number of times a gene is flagged as sig	##
## by DESeq2 iteratively in order to calculate how 		##
## often resampling yields the same genes.				##
##														##
## USAGE: Call me from the directory with the DESeq		##
## files for each resample. Then, in the terminal run:	##
##		perl deseq_resample_counter.pl <species name>	##
##														##
##	NOTE: Assumes resampling was k = 1000 as used in	##
## Favreau, Geist et al. (2021). 		. 				##
##														##
##########################################################

if (@ARGV != 1) {
	die "\nUsage: deseq_resample_counter.pl <species name>\n";
}
my $species = shift(@ARGV);

## Set the path as the current working directory
my $path = `pwd`;
chomp $path; 		# strip terminal white space from the path

## Print the outfile as:
my $out = $species."_deseq_gene_identities_counted.txt";

## Changes directory to the current wd
chdir($path) or die "$!";
## Opens directory
opendir (DIR, $path) or die "$!";

## Creates an array of all files in the directory
my @files = readdir DIR; 
closedir DIR;

my $sample = "";
my @line = ();
my %genes = ();			## Hash with gene ID as key and the number of times sign. by DESeq2
my $k = 1000;			## 1000 resampling events

foreach my $file (@files) {
	
	if ($file =~ /.*_Resampled_DESeq_Counts_(\d+)/) {
		$sample = $1;
		print "Processing...".$sample."\n";

 		open(IN, $file) or die("Could not open $file.\n\n");

			while (my $line = <IN>) {
# 				print $line;

				## If it isn't the header line
				# 			[tab]	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
				if ($line !~ /^\t/) {
				
					## Store the line into an array by splitting on white space
					@line = split(/\s+/, $line);
					## Then check to see if the gene was significant:
					if ($line[6] ne 'NA') {		## Have to do this because there are some NAs we heed to avoid
						if ($line[6] < 0.05) {
							$genes{$line[0]}++;		## Tally the gene id's as seen 
						}
					}			
				}
			}
		close(IN);

 	}
}

open(OUT, "> $out") or die "Cannot write outfile: $!";
## Print the header:
print OUT "geneID\tnum_times_seen\tprop_seen\n";

foreach my $g (sort keys %genes) {
	my $prop = $genes{$g}/$k;
    print OUT "$g\t$genes{$g}\t$prop\n";
}

