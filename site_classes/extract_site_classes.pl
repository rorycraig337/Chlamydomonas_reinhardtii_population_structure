#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script to extract site class coordinates from Rob Ness's C. reinhardtii annotation table
#will output bed files for the following site classes CDS,UTR,5UTR,3UTR,introns,intergenic
#each line will be a single base, so merge with bedtools
#intergenic & intronic must be uniquely so
#UTR must not also be CDS, if 5' or 3' must not overlap the other type of UTR
#CDS can overlap other site classes
#optionally outputs annotation table only for CDS, can be used to output 0D, 2D and 4D sites
#usage: perl extract_site_classes.pl --annotation annotation_table.txt --CDS_table

my $annotation;
my $CDS_table = 0;

GetOptions(
        'annotation=s' => \$annotation,
        'CDS_table' => \$CDS_table,
) or die "missing input\n";

open (IN, "$annotation") or die;
open (OUT1, ">intergenic.bed") or die;
open (OUT2, ">introns.bed") or die;
open (OUT3, ">UTR.bed") or die;
open (OUT4, ">5_UTR.bed") or die;
open (OUT5, ">3_UTR.bed") or die;
open (OUT6, ">CDS.bed") or die;

if ($CDS_table == 1) {
        open (OUT7, ">CDS_table.txt") or die;      
}


while (my $line = <IN>) {
        chomp $line;
        my @cols = split(/\t/, $line);
	if ($line !~ /^#/) {
        	if ($cols[18] =~ /\[\]/) {
                	my $start = $cols[1] - 1;
                	print OUT1 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
        	}
        	if ( ($cols[18] =~ /intron/) and ($cols[18] !~ /CDS/) and ($cols[18] !~ /five_prime_UTR/) and ($cols[18] !~ /three_prime_UTR/) and ($cols[18] !~ /rRNA/) ) {
                	my $start = $cols[1] - 1;
                	print OUT2 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
        	}
        	if ( ( ($cols[18] =~ /five_prime_UTR/) or ($cols[18] =~ /three_prime_UTR/) ) and ($cols[18] !~ /CDS/) ) {
                	my $start = $cols[1] - 1;
                	print OUT3 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
        	}
        	if ( ($cols[18] =~ /five_prime_UTR/ ) and ($cols[18] !~ /three_prime_UTR/) and ($cols[18] !~ /CDS/) ) {
                	my $start = $cols[1] - 1;
                	print OUT4 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
        	}
        	if ( ($cols[18] =~ /three_prime_UTR/ ) and ($cols[18] !~ /five_prime_UTR/) and ($cols[18] !~ /CDS/) ) {
                	my $start = $cols[1] - 1;
                	print OUT5 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
        	}
        	if ($cols[18] =~ /CDS/) {
                	my $start = $cols[1] - 1;
	        	print OUT6 "$cols[0]\t$start\t$cols[1]\t$cols[19]\t.\t$cols[21]\n";
        		if ($CDS_table == 1) {
                                print OUT7 "$line\n";  
                        }
		}
	}
	elsif ($CDS_table == 1) {
	       print OUT7 "$line\n";
	}
}

close IN;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;
close OUT6;
if ($CDS_table == 1) {
        close OUT7;
}

exit;
