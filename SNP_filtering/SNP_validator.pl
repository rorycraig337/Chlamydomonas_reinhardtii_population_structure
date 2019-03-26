#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#validates and corrects SNPs after indel processing, which can leave SNPs in VCF even though there are no isolates with the alt allele
#such SNPs will be corrected to invariant sites
#perl SNP_validator.pl --vcf in.vcf --out out.vcf

my $vcf;
my $out;

GetOptions(
    'vcf=s' => \$vcf,
    'out=s' => \$out,
) or die "missing input\n";

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

my $change = 0;
my $check_count = 0;

while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#CHROM/) {
	my @head_cols = split(/\t/, $line);
	$check_count = (scalar @head_cols);
    }
    unless ($line =~ /^#/) {
    	my @cols = split(/\t/, $line);
        if ($cols[4] ne ".") { #site is/was a SNP
            my $alt_count = 0;
            my $count = 0;
			foreach my $col (@cols) {
				if ($count > 8) {
					my @calls = split(/:/, $col);
					if ($calls[0] eq "1") {
						$alt_count++;
					}
					if ( ($calls[0] ne "null") and ($calls[0] ne ".") and ($calls[0] > 1) ) {
						die "There should be no multi-allelic sites: $line\n";
					}
				}
				$count++;
			}
			if ($count ne $check_count) {
				print "ERROR: line has $count columns, expected $check_count\n";
				die "$line\n";
			}
			if ($alt_count > 0) { #site is a SNP
				print OUT "$line\n";
			}
			elsif ($alt_count == 0) { #site is invariant
				$cols[4] = ".";
				print OUT join("\t", @cols), "\n";
				$change++;
			}
        }
        else {
        	print OUT "$line\n";
        }
    }
    else {
    	print OUT "$line\n";
    }
}

close IN;
close OUT;

print "changed $change SNPs to invariant\n";

exit;
