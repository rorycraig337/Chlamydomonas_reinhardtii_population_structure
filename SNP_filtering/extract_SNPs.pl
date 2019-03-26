#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#script extract SNPs based on presence of ALT allele
#usage: perl extract_SNPs.pl --vcf in.vcf --out out.vcf

my $vcf;
my $out;

GetOptions(
	'vcf=s' => \$vcf,
	'out=s' => \$out,
) or die "missing input\n";

open (IN, "$vcf") or die;
open (OUT, ">$out") or die;

while (my $line = <IN>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @cols = split(/\t/, $line);
		if ($cols[4] ne ".") { #site has an alt allele
			my $alt_count = 0;
			my $good_count = 0;
			my $count = 0;
			foreach my $col (@cols) {
				if ($count > 8) {
					unless ($col eq "null") {
						$good_count++;
						my @calls = split(/:/, $col);
						if ($calls[0] eq "1") {
							$alt_count++;
						}
					}
				}
				$count++;
			}
			unless ($good_count == $alt_count) { #possible that all good calls are alt, and therefore site is invariant
				print OUT "$line\n";	
			}
		}
	}
	else {
		print OUT "$line\n";
	}
}

close IN;
close OUT;

exit;

