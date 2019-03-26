#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qw(GetOptions);

#script produces recombination files for fineSTRUCTURE, taking list of chromosomes and recombination rate in M/bp
#requires .phase files foreach chromosome
#usage: perl phase2rec.pl --chr chr_list.txt --rate 0.00000012 --suffix rec

my $chr;
my $rate;
my $suffix;

GetOptions(
        'chr=s' => \$chr,
        'rate=f' => \$rate,
        'suffix=s' => \$suffix,
) or die "missing input\n";   


open (IN1, "$chr") or die;

while (my $chromosome = <IN1>) {
	chomp $chromosome;
	open (IN2, "$chromosome.phase") or die;
	open (OUT, ">$chromosome.$suffix") or die;
	print OUT "start.pos recom.rate.perbp\n";
	while (my $line = <IN2>) {
		if ($line =~ /^P/) {
			my @snps = split(" ", $line);
			my $snp_count = 0;
			shift @snps; #remove P
			foreach my $snp (@snps) {
				my $next = $snp_count + 1;
				if (exists $snps[$next]) {
					my $bp_dis = $snps[$next] - $snps[$snp_count];
					my $gen_dis = $bp_dis * $rate;
					print OUT "$snps[$snp_count] $gen_dis\n";
				}
				else {
					print OUT "$snps[$snp_count] 0\n";
				}
				$snp_count++;
			}
		}
	}
	close IN2;
	close OUT;
}

close IN1;

exit;
