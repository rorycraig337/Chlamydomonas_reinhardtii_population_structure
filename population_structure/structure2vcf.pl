#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long qq(GetOptions);

#script takes input used to create NJ tree and a VCF, and outputs a VCF of only those sites
#usage: perl structure2vcf.pl --structure in.tsv --vcf in.vcf --out out.vcf

my $structure;
my $vcf;
my $out;

GetOptions(
    'structure=s' => \$structure,
    'vcf=s' => \$vcf,
    'out=s' => \$out,
) or die "missing input\n";

open (IN1, "$structure") or die;

my %pos;
my $head = <IN1>;

while (my $line1 = <IN1>) {
        chomp $line1;
        my @cols1 = split(" ", $line1);
        my $site1 = "$cols1[1]\t$cols1[2]";
        $pos{$site1} = 1;
}

close IN1;

open (IN2, "$vcf") or die;
open (OUT, ">$out") or die;

while (my $line2 = <IN2>) {
        chomp $line2;
        unless ($line2 =~ /^#/) {
                my @cols2 = split(/\t/, $line2);
                my $site2 = "$cols2[0]\t$cols2[1]";
                if (exists $pos{$site2}) {
                        print OUT "$line2\n";
                }
        }
        else {
        	print OUT "$line2\n";
    	}
}

close IN2;
close OUT;

exit;
