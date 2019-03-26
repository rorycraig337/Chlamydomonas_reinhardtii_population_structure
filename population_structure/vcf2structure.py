from itertools import combinations
import sys, vcf, random

def vcf2structure(vcf_file, num_snps=1000000000000, min_fraction_called = 1.0, min_distance=20000, depth_cutoff=0, gq_cutoff=0, ploidy=1, max_heterozygotes=0):
	#vcf2structure.py all_wt.nyu.v5.3.final.vcf
	#if I allow missing data I will need to mod the script to deal with that.
		#insert -1 for missing alleles
	vcf_object = vcf.Reader(open(vcf_file, 'r'))
	#Generate hash of SNPs
	snp_data = {}
	for sample in vcf_object.samples:
		snp_data[sample] = []
	n=0
	loci = [['', 0]]
	for record in vcf_object:
		if record.POS > loci[-1][1]+min_distance or record.CHROM != loci[-1][0]:
			if record.num_called >= min_fraction_called * len(record.samples) and \
				record.is_snp and \
				len(record.get_hets()) <= max_heterozygotes:
				gts = [i['GT'] for i in record.samples]
				dps = [i['DP'] for i in record.samples]
				gqs = [i['GQ'] for i in record.samples]
				if min([d for d in dps if type(d)==int]) > depth_cutoff and \
					min([q for q in gqs if type(q)==int]) > gq_cutoff:
					for s in range(len(vcf_object.samples)):
						snp_data[vcf_object.samples[s]].append(gts[s])
					loci.append([record.CHROM, record.POS])
					n+=1
					print n, record.CHROM, record.POS, gts
					if n >num_snps:
						break
	#write outputfile
	print snp_data
	o = open(vcf_file[:-3] + 'structure.inp', 'w')
	o.write( "\t".join([l[0] + "_" + str(l[1]) for l in loci]) + "\n")
	for sample in vcf_object.samples:
		o.write("%s\t%s\n" %(sample, "\t".join(snp_data[sample])))
	o.close()

def runStructure(K,infile, outfile, structure_path = 'structure'):
	#this is a naive implementation that assumes you already have extraparams and mainparams set
	cmd_line = " %s -K %i -i %s -o %s 1>/dev/null" %(structure_path, K, infile, outfile)
	print cmd_line
	child = subprocess.Popen(cmd_line, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	child.wait()
	results = parse_structure(outfile+"_f")
	return results
	
def parse_structure(outfile):
	results = {}
	f = open(outfile, 'r')
	line = f.readline()
	while line.strip() != 'Inferred ancestry of individuals:':
		line = f.readline()
	line = f.readline()
	line = f.readline()
	while len(line) >1:
		label = line.split()[1]
		results[label] = [float(i) for i in line.split(':')[-1].split()]
		print line.strip()
		line = f.readline()
	return results
	
### MAIN ###
"""
run this script from bash commandline like this
python vcf2structure.py my_vcf.vcf.gz

vcf should be bgzipped and tabix indexed
Script will create an output that can be fed into STRUCTURE

REquirements:
You need pyVCF installed 

Technically you can use the function runStructure to sctually run the program but that is not used here
"""


vcf_file = sys.argv[1] # this reads the first argument from commandline as the VCF file name
num_snps=10000000000000 # the script will stop after this many SNPs (set to infinity here)
min_fraction_called = 1 # You need this fracttion of individuals to be called to bother with the site
min_distance=20000 # There is no sense in inputting SNPs that are too close because they are linked. This determines minimum distance between SNPs for structure
depth_cutoff=0 # Depth of reads per individual
gq_cutoff=0 # quality of SNP call per individual


vcf2structure(vcf_file,  num_snps=num_snps, min_fraction_called = min_fraction_called, min_distance=min_distance, depth_cutoff=depth_cutoff, gq_cutoff=gq_cutoff)
