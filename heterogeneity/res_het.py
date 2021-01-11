import sys
from pysam import VariantFile

vcf_file =sys.argv[1]
chr_lengths_file = sys.argv[2]
out_file = sys.argv[3]

chrm_dict = {}
tot_chrm_len = 0 
with open(chr_lengths_file,'r') as cf:
	for l in cf:
		chrm, chrm_len = l.strip().split()
		chrm_dict[chrm] = chrm_len
		tot_chrm_len += int(chrm_len)

vcf_in = VariantFile(vcf_file)

sample_list = list((vcf_in.header.samples))

het_counts ={}
hom_counts ={}
for sample in sample_list:
	het_counts[sample] = 0
	hom_counts[sample] = 0

for rec in vcf_in.fetch():
	lvals = str(rec).split()
	chrm = lvals[0]
	if chrm not in chrm_dict or len(lvals[3])>1 or len(lvals[4]) > 1:
		continue
	geno_list = lvals[9:]
	for sample,genotype in zip(sample_list, geno_list):
		geno = genotype.split(':')[0]
		if geno == '0/1':
			het_counts[sample] += 1
		elif geno == '1/1':
			hom_counts[sample] += 1

with open(out_file,'w') as outf:
	for sample in sample_list:
		if sample == "BW_01192":
			continue
		outf.write(sample+'\t'+str(het_counts[sample]/(hom_counts[sample]+het_counts[sample]))+'\n')
