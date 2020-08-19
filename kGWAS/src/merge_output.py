import sys, os

output_dir = sys.argv[1]
merged_output = sys.argv[2]

h = {}
hc = {}
for filename in os.listdir(output_dir):	
	if filename.startswith('.'):
		continue

	with open(output_dir + '/' +filename, 'r') as f:
		for l in f:
			lvals = l.strip().split()
			coord = lvals[0]+':'+lvals[1]+'_'+lvals[2]
			pval = float(lvals[3])
			cor = float(lvals[4])
			if cor < 0 :
				pval = -1*pval
			n_kmers = int(lvals[5])
			if coord in hc:
				hc[coord] += n_kmers
				if pval in h[coord]:
					h[coord][pval] = [cor,h[coord][pval][1]+n_kmers]
				else:
					h[coord][pval] = [cor,n_kmers]
			else:
				hc[coord] = n_kmers
				h[coord] = {}
				h[coord][pval] = [cor,n_kmers]

with open(merged_output, 'w') as out:
	for coord in hc:
		if hc[coord] > 50:
			chrm = coord.split(':')[0]
			start,end = coord.split(':')[1].split('_')
			for pval in h[coord]:
				cor = h[coord][pval][0]
				n_kmers = h[coord][pval][1]
				out.write(chrm + "\t" + start + "\t" + end + "\t" + str(abs(pval)) + "\t" + str(cor) + "\t" + str(n_kmers) + "\n")
