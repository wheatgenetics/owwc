import sys 
import gzip
import io
import argparse
import time
from bitarray import bitarray 

if __name__ == '__main__':
        
    parser = argparse.ArgumentParser(description = "Number of valid kmers for Bonferroni correction")
    
    inputMatrix = parser.add_argument_group('Kmer presence/absence matrix')
    inputMatrix.add_argument('-i', '--inputmatrix',  required=True, help='Path to file containing the gzipped presence/absence matrix of kmers')
    inputMatrix.add_argument('-hd', '--header',  required=True, help='Path to file containing the header for the presence/absence matrix of kmers')
    
    parser.add_argument("-u", "--usable", required=True, help="file containing the list of usable accessions for association mapping")

    parser.add_argument("-mc", "--mincount", type=int, default = 4, help="only those k-mers which are present/absent in more than this number of accessions are used for Bonferroni correction.")

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', required=True, help='Path to output file to store the number of valid kmers in the matrix for Bonferroni correction')
    
    args = parser.parse_args()
    
    gz = gzip.open(args.inputmatrix, 'r')
    f = io.BufferedReader(gz)

    with open(args.header, 'r') as h:
        header = [l.strip() for l in h]
   
    accessions_dict = {}
    for i in range(len(header)):
        accessions_dict[header[i]] = i

    with open(args.usable, 'r') as u:
        usable = [l.strip() for l in u]

    accessions = list(set(header) & set(usable))
    n_accessions = len(accessions)

    progress = 0
    t_init = time.time()
    k_count = 0

    for inputline in f:
        # track progress
        progress += 1
        if progress%1000000 == 0:
            print (str(progress)+" k-mers parsed in time: " +str(time.time() - t_init))

        split = inputline.decode("utf-8").split()
        presence_bitarray = bitarray(split[1])
        
        presence_sum = 0
        for accession in accessions:
            index = accessions_dict[accession]
            if presence_bitarray[index]:
                presence_sum += 1

        if presence_sum < args.mincount or presence_sum > n_accessions-args.mincount:
            continue

        k_count += 1
                                    
    gz.close()

    with open(args.output, 'w') as out:
        out.write(str(k_count)+" out of "+ str(progress)+" k-mers are usable for Bonferroni correction")
    