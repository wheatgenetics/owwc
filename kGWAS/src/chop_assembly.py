from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Chop assembly file into files of approximately equal size.")

    chop_params = parser.add_argument_group('Chop parameters.')
    chop_params.add_argument('-a', '--assembly',  required=True, help='Path to reference assembly.')
    chop_params.add_argument('-c', '--chopsize',  required=True, help='Assembly file will be chopped into different files, each containing sequence of around this length.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--outputprefix', required = True, help='Path and prefix of output files.')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.assembly),'fasta')
    
    cum_len = 0 
    count = 0
    output = args.outputprefix+"_0.fa"
    out = open(output,'w')

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        cum_len += len(sequence)
        out.write('>'+name+'\n')
        out.write(sequence+'\n')
        if cum_len > int(args.chopsize):
            out.close()
            count += 1
            cum_len = 0
            output = args.outputprefix+"_"+str(count)+".fa"
            out = open(output,'w')
        