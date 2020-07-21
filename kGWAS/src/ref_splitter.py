from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Split reference assembly.")

    split_params = parser.add_argument_group('Split parameters.')
    split_params.add_argument('-a', '--assembly',  required=True, help='Path to reference assembly.')
    split_params.add_argument('-s', '--splitsize',  required=True, help='Splits of each scaffold.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', required = True, help='Output fasta file.')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.assembly),'fasta')
    out = open(args.output,'w')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        len_seq = len(sequence)
        chunk_size = int(args.splitsize)
        n_splits = int(len_seq/chunk_size)+1

        for i in range(n_splits):
            if i == n_splits-1:
                sequence_chunk = sequence[i*chunk_size:len_seq]
                index = str(i*chunk_size)+"_"+str(len_seq)
            else:
                sequence_chunk = sequence[i*chunk_size:(i+1)*chunk_size]
                index = str(i*chunk_size)+"_"+str((i+1)*chunk_size)

            chunk_name = name+':'+index
            out.write('>'+chunk_name+'\n')
            out.write(sequence_chunk+'\n')

            