from Bio import SeqIO, Seq
import os
import io
import sys
import time
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Assembly.")

    anchor_params = parser.add_argument_group('Anchoring parameters.')
    anchor_params.add_argument('-a', '--assembly',  required=True, help='Path to assembly being anchored.')
    anchor_params.add_argument('-m', '--mapping',  required=True, help='File containing the unique reference mapping coordiantes of each scaffold.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', required = True, help='Output fasta file.')

    args = parser.parse_args()

    mapping_dict = {}
    with open(args.mapping) as f:
        for l in f:
            lvals = l.strip().split()
            scaf = lvals[0]
            map_coord = lvals[2]+':'+lvals[3]+'_'+lvals[4]
            mapping_dict[scaf] = map_coord

    fasta_sequences = SeqIO.parse(open(args.assembly),'fasta')
    out = open(args.output,'w')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name in mapping_dict:
            map_name = mapping_dict[name]
            out.write('>'+map_name+'\n')
            out.write(sequence+'\n')
