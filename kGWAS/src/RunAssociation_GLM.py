
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@authors: kumar gaurav

"""

from KmerProjection_GLM import KmerProjection
from Phenotype_GLM import Phenotype
import sys 
import argparse
import os.path
from pathlib import Path

if __name__ == '__main__':
    
    if sys.version_info[0] < 3:
        sys.stderr.write('Python version 3 or above required')
        
    parser = argparse.ArgumentParser(description = "AgRenSeq with Generalized Linear Model")
    
    inputMatrix = parser.add_argument_group('Kmer presence/absence matrix')
    inputMatrix.add_argument('-i', '--inputmatrix',  required=True, help='Path to file containing the gzipped presence/absence matrix of kmers')
    inputMatrix.add_argument('-hd', '--header',  required=True, help='Path to file containing the header for the presence/absence matrix of kmers')
       
    assembly = parser.add_argument_group('Assembly for kmer projection')
    assembly.add_argument('-a','--assembly', required=True, help='Path to assembly file of accession onto which kmers are mapped for plotting')

    phenotype = parser.add_argument_group('Phenotype')
    phenotype.add_argument('-p','--phenotype', required=True, help='Path to phenotype file')

    SNPmarker = parser.add_argument_group('SNP markers matrix')
    SNPmarker.add_argument('-s','--snp', required=True, help='Path to file containing matrix of SNP markers to compute PCA and correct for population structure')
    
    SNPmarker.add_argument("-dim", "--pcadimensions", type=int, default = 3, help="Number of significant PCA dimensions retained for regression analysis")

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', required=True, help='Path to output file to store the association values corresponding to contigs in the given assembly')
    
    parser.add_argument("-u", "--usable", help="file containing the list of usable accessions")

    parser.add_argument("-c", "--correlationthreshold", type=float, default = 0.2, help="only those k-mers whose correlation with phenotype is greater than this value are retained for  regression analysis")

    parser.add_argument("-pv", "--pvalthreshold", type=float, default = 15.0, help="only those k-mers with p-value greater than this value are retained")

    parser.add_argument("-st", "--stackman", help="convert phenotype scores from Stackman's IT to AgRenSeq scores", action="store_true")

    parser.add_argument("-per", "--permute", help="permute phenotype", action="store_true")
    
    parser.add_argument("-sub", "--subsample", type=int, help="random subsample of this size")
    
    
    args = parser.parse_args()
    
    inputMatrix_filename = args.inputmatrix
    try:
        if os.path.getsize(inputMatrix_filename) == 0:
            print("\nFile " + inputMatrix_filename + " is empty.")
            sys.exit()
    except OSError:
        print("\nFile " + inputMatrix_filename + " does not exist or is not accessible.")
        sys.exit()


    header_filename = args.header
    try:
        if os.path.getsize(header_filename) == 0:
            print("\nFile " + header_filename + " is empty.")
            sys.exit()
    except OSError:
        print("\nFile " + header_filename + " does not exist or is not accessible.")
        sys.exit()


    assembly_filename = args.assembly
    try:
        if os.path.getsize(assembly_filename) == 0:
            print("\nFile " + assembly_filename + " is empty.")
            sys.exit()
    except OSError:
        print("\nFile " + assembly_filename + " does not exist or is not accessible.")
        sys.exit()
        
        
    
    phenotype_filename = args.phenotype
    
    try:
        p = Phenotype(phenotype_filename, args.stackman)
    except FileNotFoundError:
        print("\nFile " + phenotype_filename + " does not exist.")
        sys.exit()

    if args.usable is not None:
        try:
            p.selectAccessions(args.usable)
        except FileNotFoundError:
            print("\nFile " + args.usable + " does not exist.")
            sys.exit()
    
    if args.subsample is not None:        
        p.selectRandomAccessions(args.subsample)
    
    if args.permute:
        p.permutePhenotype()
    
    print("\nThe number of accessions used for association are: " + str(len(p.phenoScores_dict)))


    snp_filename = args.snp
    try:
        if os.path.getsize(snp_filename) == 0:
            print("\nFile " + snp_filename + " is empty.")
            sys.exit()
    except OSError:
        print("\nFile " + snp_filename + " does not exist or is not accessible.")
        sys.exit()
        
    pca_dimensions = args.pcadimensions

    cor_threshold = args.correlationthreshold
    pval_threshold = args.pvalthreshold


    output_filename = args.output
    if os.path.isfile(output_filename):
        print('\nFile ' + output_filename + ' already exists. \nTerminate the program if you do not want to overwrite.')
             
				
    projection = KmerProjection( p, assembly_filename, inputMatrix_filename, header_filename, snp_filename, pca_dimensions, cor_threshold, pval_threshold )
			
    projection.readAssembly_GLM()
    projection.readMatrix_GLM()
    projection.writeAssociationScore_GLM( output_filename )