#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import random


"""

This data structure holds phenotype scores for use in correlation pre-filtering and regression analysis.

@authors: kumar gaurav

"""

class Phenotype(object):

    # dictionary mapping phenotype scores in Stackman's IT to AgRenSeq scores
    stackman_to_agrenseq_dict =  {'3+': '-1.33', '1': '1', '0': '2', '3': '-1', '2': '0', '4': '-2', '3-': '-0.67', '2-': '0.33', ';' : '1.67', '1-': '1.33', '1+': '0.67', '2+': '-0.33', '0;' : '1.67', '1;' : '1'}

    """
    @param phenotypeFile: name of the file containing the phenotype scores.
    @param stackman_conversion: 
                        If False(default value), the input phenotype scores are assumed to be numeric values, with a higher score representing a higher resistance. 
                        If True, the input phenotype scores are assumed to be in the Stackman's IT, which are then converted to AgRenSeq scores as specified at https://github.com/steuernb/AgRenSeq.
    """    
    def __init__(self, phenotypeFile, stackman_conversion = False):

        self.phenoScores_dict = {}

        # store the phenotype scores in the form of a dictionary as well as pandas series
        self.phenoScores_series = self.readScores(phenotypeFile, stackman_conversion)


    def readScores(self, phenotypeFile, stackman_conversion):

        with open(phenotypeFile) as f:
            for l in f:
                lvals = l.split()
                if stackman_conversion:
                    try: 
                        average_phenotype = np.mean(np.array([self.stackman_to_agrenseq_dict[pheno] for pheno in lvals[1:]]).astype(np.float)) 
                    except KeyError:
                        print("Phenotype scores are not in Stackman's IT.")
                        sys.exit()
                else:
                    try:
                        average_phenotype = np.mean(np.array(lvals[1:]).astype(np.float))
                    except KeyError:
                        print("Phenotype scores are not in numeric format. Please use the flag --stackman if the scores are in Stackman's IT.")
                        sys.exit()

                assert not np.isnan(average_phenotype)
                self.phenoScores_dict[lvals[0]] = average_phenotype

        return pd.Series(self.phenoScores_dict)

                     
    """
    Only the phenotype scores of those accessions in the usable file, which also exist in the original phenotype file, are retained.
    @param usableFile
    """           
    def selectAccessions(self, usableFile):

        with open(usableFile) as f:
            usable_accessions = f.readlines()
        
        self.phenoScores_dict = dict((accession.strip(), self.phenoScores_dict[accession.strip()]) for accession in usable_accessions if accession.strip() in self.phenoScores_dict)
        self.phenoScores_series = pd.Series(self.phenoScores_dict)
        
        
    def permutePhenotype(self):

        permutation = np.random.permutation(self.phenoScores_series.values)
        self.phenoScores_series = pd.Series(permutation, index = self.phenoScores_series.index)
        self.phenoScores_dict = self.phenoScores_series.to_dict()
        

    """
    A random subsample of the phenotype scores is retained.
    @param n_accessions: size of subsample
    """
    def selectRandomAccessions(self, n_accessions):

        self.phenoScores_dict = dict(random.sample(self.phenoScores_dict.items(), n_accessions))
        self.phenoScores_series = pd.Series(self.phenoScores_dict)
        
        