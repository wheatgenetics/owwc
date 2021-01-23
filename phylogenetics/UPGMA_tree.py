"""
 make UPGMA tree 

"""
import sys
import random
random.seed(123)
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.Consensus import *

#--------------- Convert to fasta -------------------------------
def format_conv(csv_file, fasta_file):
	with open(csv_file, 'r') as inp, open(fasta_file,'w') as fas:
		inp.readline()
		for line in inp:
			A=line.strip().split('\t')
			acc = A[0]
			seq = ''.join(A[1:])
			fas.write('>'+ acc + '\n')
			fas.write(     seq + '\n')


#-------------- patching section for iTOL visualisation------------------------------
## Following are the parts from Biopython Consensus code 
from Bio.Phylo.NewickIO import Writer

def _get_comment(clade):
    if hasattr(clade, 'comment') and clade.comment:
        return _format_comment(str(clade.comment))
    else:
        return ''

def _info_factoryJ(self, plain, confidence_as_branch_length,
                      branch_length_only, max_confidence, format_confidence,
                      format_branch_length):
        """Return a function that creates a nicely formatted node tag."""
        if plain:
            # Plain tree only. That's easy.
            def make_info_string(clade, terminal=False):
                return _get_comment(clade)

        elif confidence_as_branch_length:
            # Support as branchlengths (eg. PAUP), ignore actual branchlengths
            def make_info_string(clade, terminal=False):
                if terminal:
                    # terminal branches have 100% support
                    return (':' + format_confidence % max_confidence) + _get_comment(clade)
                else:
                    return (':' + format_confidence % clade.confidence) + _get_comment(clade)

        elif branch_length_only:
            # write only branchlengths, ignore support
            def make_info_string(clade, terminal=False):
                return (':' + format_branch_length % clade.branch_length) + _get_comment(clade)

        else:
            # write support and branchlengths (e.g. .con tree of mrbayes)
            def make_info_string(clade, terminal=False):
                if (terminal or
                        not hasattr(clade, 'confidence') or
                        clade.confidence is None):
                    return (':' + format_branch_length
                            ) % (clade.branch_length or 0.0) + _get_comment(clade)
                else:
                    return (':' + format_branch_length + '[' + format_confidence + ']'   ## this  modify
                            ) % (clade.branch_length or 0.0, clade.confidence ) + _get_comment(clade) ## and this modify

        return make_info_string

#---------------------------------------
def bootstrap(msa, times=10):
    """Generate bootstrap replicates from a multiple sequence alignment object

    :Parameters:
        msa : MultipleSeqAlignment
            multiple sequence alignment to generate replicates.
        times : int
            number of bootstrap times.
    """

    length = len(msa[0]) 
    i = 0
    while i < times:
        i += 1
        item = None
        for j in range(length):
            col = random.randint(0, length - 1)
            if not item:
                item = msa[:, col:col + 1]
            else:
                item += msa[:, col:col + 1]
        yield item

def bootstrap_trees(msa, times, tree_constructor=None, distance_calculator = None):
    """Generate bootstrap replicate trees from a multiple sequence alignment.

    :Parameters:
        msa : MultipleSeqAlignment
            multiple sequence alignment to generate replicates.
        times : int
            number of bootstrap times.
        tree_constructor : TreeConstructor
            tree constructor to be used to build trees.
    """

    msas = bootstrap(msa, times)
    count = 0
    for aln in msas:
        count += 1
        print(count)
        dmat = distance_calculator.get_distance(aln)
        tree = tree_constructor.upgma(dmat)
        yield tree


class KmerDistanceCalculator:

    def __init__(self,skip_letters=None):
        if skip_letters:
            self.skip_letters = skip_letters


    def _pairwise(self, seq1, seq2):
        score = 0
        max_score = 0
        for l1, l2 in zip(seq1, seq2):
            if l1  in self.skip_letters and l2  in self.skip_letters:
                continue
            score += int(l1 == l2)
            max_score += 1

        if max_score == 0:
            return 1  # max possible scaled distance
        return 1 - (score * 1.0 / max_score)

    def get_distance(self, msa):
        """Return a DistanceMatrix for MSA object.
        :Parameters:
            msa : MultipleSeqAlignment
        """
        if not isinstance(msa, MultipleSeqAlignment):
            raise TypeError("Must provide a MultipleSeqAlignment object.")

        names = [s.id for s in msa]
        dm = DistanceMatrix(names)
        for seq1, seq2 in itertools.combinations(msa, 2):
             pairwise_dist = self._pairwise(seq1, seq2)
             dm[seq1.id, seq2.id] = pairwise_dist
        return dm
        
## Until here from Biopython Consensus code
#-------------------------------------------------



if __name__ == "__main__":
   
  ## input
  input_file = sys.argv[1]
  output_original_tree = sys.argv[2]
  output_bootstrapped_tree = sys.argv[3]
  bootstrap_times = int(sys.argv[4])

  format_conv(csv_file = sys.argv[1], fasta_file = 'kmers_phylo.fa')

  aln = AlignIO.read(open('kmers_phylo.fa'), 'fasta')
  dcalc = KmerDistanceCalculator('0')
  dmat = dcalc.get_distance(aln)
  dtc = DistanceTreeConstructor(method = "upgma")

  orig_tree = dtc.upgma(dmat) ## original Tree

  Phylo.write(orig_tree, output_original_tree, "newick") #

  bootrees = bootstrap_trees(aln, times = bootstrap_times, tree_constructor = dtc, distance_calculator = dcalc)

  trees = list(bootrees)

  support_tree = get_support(orig_tree, trees)

  Writer._info_factory = _info_factoryJ

  Phylo.write(support_tree, output_bootstrapped_tree, "newick")
