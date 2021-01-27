from Bio import SeqIO, Seq
import argparse
import dna_jellyfish as jellyfish

def get_kmer_presence(kmer,jelly_qlist):

    mer = jellyfish.MerDNA(kmer)
    mer.canonicalize()
    kmer_pres = []
    for i in range(len(jelly_qlist)):
        pres = int(jelly_qlist[i][mer]>0)
        kmer_pres.append(pres)

    sum_kmer_pres = sum(kmer_pres)
    return int(sum_kmer_pres>0)


def sequence_parser(l1jelly_qlist, l2jelly_qlist, l3jelly_qlist, lrjelly_qlist, assembly_kmers):

    pres_kmers = {}
    for kmer in assembly_kmers:
        if assembly_kmers[kmer] > 1:
            continue
        kmerR = str(Seq.Seq(kmer).reverse_complement())
        if kmerR != kmer and kmerR in assembly_kmers:
            continue

        if len(lrjelly_qlist)>0:
            kmer_pres_lr = get_kmer_presence(kmer,lrjelly_qlist)
            if kmer_pres_lr == 0:
                continue

        kmer_pres_l1 = get_kmer_presence(kmer,l1jelly_qlist)
        kmer_pres_l2 = get_kmer_presence(kmer,l2jelly_qlist)
        kmer_pres_l3 = get_kmer_presence(kmer,l3jelly_qlist)

        pres_kmers[kmer] = [kmer_pres_l1, kmer_pres_l2, kmer_pres_l3]

    return pres_kmers
       

def read_assembly(assembly, kmersize):
    
    fasta_sequences = SeqIO.parse(open(assembly),'fasta')
    assembly_kmers = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        for i in range(len(sequence) - kmersize + 1):
            kmer =  sequence[i:i+kmersize].upper()
            if 'N' in kmer:
                continue
            if kmer in assembly_kmers:
                assembly_kmers[kmer] += 1
            else:
                assembly_kmers[kmer] = 1
    return assembly_kmers


def writeOutput(assembly, kmersize, pres_kmers, outputFile):
         
    fasta_sequences = SeqIO.parse(open(assembly),'fasta')
    out = open(outputFile, 'w')
    
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        chrm = name.split(':')[0]
        start,end = name.split(':')[1].split('_')
        win_counts = [0,0,0,0]
        for i in range(len(sequence) - kmersize + 1):                
            kmer = sequence[i:i+kmersize].upper()
            if kmer not in pres_kmers:
                continue
            kmer_pres = pres_kmers[kmer]
            if kmer_pres[0] and not kmer_pres[1] and not kmer_pres[2]:
                win_counts[0] += 1   
            if kmer_pres[1] and not kmer_pres[0] and not kmer_pres[2]:
                win_counts[1] += 1
            if kmer_pres[2] and not kmer_pres[0] and not kmer_pres[1]:
                win_counts[2] += 1  
            if kmer_pres[0] or kmer_pres[1] or kmer_pres[2]:
                win_counts[3] += 1  
                             
        out.write(chrm + "\t" + start + '\t' + end +'\t' + str(win_counts[0]) + "\t" + str(win_counts[1]) + "\t" + str(win_counts[2]) + "\t" + str(win_counts[3])+ "\n")

    out.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Parse fasta file to check the origin of non-overlapping windows from different lineages.")

    sequence_params = parser.add_argument_group('Parameters of file to be parsed.')
    sequence_params.add_argument('-a', '--assembly',  required=True, help='Fasta file to parse.')

    jelly_params = parser.add_argument_group('Parameters of jellyfish files used for checking presence/absence.')
    jelly_params.add_argument('-k', '--kmersize', type=int, default = 51, help='Kmer size specified while making jellyfish.')
    jelly_params.add_argument('-l1c', '--lin1config',  required=True, help='Configuration file containing the path to the jellyfish file of each lineage 1 accession.')
    jelly_params.add_argument('-l2c', '--lin2config',  required=True, help='Configuration file containing the path to the jellyfish file of each lineage 2 accession.')
    jelly_params.add_argument('-l3c', '--lin3config',  required=True, help='Configuration file containing the path to the jellyfish file of each lineage 3 accession.')
    jelly_params.add_argument('-lrc', '--landraceconfig', help='Configuration file containing the path to the jellyfish file of each landrace.')

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', required=True, help='Path of output file.')

    args = parser.parse_args()

    l1jelly_qlist = []
    with open(args.lin1config,'r') as f:
        for l in f:
            l1jelly_qlist.append(jellyfish.QueryMerFile(l.strip()))

    l2jelly_qlist = []
    with open(args.lin2config,'r') as f:
        for l in f:
            l2jelly_qlist.append(jellyfish.QueryMerFile(l.strip()))

    l3jelly_qlist = []
    with open(args.lin3config,'r') as f:
        for l in f:
            l3jelly_qlist.append(jellyfish.QueryMerFile(l.strip()))

    lrjelly_qlist = []
    if args.landraceconfig is not None:      
        with open(args.landraceconfig,'r') as f:
            for l in f:
                lrjelly_qlist.append(jellyfish.QueryMerFile(l.strip()))

    assembly_kmers = read_assembly(args.assembly, args.kmersize)

    pres_kmers = sequence_parser(l1jelly_qlist, l2jelly_qlist, l3jelly_qlist, lrjelly_qlist, assembly_kmers)

    writeOutput(args.assembly, args.kmersize, pres_kmers, args.output)
