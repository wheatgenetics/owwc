import argparse
import dna_jellyfish as jellyfish

def get_kmer_presence(kmerF,jelly_qlist,idx,minCount):

    mer = jellyfish.MerDNA(kmerF)
    mer.canonicalize()
    kmer_pres = []
    for i in range(len(jelly_qlist)):
        pres = int(jelly_qlist[i][mer]>0)
        kmer_pres.append(pres)
        if i < idx:
            if pres:
                return None
        if i == idx:
            if not pres:
                return None

    sum_kmer_pres = sum(kmer_pres)
    if sum_kmer_pres > len(kmer_pres)-minCount or sum_kmer_pres < minCount:
        return None

    return kmer_pres


def jfdump_parser(jelly_qlist, idx, jfdump, kmerSize, minCount, output):
    out = open(output, 'w')
    kmers_dump = open(jfdump, 'r')
    for l in kmers_dump: 
        kmerF = l.split()[0].upper()
        if 'N' in kmerF:
            continue
        kmer_pres=get_kmer_presence(kmerF,jelly_qlist,idx,minCount)
        if kmer_pres is not None: 
            kmer_pres_str=[str(pres) for pres in kmer_pres]
            out.write(kmerF+'\t'+''.join(kmer_pres_str)+'\n')

    kmers_dump.close()
    out.close()
       



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Parse jellyfish dump file to check presence/absence of k-mers in diversity panel.")

    parser.add_argument('-a', '--accname',  required=True, help='Name of the accession containing the sequence.')
    parser.add_argument('-j', '--jfdump',  required=True, help='Jellyfish dump file to parse.')
    parser.add_argument('-c', '--config',  required=True, help='Configuration file containing the name of each accession tab-separated from the path to the jellyfish of that line.')
    parser.add_argument('-k', '--kmersize', type=int, default = 51, help='Kmer size specified while making jellyfish.')
    parser.add_argument("-mc", "--mincount", type=int, default = 2, help="only those k-mers are retained which are present/absent in more than this number of accessions.")

    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', help='Path of output file.')

    args = parser.parse_args()

    jelly_qlist = []
    count = 0
    with open(args.config,'r') as f:
        for l in f:
            lvals = l.strip().split()
            if lvals[0]==args.accname:
                idx = count
            jelly_qlist.append(jellyfish.QueryMerFile(lvals[1]))
            count += 1

    jfdump_parser(jelly_qlist, idx, args.jfdump, args.kmersize, args.mincount, args.output)
