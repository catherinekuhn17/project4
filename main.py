# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")
    
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10,-1)
    align_dict = {
       'gg' : nw.align(hs_seq, gg_seq),
       'mm' : nw.align(hs_seq, mm_seq),
       'br' : nw.align(hs_seq, br_seq),
       'tt' : nw.align(hs_seq, tt_seq)
    }
    ordering = np.argsort([e[2] for e in align_dict.values()])[::-1]
    for idx in ordering:
        print(f'species: {list(align_dict.keys())[idx]}')
        print(f'score: {list(align_dict.values())[idx][2]}')
        

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    pass

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    pass

if __name__ == "__main__":
    main()
