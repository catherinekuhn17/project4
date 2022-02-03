# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    This unit test compares the computed align and gap matrices to these same matrices 
    determined by hand, to see if they were filled out correctly in the algorithm
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10,-1)
    alignment = nw.align(seq1, seq2)
    # These matrices were determined by hand to have the correct values to check against computed matrices
    assert (nw.align_matrix == [[  0., -np.inf, -np.inf, -np.inf],
                               [-np.inf,   5., -11., -13.],
                               [-np.inf, -12.,   4.,  -8.],
                               [-np.inf, -12.,  -1.,   5.],
                               [-np.inf, -14.,  -6.,   4.]]).all()
    assert (nw.gapA_matrix == [[-10., -np.inf, -np.inf, -np.inf],
                               [-11., -22., -23., -24.],
                               [-12.,  -6., -17., -18.],
                               [-13.,  -7.,  -7., -18.],
                               [-14.,  -8.,  -8.,  -6.]]).all()
    assert (nw.gapB_matrix == [[-10., -11., -12., -13.],
                               [-np.inf, -22.,  -6.,  -7.],
                               [-np.inf, -23., -17.,  -7.],
                               [-np.inf, -24., -18., -12.],
                               [-np.inf, -25., -19., -17.]]).all()

def test_nw_backtrace():
    """
    This unit tests compares the output of the algorithm in aligning two 
    sequences, against the correct output.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10,-1)
    assert (17.0, 'MAVHQLIRRP', 'M---QLIRHP') == nw.align(seq3, seq4) # this is based off of a known example
   
    pass




