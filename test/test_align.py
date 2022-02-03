# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10,-1)
    alignment = nw.align(seq1, seq2)
    assert nw.align_matrix == array([[  0., -inf, -inf, -inf, -inf, -inf, -inf, -inf],
                                     [-inf,   5., -11., -10., -12., -15., -17., -18.],
                                     [-inf, -12.,   4.,  -7.,  -8.,  -9., -11., -11.],
                                     [-inf, -11.,  -8.,   5.,  -4., -11., -12., -12.],
                                     [-inf, -15.,  -7., -10.,   2.,  -4.,   1., -10.],
                                     [-inf, -14.,  -3.,  -9.,  -9.,   3.,  -4.,   0.],
                                     [-inf, -13., -11.,   1.,  -5., -11.,   0.,  -7.],
                                     [-inf, -15., -13.,  -8.,   5.,  -8., -11.,  -3.],
                                     [-inf, -18., -10., -13., -11.,  10.,  -6.,  -9.],
                                     [-inf, -19., -11., -12., -13.,  -1.,  10.,  -3.],
                                     [-inf, -21., -14., -14., -14.,  -9.,  -3.,  17.]])
    
    assert nw.gapA_matrix == array([[-10., -inf, -inf, -inf, -inf, -inf, -inf, -inf],
                                    [-11., -22., -23., -24., -25., -26., -27., -28.],
                                    [-12.,  -6., -17., -18., -19., -20., -21., -22.],
                                    [-13.,  -7.,  -7., -18., -19., -20., -21., -22.],
                                    [-14.,  -8.,  -8.,  -6., -15., -18., -19., -20.],
                                    [-15.,  -9.,  -9.,  -7.,  -9., -15., -10., -21.],
                                    [-16., -10., -10.,  -8., -10.,  -8., -11., -11.],
                                    [-17., -11., -11.,  -9., -11.,  -9., -11., -12.],
                                    [-18., -12., -12., -10.,  -6., -10., -12., -13.],
                                    [-19., -13., -13., -11.,  -7.,  -1., -12., -13.],
                                    [-20., -14., -14., -12.,  -8.,  -2.,  -1., -12.]])
    assert nw.gapB_matrix == array([[-10., -11., -12., -13., -14., -15., -16., -17.],
                                    [-inf, -22.,  -6.,  -7.,  -8.,  -9., -10., -11.],
                                    [-inf, -23., -17.,  -7.,  -8.,  -9., -10., -11.],
                                    [-inf, -24., -18., -18.,  -6.,  -7.,  -8.,  -9.],
                                    [-inf, -25., -19., -18., -17.,  -9., -10., -10.],
                                    [-inf, -26., -20., -14., -15., -16.,  -8.,  -9.],
                                    [-inf, -27., -21., -21., -10., -11., -12., -11.],
                                    [-inf, -28., -22., -22., -19.,  -6.,  -7.,  -8.],
                                    [-inf, -29., -23., -21., -21., -17.,  -1.,  -2.],
                                    [-inf, -30., -24., -22., -22., -18., -12.,  -1.],
                                    [-inf, -31., -25., -25., -23., -19., -13., -12.]])

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10,-1)
    assert (17.0, 'MAVHQLIRRP', 'M---QLIRHP') == nw.align(seq3, seq4) # this is based off of a known example
   
    pass




