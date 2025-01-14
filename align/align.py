# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    
    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        This method aligns two sequences using the Needleman-Wunsch algorithm. Match score is determined by a scoring matrix,
        and gap opening/extension scores were determined when the NeedlemanWunsch object was created. Both the alignment score
        and sequence alignments are ultimately determined
        
        Parameters
        ----------
            seqA : str
                The first sequence to be aligned
            seqB : str
                The second sequence to be aligned
        
        Returns
        -------
            alignment : tuple
                A tuple containing the final alignment score, as determined by the algorithm, as well as the alignment of 
                the two sequences, formated (alignment score, alignment seq 1, alignment seq 2)
                
        """
        self._seqA = seqA
        self._seqB = seqB
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # set score matrix 0,0 to 0
        self._align_matrix[0][0] = 0 
        
        # set first columns as gap A matrix to extending gaps, and first row of gap B matrix to extending gaps
        self._gapA_matrix[:,0] = [self.gap_open + i*self.gap_extend 
                                      for i in range(len(self._align_matrix))] 
        self._gapB_matrix[0] = [self.gap_open + i*self.gap_extend 
                                      for i in range(len(self._align_matrix[0]))] 
        
        # creation of all matrices
        for i in np.arange(1, len(self._align_matrix)): # going by rows (so A)
            for j in np.arange(1, len(self._align_matrix[0])):# going by columns (so B)
                # filling out score:
                matrix_vals = [self._align_matrix[i-1][j-1], # M
                               self._gapA_matrix[i-1][j-1], # A
                               self._gapB_matrix[i-1][j-1]] # B
                max_val = max(matrix_vals) # find which of M, A, or B is max
                self._back[i][j] = matrix_vals.index(max_val) # assign which one gives max (0 for M, 1 for A, 2 for B) 
                # --> additional note, if the max is the same for more than one (say A and M) this index call returns 
                # the index that occurs first, which would be M in that case (as M is first in the list). This assures 
                # that we follow the highroad alignment as alignment is prioritized.
                self._align_matrix[i][j] = self.sub_dict[(seqA[i-1], seqB[j-1])] + max_val # fill in new M value

                # A is in the rows and so we want to look at i-1 (the previous row!)
                matrix_vals_A = [self.gap_open + self.gap_extend + self._align_matrix[i-1][j], # M
                                 self.gap_extend + self._gapA_matrix[i-1][j], # A
                                 self.gap_open + self.gap_extend + self._gapB_matrix[i-1][j]] # B      
                max_val_A = max(matrix_vals_A) # score is the max of these
                self._back_A[i][j] = matrix_vals_A.index(max_val_A) # assign which one gives max (0 for M, 1 for A, 2 for B)
                self._gapA_matrix[i][j] = max_val_A # fill in new A value

                # B is in the columns, so we want to look at j-1 (the previous column!)
                matrix_vals_B = [self.gap_open + self.gap_extend + self._align_matrix[i][j-1], # M
                                 self.gap_open + self.gap_extend + self._gapA_matrix[i][j-1], # A
                                 self.gap_extend + self._gapB_matrix[i][j-1]] # B
                max_val_B = max(matrix_vals_B) # score is the max of these
                self._back_B[i][j] = matrix_vals_B.index(max_val_B) # # assign which one gives max (0 for M, 1 for A, 2 for B)
                self._gapB_matrix[i][j] = max_val_B # fill in new B value
                
        final_score_all = [self._align_matrix[-1][-1], self._gapA_matrix[-1][-1], self._gapB_matrix[-1][-1]]
        final_score = max(final_score_all) # final score is the max of the last box
        self._final_score_idx = final_score_all.index(final_score) # determine if final score came from M, A, or B
        self._alignment_score = final_score

        #some additional setup of the backtrace matrix (filling in what we haven't so far)
        self._back_A[:,0] = np.ones(len(self._back_B[:,0]), dtype=int)*1 # can only come from A (no more left of B) 
        # --> also from how we filled this column (gap open + i*(gap extend))
        self._back_B[0] =  np.ones(len(self._back_A[0]), dtype=int)*2 # can only come from B (no more left of A) 
        # --> also from how we filled this row (gap open + i*(gap extend))
        
        # making these callable outside for unit tests
        self.align_matrix = self._align_matrix 
        self.gapA_matrix = self._gapA_matrix
        self.gapB_matrix = self._gapB_matrix 

        alignment = self._backtrace()
        return alignment

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        This method backtraces through the three backtrace matrices created in the method align to 
        create the final alignment.
        
        Returns:
        -------
            alignment_tup : tuple
                a tuple formatted (alignment score, alignment for seq A, alignment for seq B)
        
        """        
        i = -1 # we want to go backwards with i and j until we reach the beginning (start at the end of the sequences)
        j = -1        
        idx = self._final_score_idx # we are starting here (M, A, or B), determined from where final score came from.
        while abs(i)<=len(self._seqA)+1 and abs(j)<=len(self._seqB)+1:
            if abs(j) == len(self._seqB)+1 and abs(i) == len(self._seqA)+1: # cond'n that we have completed the alignment!
                break
            else:
                if idx == 0: # max came from M
                    self.seqA_align = self._seqA[i] + self.seqA_align # continue in A
                    self.seqB_align = self._seqB[j] + self.seqB_align # continue in B
                    idx = self._back[i][j] # where to go next
                    i-=1 # go back 1 in i direction (A)
                    j-=1 # go back 1 in j direction (B)

                elif idx == 1: # max came from A
                    self.seqB_align = '-' + self.seqB_align # inserting gap into B alignment
                    self.seqA_align = self._seqA[i] + self.seqA_align # continue in A
                    idx = self._back_A[i][j] # where to go next
                    i-=1 # go back in A direction

                elif idx == 2: # max came from B
                    self.seqA_align = '-' + self.seqA_align  # inserting gap into A alignment
                    self.seqB_align = self._seqB[j]+ self.seqB_align # continue in B
                    idx = self._back_B[i][j] # where to go next
                    j-=1 # go back in B direction
        alignment_tup =  (self._alignment_score, self.seqA_align, self.seqB_align)
        return alignment_tup


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header





