from Bio.Align import substitution_matrices


def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    gap_pen = 1
    

    matrix = []

    for i in range(len(seq1) + 1):
        row = []

        for j in range(len(seq2) + 1):
            row.append(0)

        matrix.append(row)


    for i in range(0,len(matrix[0])):
        matrix[0][i] = i*-gap_pen

    for i in range(0,len(matrix)):
        matrix[i][0] = i*-gap_pen

    

    for i in range(1,len(matrix)):
        for j in range(1,len(matrix[i])):
            up = matrix[i-1][j] - gap_pen
            left = matrix[i][j-1] - gap_pen
            diagonal = scoring_function(seq1[i-1], seq2[j-1]) + matrix[i-1][j-1]
            

            matrix[i][j] = max(up,left,diagonal)


    i = len(seq1)
    j = len(seq2)

 
    


    aligned_seq1 = ""
    aligned_seq2 = ""

    while i > 0 or j > 0:

        if i > 0 and j > 0:

            diagonal = (
                matrix[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            if matrix[i][j] == diagonal:
                aligned_seq1 += seq1[i - 1]
                aligned_seq2 += seq2[j - 1]

                i -= 1
                j -= 1
                continue

        if i > 0:

            up = matrix[i - 1][j] - gap_pen

            if matrix[i][j] == up:
                aligned_seq1 += seq1[i - 1]
                aligned_seq2 += "-"

                i -= 1
                continue

   
        if j > 0:
            aligned_seq1 += "-"
            aligned_seq2 += seq2[j - 1]

            j -= 1

    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    final_score = float(matrix[-1][-1])

    return aligned_seq1, aligned_seq2, final_score
    
   


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment(
    ...     "pending itch",
    ...     "unending glitch",
    ...     lambda x, y: [-1, 1][x == y]
    ... )
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.
    """

    gap_pen = 1

    
    matrix = []
    
    for i in range(len(seq1) + 1):
        row = []

        for j in range(len(seq2) + 1):
            row.append(0)

        matrix.append(row)

    for i in range(1, len(matrix)):
        for j in range(1, len(matrix[i])):

            up = matrix[i - 1][j] - gap_pen

            left = matrix[i][j - 1] - gap_pen

            diagonal = (
                matrix[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            
            matrix[i][j] = max(0, up, left, diagonal)

   
    max_val = 0
    i, j = 0, 0

    for k in range(len(matrix)):
        for l in range(len(matrix[k])):

            if matrix[k][l] > max_val:
                max_val = matrix[k][l]
                i, j = k, l

    
    aligned_seq1 = ""
    aligned_seq2 = ""

    
    while matrix[i][j] != 0:

       
        if i > 0 and j > 0:

            diagonal = (
                matrix[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            if matrix[i][j] == diagonal:
                aligned_seq1 += seq1[i - 1]
                aligned_seq2 += seq2[j - 1]

                i -= 1
                j -= 1
                continue

      
        if i > 0:

            up = matrix[i - 1][j] - gap_pen

            if matrix[i][j] == up:
                aligned_seq1 += seq1[i - 1]
                aligned_seq2 += "-"

                i -= 1
                continue

        
        if j > 0:

            left = matrix[i][j - 1] - gap_pen

            if matrix[i][j] == left:
                aligned_seq1 += "-"
                aligned_seq2 += seq2[j - 1]

                j -= 1
                continue

   
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    
    final_score = float(max_val)

    return aligned_seq1, aligned_seq2, final_score


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)



_blosum62_table = {}
_blosum62 = substitution_matrices.load("BLOSUM62")
for a in _blosum62.alphabet:
    for b in _blosum62.alphabet:
        _blosum62_table[(a, b)] = _blosum62[a][b]

def scoring_function_blosum62(aa_i, aa_j):
    return _blosum62_table[(aa_i, aa_j)]


print(global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y]))


print(local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y]))