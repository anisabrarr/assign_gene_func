import numpy as np
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

    blosum62 = substitution_matrices.load("BLOSUM62")
    n = len(seq1)
    m = len(seq2)
    gap_penalty = -4

    score_matrix = np.zeros((n + 1, m + 1))
    for i in range(n + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(m + 1):
        score_matrix[0][j] = j * gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair = (seq1[i-1], seq2[j-1])
            match_score = score_matrix[i-1][j-1] + blosum62.get(pair, blosum62.get((pair[1], pair[0]), -4))
            delete_score = score_matrix[i-1][j] + gap_penalty
            insert_score = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match_score, delete_score, insert_score)

    i, j = n, m
    align1 = []
    align2 = []

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]
        if i > 0 and j > 0:
            pair = (seq1[i-1], seq2[j-1])
            m_score = blosum62.get(pair, blosum62.get((pair[1], pair[0]), -4))
            if abs(current_score - (score_matrix[i-1][j-1] + m_score)) < 1e-5:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
                continue

        if i > 0 and abs(current_score - (score_matrix[i-1][j] + gap_penalty)) < 1e-5:
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1
        else:
            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1

    final_align1 = "".join(reversed(align1))
    final_align2 = "".join(reversed(align2))
    final_score = score_matrix[n][m]

    return final_align1, final_align2, final_score
    

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
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    blosum62 = substitution_matrices.load("BLOSUM62")
    n = len(seq1)
    m = len(seq2)
    gap_penalty = -4

    score_matrix = np.zeros((n + 1, m + 1))
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            pair = (seq1[i-1],seq2[j-1])
            match_score = score_matrix[i-1][j-1] + blosum62.get(pair, blosum62.get((pair[1], pair[0]), -4))
            delete_score = score_matrix[i-1][j] + gap_penalty
            insert_score = score_matrix[i][j-1] + gap_penalty

            score = max(0, match_score, delete_score, insert_score)
            score_matrix[i][j] = score

            if score > max_score:
                max_score = score
                max_pos = (i, j)

    i, j = max_pos
    align1 = []
    align2 = []

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current_score = score_matrix[i][j]

        pair = (seq1[i-1], seq2[j-1])
        m_score = blosum62.get(pair, blosum62.get((pair[1], pair[0]), -4))
        if abs(current_score - (score_matrix[i-1][j-1] + m_score)) < 1e-5:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
            continue

        if abs(current_score - (score_matrix[i-1][j] + gap_penalty)) < 1e-5:
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1
        else:
            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1

    final_align1 = "".join(reversed(align1))
    final_align2 = "".join(reversed(align2))

    return final_align1, final_align2, max_score


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
