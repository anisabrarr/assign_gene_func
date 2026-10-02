import numpy as np


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

    """
    n = len(seq1)
    m = len(seq2)

    score_matrix = np.zeros((n + 1, m + 1))

    for i in range(1, n + 1):
        score_matrix[i][0] = score_matrix[i-1][0] + scoring_function(seq1[i-1], "-")

    for j in range(1, m + 1):
        score_matrix[0][j] = score_matrix[0][j-1] + scoring_function("-", seq2[j-1])

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = score_matrix[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            delete_score = score_matrix[i-1][j] + scoring_function(seq1[i-1], "-")
            insert_score = score_matrix[i][j-1] + scoring_function("-", seq2[j-1])
            score_matrix[i][j] = max(match_score, delete_score, insert_score)

    i, j = n, m
    align1 = []
    align2 = []

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        if i > 0 and j > 0:
            diag_score = score_matrix[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            if abs(current_score - diag_score) < 1e-5:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
                continue

        if i > 0:
            up_score = score_matrix[i-1][j] + scoring_function(seq1[i-1], "-")
            if abs(current_score - up_score) < 1e-5:
                align1.append(seq1[i-1])
                align2.append("-")
                i -= 1
                continue

        align1.append("-")
        align2.append(seq2[j-1])
        j -= 1

    final_align1 = "".join(reversed(align1))
    final_align2 = "".join(reversed(align2))
    final_score = float(score_matrix[n][m])

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


    """
    n = len(seq1)
    m = len(seq2)

    score_matrix = np.zeros((n + 1, m + 1))
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = score_matrix[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            delete_score = score_matrix[i-1][j] + scoring_function(seq1[i-1], "-")
            insert_score = score_matrix[i][j-1] + scoring_function("-", seq2[j-1])

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

        diag_score = score_matrix[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
        if abs(current_score - diag_score) < 1e-5:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
            continue

        up_score = score_matrix[i-1][j] + scoring_function(seq1[i-1], "-")
        if abs(current_score - up_score) < 1e-5:
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1
            continue

        align1.append("-")
        align2.append(seq2[j-1])
        j -= 1

    final_align1 = "".join(reversed(align1))
    final_align2 = "".join(reversed(align2))

    return final_align1, final_align2, float(max_score)


def scoring_function_simple(aa_i, aa_j):
    score = [-1, 5][aa_i == aa_j] if "-" not in (aa_i, aa_j) else -4
    return (score)