# https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
def Alignment(seq1, seq2, gap=1, substitution=1):
    F = [[0 for col in range(len(seq2) + 1)] for row in range(len(seq1) + 1)]

    for i in range(len(seq1) + 1):
        F[i][0] = gap * i

    for j in range(len(seq2) + 1):
        F[0][j] = gap * j

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            penalty = 0 if seq1[i - 1] == seq2[j - 1] else substitution
            match = F[i - 1][j - 1] + penalty
            delete = F[i - 1][j] + gap
            insert = F[i][j - 1] + gap
            F[i][j] = min(match, insert, delete)

    ali1, ali2 = "", ""
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        try:
            penalty = 0 if seq1[i - 1] == seq2[j - 1] else substitution
        except:
            pass

        if j > 0 and i > 0 and F[i][j] == F[i - 1][j - 1] + penalty:
            ali1 = seq1[i - 1] + ali1
            ali2 = seq2[j - 1] + ali2
            i -= 1
            j -= 1
        elif i > 0 and F[i][j] == F[i - 1][j] + gap:
            ali1 = seq1[i - 1] + ali1
            ali2 = "_" + ali2
            i -= 1
        else:
            ali1 = "_" + ali1
            ali2 = seq2[j - 1] + ali2
            j -= 1

    return F[len(seq1)][len(seq2)], [ali1, ali2]
