'''
Project 4
'''

def build_scoring_matrix(alphabet='AGCT',
                         diag_score=10, off_diag_score=4, dash_score=6):
    '''
    input:
        string of alphabet         alphabet
        score of matches           diag_score
        score of all mismatches    off_diag_score
        score of gaps              dash_score

    function returns scoring matrix as a dict of dicts
    '''
    alpha = ['-'] + list(alphabet)
    res = {dmy_i: {dmy_j: off_diag_score for dmy_j in alpha} for dmy_i in alpha}

    for idx in alpha:
        res[idx][idx] = diag_score
        res['-'][idx] = dash_score
        res[idx]['-'] = dash_score

    return res


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix=build_scoring_matrix(),
                             global_flag=True):
    '''
    input:
        sequences x and y          seq_x, seq_y
        scoring_matrix built by function build_scoring_matrix
        global or local alignment  global_flag

    function returns alignment matrix as list of lists
    '''
    num_x = len(seq_x) + 1
    num_y = len(seq_y) + 1
    matrix = [[0 for _ in range(num_y)] for _ in range(num_x)]

    for idx_x in range(1, num_x):
        matrix[idx_x][0] = matrix[idx_x-1][0] + scoring_matrix[ seq_x[idx_x-1] ]['-']
        if not global_flag:
            matrix[idx_x][0] = max(matrix[idx_x][0], 0)

    for idx_y in range(1, num_y):
        matrix[0][idx_y] = matrix[0][idx_y-1] + scoring_matrix['-'][ seq_y[idx_y-1] ]
        if not global_flag:
            matrix[0][idx_y] = max(matrix[0][idx_y], 0)

    for idx_x in range(1, num_x):
        for idx_y in range(1, num_y):
            vert = matrix[idx_x-1][idx_y] + scoring_matrix[ seq_x[idx_x-1] ]['-']
            horiz = matrix[idx_x][idx_y-1] + scoring_matrix['-'][ seq_y[idx_y-1] ]
            diag = matrix[idx_x-1][idx_y-1] + scoring_matrix[ seq_x[idx_x-1] ][ seq_y[idx_y-1] ]
            matrix[idx_x][idx_y] = max(vert, horiz, diag)
            if not global_flag:
                matrix[idx_x][idx_y] = max(matrix[idx_x][idx_y], 0)

    return matrix


def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    '''
    function returns alignment of input sequences
    '''
    idx_x = len(seq_x)
    idx_y = len(seq_y)
    score = alignment_matrix[idx_x][idx_y]

    alignment_x = ''
    alignment_y = ''

    while idx_x != 0 and idx_y != 0:
        current_score = alignment_matrix[idx_x][idx_y]
        if current_score == alignment_matrix[idx_x-1][idx_y-1] +\
                            scoring_matrix[ seq_x[idx_x-1] ][ seq_y[idx_y-1] ]:
            alignment_x = seq_x[idx_x-1] + alignment_x
            alignment_y = seq_y[idx_y-1] + alignment_y

            idx_x -= 1
            idx_y -= 1

        elif current_score == alignment_matrix[idx_x-1][idx_y] +\
                              scoring_matrix[ seq_x[idx_x-1] ]['-']:
            alignment_x = seq_x[idx_x-1] + alignment_x
            alignment_y = '-' + alignment_y

            idx_x -= 1

        else:
            alignment_x = '-' + alignment_x
            alignment_y = seq_y[idx_y-1] + alignment_y

            idx_y -= 1

    while idx_x != 0:
        alignment_x = seq_x[idx_x-1] + alignment_x
        alignment_y = '-' + alignment_y
        idx_x -= 1

    while idx_y != 0:
        alignment_x = '-' + alignment_x
        alignment_y = seq_y[idx_y-1] + alignment_y
        idx_y -= 1

    return (score, alignment_x, alignment_y)


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    '''
    function returns alignment of input sequences
    '''
    max_score = [-float('Inf'), 0, 0]

    for idx_x in range(len(alignment_matrix)):
        for idx_y in range(len(alignment_matrix[0])):
            if alignment_matrix[idx_x][idx_y] > max_score[0]:
                max_score = [alignment_matrix[idx_x][idx_y], idx_x, idx_y]

    score, idx_x, idx_y = max_score
    alignment_x = ''
    alignment_y = ''

    while idx_x != 0 and idx_y != 0:
        current_score = alignment_matrix[idx_x][idx_y]

        if current_score <= 0:
            break

        if current_score == alignment_matrix[idx_x-1][idx_y-1] +\
                            scoring_matrix[ seq_x[idx_x-1] ][ seq_y[idx_y-1] ]:
            alignment_x = seq_x[idx_x-1] + alignment_x
            alignment_y = seq_y[idx_y-1] + alignment_y
            idx_x -= 1
            idx_y -= 1
        elif current_score == alignment_matrix[idx_x-1][idx_y] +\
                              scoring_matrix[ seq_x[idx_x-1] ]['-']:
            alignment_x = seq_x[idx_x-1] + alignment_x
            alignment_y = '-' + alignment_y
            idx_x -= 1
        else:
            alignment_x = '-' + alignment_x
            alignment_y = seq_y[idx_y-1] + alignment_y
            idx_y -= 1

    return (score, alignment_x, alignment_y)
