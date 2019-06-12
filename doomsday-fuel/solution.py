"""
Doomsday Fuel

A cleverly hidden Markov Chain problem.
Shout-out https://brilliant.org/wiki/absorbing-markov-chains/ for help with the math.
"""
import functools
from fractions import Fraction


def log(mat):
    for row in mat:
        print(row)
    pass


def transposeMatrix(m):
    return list(map(list, zip(*m)))


def getMatrixMinor(m, i, j):
    return [row[:j] + row[j + 1:] for row in (m[:i] + m[i + 1:])]


def get_determinant(m):
    if len(m) == 1:
        return m[0][0]
    if len(m) == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1) ** c) * m[0][c] * get_determinant(getMatrixMinor(m, 0, c))
    return determinant


def get_inverse(m):
    if len(m) == 1:
        return [[1 / m[0][0]]]
    determinant = get_determinant(m)
    if len(m) == 2:
        return [[m[1][1] / determinant, -1 * m[0][1] / determinant],
                [-1 * m[1][0] / determinant, m[0][0] / determinant]]

    # find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMatrixMinor(m, r, c)
            cofactorRow.append(((-1) ** (r + c)) * get_determinant(minor))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c] / determinant
    return cofactors


def identity_matrix(rank):
    assert (rank >= 1)
    rval = []
    for r in range(rank):
        row = [Fraction(0, 1)] * rank
        row[r] = Fraction(1, 1)
        rval.append(row)
    return rval


def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    """Return lowest common multiple."""
    return (a * b) // gcd(a, b)


def lcmm(*args):
    """Return lcm of args."""
    return functools.reduce(lcm, args)


def zeroMatrix(x, y):
    """ Matrix full of zero-valued Fractions """

    def newRow(n):
        return [Fraction(0, 1) for _ in range(n)]

    return [newRow(y) for _ in range(x)]


def matrix_mult(X, Y):
    rv = zeroMatrix(len(X), len(Y[0]))
    for i in range(len(X)):
        for j in range(len(Y[0])):
            for k in range(len(Y)):
                rv[i][j] += X[i][k] * Y[k][j]
    return rv


def state_matrices(m):
    """ Parse input `m` into Markov form """
    def parse_row(n, row):
        row_type = 'trans'
        denom = sum(row)
        if denom == 1 and row[n] == 1:
            row_type = 'terminal'
        if denom == 0:
            denom = 1
            row_type = 'terminal'
        return dict(
            input_row_num=n,
            row_type=row_type,
            denom=denom,
        )

    parse_results = [parse_row(k, row) for (k, row) in enumerate(m)]

    # Figure out our order. Transient states all first, then terminal states.
    rows_in_order = []
    for row_parse in parse_results:
        if row_parse['row_type'] == 'trans':
            row_parse['our_row_num'] = len(rows_in_order)
            rows_in_order.append(row_parse)
    num_trans_states = len(rows_in_order)

    for row_parse in parse_results:
        if row_parse['row_type'] != 'trans':
            row_parse['our_row_num'] = len(rows_in_order)
            rows_in_order.append(row_parse)
    log(rows_in_order)

    row_map = {}
    for row in rows_in_order:
        if row['input_row_num'] in row_map:
            raise Exception
        row_map[row['input_row_num']] = row['our_row_num']

    S = zeroMatrix(len(m), len(m))
    for (x, row_data) in enumerate(rows_in_order):
        m_data = m[row_data['input_row_num']]
        den = row_data['denom']
        for (y, _) in enumerate(m_data):
            num = m_data[y]
            S[x][row_map[y]] = Fraction(num, den)

    Q = [row[:num_trans_states] for row in S[:num_trans_states]]
    R = [row[num_trans_states:] for row in S[:num_trans_states]]
    return S, Q, R


def matrix_subtract(lhs, rhs):
    assert (len(lhs) == len(rhs))
    rv = zeroMatrix(len(lhs), len(lhs[0]))
    # Very manual matrix subtraction
    for (x, row) in enumerate(lhs):
        for (y, val) in enumerate(row):
            rv[x][y] = lhs[x][y] - rhs[x][y]
    return rv


def format(results):
    """ Take our preferred result-format - a list of Fractions -
    and morph it into Google's - a list of numerators, followed by a shared denominator. """
    result_lcm = functools.reduce(lcm, [f.denominator for f in results])
    result_nums = [f.numerator * (result_lcm // f.denominator) for f in results]
    rval = result_nums + [result_lcm]
    log("RVAL")
    log(rval)
    return rval


def is_terminal_state(n, row):
    """ Boolean indication of whether row `n` represents a terminal/ final state. """
    if sum(row) == 0:
        return True
    if sum(row) == 1 and row[n] == 1:
        return True
    return False


def initial_terminal_state_method(m):
    """ Break out a special side-method if the initial state is also a terminal state.  Hmph. """
    num_final_states = 0
    for (n, row) in enumerate(m):
        if is_terminal_state(n, row):
            num_final_states += 1
    results = [Fraction(0, 1)] * num_final_states
    results[0] = Fraction(1, 1)
    return results


def solve(m):
    """ Markov-Style Solver """
    S, Q, R = state_matrices(m)
    I = identity_matrix(len(Q))
    ImQ = matrix_subtract(I, Q)
    N = get_inverse(ImQ)
    # Double-check that!
    Iback = matrix_mult(ImQ, N)
    assert (Iback == I)
    M = matrix_mult(N, R)

    log("m:")
    log(m)
    log("S:")
    log(S)
    log("Q:")
    log(Q)
    log("R:")
    log(R)
    # log(I)
    log("ImQ:")
    log(ImQ)
    log("N=ImQ**-1:")
    log(N)
    log("M=N*R:")
    log(M)

    results = M[0]
    return results


def solution(m):
    if (not len(m)) or (not len(m[0])):
        return []
    if is_terminal_state(0, m[0]):
        # Special case for the initial state, also being a final state
        results = initial_terminal_state_method(m)
    else:
        results = solve(m)

    log("RESULTS")
    log(results)

    return format(results)


case = [
    [0, 1, 1, 0, 0],
    [0, 0, 0, 1, 1],
    [0, 0, 0, 1, 1],
    [0, 0, 0, 0, 0],  # Terminal
    [0, 0, 0, 0, 0],  # Terminal
]
result = [1, 1, 2]
assert (solution(case) == result)

case = [
    [0, 2, 1, 0, 0],
    [0, 0, 0, 3, 4],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0]
]
assert (solution(case) == [7, 6, 8, 21])

case = [
    [0, 1, 1, 0, 0],
    [0, 0, 1, 1, 0],
    [0, 1, 0, 0, 1],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
]
result = [1, 1, 2]
assert (solution(case) == result)

case = [
    [0, 1, 1, 0],
    [0, 0, 1, 1],
    [0, 1, 0, 1],
    [0, 0, 0, 0],
]
result = [1, 1]
assert (solution(case) == result)

case = [
    [0, 1, 0, 0],
    [0, 0, 1, 1],
    [0, 1, 0, 1],
    [0, 0, 0, 0],
]
result = [1, 1]
assert (solution(case) == result)

case = [
    [0, 1, 1, 0],
    [1, 0, 0, 1],
    [0, 0, 0, 0],  # Terminal
    [0, 0, 0, 0],  # Terminal
]
result = [2, 1, 3]
assert (solution(case) == result)

case = [
    [0, 1, 0, 0, 0, 1],
    [4, 0, 0, 3, 2, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
]
assert (solution(case) == [0, 3, 2, 9, 14])

case = [
    [0, 1],
    [0, 0],
]
result = [1, 1]
assert (solution(case) == result)

case = [
    [1, 1],
    [0, 0],
]
result = [1, 1]
assert (solution(case) == result)

case = [
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
]
result = [1, 0, 0, 0, 1]
assert (solution(case) == result)

case = [
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
]
result = [1, 0, 0, 0, 1]
assert (solution(case) == result)

case = [
    [0, 0, 0, 0, 0, 0],
    [0, 11, 0, 0, 0, 71],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [4, 0, 0, 3, 2, 0],
    [0, 0, 0, 0, 0, 0]
]
assert (solution(case) == [1, 0, 0, 0, 1])

case = [
    [0, 0],
    [1, 1],
]
result = [1, 1]
assert (solution(case) == result)

case = [
    [0, 0, 0],
    [1, 1, 1],
    [1, 1, 1],
]
result = [1, 1]
assert (solution(case) == result)


case = [
    [1, 0, 0],
    [1, 1, 1],
    [1, 1, 1],
]
result = [1, 1]
assert (solution(case) == result)


case = [
    [0, 511, 511, 511],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
]
result = [1, 1, 1, 3]
assert (solution(case) == result)

case = [
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [1, 1, 1, 1],
    [1, 1, 1, 1],
]
result = [1, 0, 1]
assert (solution(case) == result)

case = [
    [0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
]
result = [0, 1, 0, 0, 1]
assert (solution(case) == result)

case = [
    [0, 0],
    [0, 0],
]
result = [1, 0, 1]
assert (solution(case) == result)

case = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
]
result = [1, 0, 0, 1]
assert (solution(case) == result)

case = [
    [0],
]
result = [1, 1]
assert (solution(case) == result)

case = [
    [1],
]
result = [1, 1]
assert (solution(case) == result)

case = [
    [],
]
result = []
assert (solution(case) == result)

case = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
]
result = [1, 0, 0, 1]
assert (solution(case) == result)

case = [
    [0, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
]
result = [1, 0, 0, 1]
assert (solution(case) == result)

assert (
        solution([
            [0, 2, 1, 0, 0],
            [0, 0, 0, 3, 4],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0]
        ]) == [7, 6, 8, 21]
)

assert (
        solution([
            [0, 1, 0, 0, 0, 1],
            [4, 0, 0, 3, 2, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0]
        ]) == [0, 3, 2, 9, 14]
)

assert (
        solution([
            [1, 2, 3, 0, 0, 0],
            [4, 5, 6, 0, 0, 0],
            [7, 8, 9, 1, 0, 0],
            [0, 0, 0, 0, 1, 2],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0]
        ]) == [1, 2, 3]
)
assert (
        solution([
            [0]
        ]) == [1, 1]
)

assert (
        solution([
            [0, 0, 12, 0, 15, 0, 0, 0, 1, 8],
            [0, 0, 60, 0, 0, 7, 13, 0, 0, 0],
            [0, 15, 0, 8, 7, 0, 0, 1, 9, 0],
            [23, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [37, 35, 0, 0, 0, 0, 3, 21, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]) == [1, 2, 3, 4, 5, 15]
)

assert (
        solution([
            [0, 7, 0, 17, 0, 1, 0, 5, 0, 2],
            [0, 0, 29, 0, 28, 0, 3, 0, 16, 0],
            [0, 3, 0, 0, 0, 1, 0, 0, 0, 0],
            [48, 0, 3, 0, 0, 0, 17, 0, 0, 0],
            [0, 6, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]) == [4, 5, 5, 4, 2, 20]
)

assert (
        solution([
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]) == [1, 1, 1, 1, 1, 5]
)

assert (
        solution([
            [1, 1, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]) == [2, 1, 1, 1, 1, 6]
)

assert (
        solution([
            [0, 86, 61, 189, 0, 18, 12, 33, 66, 39],
            [0, 0, 2, 0, 0, 1, 0, 0, 0, 0],
            [15, 187, 0, 0, 18, 23, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]) == [6, 44, 4, 11, 22, 13, 100]
)

assert (
        solution([
            [0, 0, 0, 0, 3, 5, 0, 0, 0, 2],
            [0, 0, 4, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 4, 4, 0, 0, 0, 1, 1],
            [13, 0, 0, 0, 0, 0, 2, 0, 0, 0],
            [0, 1, 8, 7, 0, 0, 0, 1, 3, 0],
            [1, 7, 0, 0, 0, 0, 0, 2, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]) == [1, 1, 1, 2, 5]
)