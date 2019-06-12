"""
Doomsday Fuel

A cleverly hidden Markov Chain problem.
Shout-out https://brilliant.org/wiki/absorbing-markov-chains/ for help with the math.
"""
import functools
from fractions import Fraction


def log(arg):
    print(arg)
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


#
# I = [[Fraction(1, 1), Fraction(0, 1)],
#      [Fraction(0, 1), Fraction(1, 1)]]
#
# for r in range(2,10):
#     print(r)
#     I = identity_matrix(r)
#     print(I)
#     inv = getMatrixInverse(I)
#     print(inv)
#     assert (inv == I)


def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    """Return lowest common multiple."""
    return a * b // gcd(a, b)


def lcmm(*args):
    """Return lcm of args."""
    return functools.reduce(lcm, args)


def zeroMatrix(x, y):
    """ Matrix full of zero-valued Fractions """

    def newRow(n):
        return [Fraction(0, 1) for _ in range(n)]

    return [newRow(y) for _ in range(x)]


def matrixMult(X, Y):
    rv = zeroMatrix(len(X), len(Y[0]))
    for i in range(len(X)):
        for j in range(len(Y[0])):
            for k in range(len(Y)):
                rv[i][j] += X[i][k] * Y[k][j]
    return rv


def solution(m):
    num_end_states = 0
    for row in m:
        if sum(row) == 0:
            num_end_states += 1
    num_trans_states = len(m) - num_end_states

    def den(row):
        return max(1, sum(row))

    S = zeroMatrix(len(m), len(m))
    for (x, row) in enumerate(m):
        d = den(row)
        S[x] = [Fraction(k, d) for k in row]

    # FIXME: might need some row/col re-orders
    Q = [row[:num_trans_states] for row in S[:num_trans_states]]
    R = [row[num_trans_states:] for row in S[:num_trans_states]]
    I = identity_matrix(num_trans_states)
    ImQ = zeroMatrix(num_trans_states, num_trans_states)
    for (x, row) in enumerate(Q):
        for (y, val) in enumerate(row):
            ImQ[x][y] = I[x][y] - Q[x][y]
    N = get_inverse(ImQ)
    M = matrixMult(N, R)

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
    assert (matrixMult(ImQ, N) == I)
    log("M=N*R:")
    log(M)

    results = M[0]

    log("RESULTS")
    log(results)

    result_lcm = functools.reduce(lcm, [f.denominator for f in results])
    result_nums = [f.numerator * (result_lcm // f.denominator) for f in results]
    rval = result_nums + [result_lcm]
    log("RVAL")
    log(rval)
    return rval


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
    [0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
]
result = [1, 0, 0, 0, 1]
assert (solution(case) == result)
