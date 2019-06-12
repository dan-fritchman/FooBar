import functools
from fractions import Fraction


def log(arg):
    print(arg)
    pass


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


class State:
    def __init__(self, graph, n, prob_to, prob_from):
        self.graph = graph
        self.n = n
        self.prob_to = prob_to
        self.prob_from = prob_from

    def walk(self):
        """ Probability of landing in this state.
        Valid only if *starting from* terminal states. """
        if self.n == 0:
            return Fraction(1, 1)
        prob = Fraction(0, 1)
        for (k, p) in enumerate(self.prob_from):
            if p > 0:
                from_state = self.graph.states[k]
                prob_from_state = from_state.walk()
                log(from_state)
                log(prob_from_state)
                contribution = Fraction(p * prob_from_state)
                log(contribution)
                prob += contribution
        return prob

    @property
    def terminal(self):
        return sum(self.prob_to) == 0

    def __repr__(self):
        return 'State ' + str(self.n) + str(self.prob_from)


class Graph:
    """ Graph model of the state-space.
    Most importantly sets up State objects with `prob_from` attributes. """

    def __init__(self, m):
        self.m = m
        self.rank = len(m)
        self.dens = [max(sum(row), 1) for row in m]
        self.states = []
        self.cols = []
        for k in range(len(m)):
            self.cols.append([row[k] for row in m])
        for n, row in enumerate(m):
            prob_from = [Fraction(self.cols[n][_], self.dens[_]) for _ in range(self.rank)]
            s = State(graph=self, n=n, prob_to=row, prob_from=prob_from)
            self.states.append(s)

    def probs(self):
        """ Probabilities of ending in each of our terminal states
        Returned as a list of fractions. """
        return [s.walk() for s in self.states if s.terminal]
        # results = []
        # for s in self.states:
        #     if s.terminal:
        #         results.append(s.walk())
        # return results


def solution(m):
    g = Graph(m)
    results = g.probs()

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

# case = [
#     [0, 1, 1, 0],
#     [1, 0, 0, 1],
#     [0, 0, 0, 0],  # Terminal
#     [0, 0, 0, 0],  # Terminal
# ]
# result = [1, 2, 3]
# assert (solution(case) == result)

# case = [
#     [0, 1, 0, 0, 0, 1],
#     [4, 0, 0, 3, 2, 0],
#     [0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0]
# ]
# assert (solution(case) == [0, 3, 2, 9, 14])
