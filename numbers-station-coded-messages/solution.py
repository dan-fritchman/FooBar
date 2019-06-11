def log(arg):
    print(arg)


def solution(l, t):
    """ Coded Messages
    Find first sub-list of `l` which sums up to `t`. """
    log(l)
    log(t)
    for start_index in range(len(l)):
        for end_index in range(start_index, len(l) + 1):
            sublist = l[start_index:end_index]
            log(sublist)
            if sum(sublist) == t:
                rval = [start_index, end_index - 1]
                log("Returning: " + str(rval))
                return rval
    rval = [-1, -1]
    log("Returning: " + str(rval))
    return rval


assert (solution(l=[4, 3, 5, 7, 8], t=12) == [0, 2])
assert (solution(l=[1, 2, 3, 4], t=15) == [-1, -1])
assert (solution(l=[1], t=1) == [0, 0])
assert (solution(l=[1, 2], t=3) == [0, 1])
assert (solution(l=[1, 2], t=4) == [-1, -1])
assert (solution(l=[], t=10) == [-1, -1])
assert (solution(l=[1, 2, 1, 2], t=3) == [0, 1])
