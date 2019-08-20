def solution(x, y):
    """ Prison Labor Dodgers
    Find the element in one list, but not in the other. """
    # Search the longer of the two lists.
    # If lengths are equal, things break.
    search = y
    compare = x
    if len(x) > len(y):
        search = x
        compare = y

    # Loop over the items of the longer list.
    # Keep an iteration "time-out", in case things run off the rails.
    iters = 999
    while len(search) and iters > 0:
        # The FooBar runtime apparently doesn't like "print".  Fine.
        # print(search)
        # print(compare)
        iters = iters - 1
        k = search.pop()
        # print(k)
        # If the last/ popped item is not in the other list, return it.
        if k not in compare:
            return k
        # Otherwise splice that item out of the other list.
        compare_index = compare.index(k)
        updates = compare[0:compare_index]
        updates.extend(compare[compare_index + 1:])
        compare = updates
    if not iters:
        raise Exception("Iterations Timed out for: " + str(x) + str(y))


# Test Cases

assert (solution(x=[13, 5, 6, 2, 5], y=[5, 2, 5, 13]) == 6)
assert (solution(x=[14, 27, 1, 4, 2, 50, 3, 1], y=[2, 4, -4, 3, 1, 1, 14, 27, 50]) == -4)
assert (solution(x=[1, 2], y=[1, 2, 3]) == 3)
assert (solution(x=[1, 2, 3], y=[1, 2]) == 3)
assert (solution(x=[1, 1, 1, 2, 3], y=[1, 1, 1, 2]) == 3)
assert (solution(x=[1], y=[]) == 1)
assert (solution(x=[], y=[1]) == 1)
assert (solution(x=[-1], y=[]) == -1)
assert (solution(x=[], y=[-1]) == -1)
