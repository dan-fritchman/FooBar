"""
"Fuel Injection Perfection"

Figure out the min # steps from `n` to 1,
given 2-3 available moves:
* Increase by one
* Decrease by one
* If `n` is even, divide by two
"""


def log(msg):
    print(msg)
    pass


def solution(n):
    """
    Not sure we could prove this is the shortest path.
    But it seems to work.
    Something like:
    * If we can, divide by 2.
    * If incrementing will get us to a multiple of 4, increment,
    * Otherwise, decrement.
    """

    try:
        nn = int(n)
    except:
        return 0
    if nn <= 1:
        return 0

    steps = 0
    while nn > 1:
        log(nn)
        steps += 1
        if not nn % 2:
            nn >>= 1
        elif nn != 3 and nn % 4 == 3:
            nn += 1
        else:
            nn -= 1
    return steps

