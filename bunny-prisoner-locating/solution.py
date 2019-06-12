"""
The "Cell Layout":
| 11
| 7  12
| 4  8  13
| 2  5  9  14
| 1  3  6  10 15
"""


def log(arg):
    # Reminder foo.bar don't like print!
    print(arg)
    pass


def solution(x, y):
    """ Retrieve the value at location (x,y) in that cell layout above. """
    # Use something like the radial axis, or what we'll call a "diagonal stripe"
    # First create the zero-referenced row & column
    x0 = x - 1
    y0 = y - 1
    # x0, y0 = x - 1, y - 1
    # The "stripe", or "radius", is just the sum of the two.
    stripe = x0 + y0
    log(stripe)
    # Sort out the first/ leftmost entry in the stripe
    stripe_start = 1 + sum([(k + 1) for k in range(stripe)])
    log(stripe_start)
    # The answer just is that leftmost entry, plus the x-offset
    # And remember they want it as a string!
    return str(stripe_start + x0)


def solution(x, y):
    """ Retrieve the value at location (x,y) in that cell layout above. """
    # One-liner edition:
    return str(x + sum([(k + 1) for k in range(x + y - 2)]))


# Google's (Public) Tests
assert (solution(5, 10) == '96')

# My Tests
assert (solution(1, 1) == "1")
assert (solution(1, 2) == "2")
assert (solution(2, 1) == "3")
assert (solution(1, 3) == "4")
assert (solution(2, 2) == "5")
assert (solution(3, 1) == "6")
assert (solution(1, 4) == "7")
assert (solution(2, 3) == "8")
assert (solution(3, 2) == "9")
assert (solution(4, 1) == "10")
assert (solution(1, 5) == "11")
assert (solution(2, 4) == "12")
assert (solution(3, 3) == "13")
assert (solution(4, 2) == "14")
assert (solution(5, 1) == "15")
