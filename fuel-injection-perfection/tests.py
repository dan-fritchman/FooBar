import solution

# Not confusing at all.
solution = solution.solution

# Test Time!
assert (solution('4') == 2)
assert (solution('15') == 5)

# My Tests
assert (solution('1') == 0)
assert (solution('2') == 1)
assert (solution('3') == 2)
assert (solution('4') == 2)

assert (solution('8') == 3)
assert (solution('16') == 4)
assert (solution('32') == 5)
assert (solution('32') == 5)

for k in range(1100):
    assert (solution(str(2 ** k)) == k)
    print(f'2**{k} == {2 ** k} PASSED')
