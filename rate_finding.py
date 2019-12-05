import read_class
import math

from read_class import Read
from typing import Callable, dict, list, tuple, Union


def find_rates(snps: list[int], reads: list[Read], rate: float) -> dict[int, float]:
    ##if sufficient coverage, approximation works well
    ##snps must be sorted
    lb_for_approx = 200
    counts = get_counts(snps, reads)
    m = max(map(max, counts.values()))

    if m == 0:
        print("should not happen", snps, reads)
    elif m > lb_for_approx:
        r = approx_rate(list(counts.values()))
    else:
        r = NW(counts, rate)

    return {0: r, 1: 1 - r}  ##had swapped


def approx_rate(counts: list[list[int]]) -> float:
    l = len(counts)
    maxes = map(max, counts)
    sums = map(sum, counts)
    # if sum(sums) == 0:
    #        print counts
    return float(sum(maxes)) / float(sum(sums))


def get_counts(
    snps: list[int], reads: list[Read], shuffle: bool = True
) -> dict[int, list[int]]:
    if shuffle:
        reads = read_class.sample_from_reads(reads)
    counts = {i: [0, 0] for i in range(len(snps))}
    back = {}
    for i in range(len(snps)):
        back[snps[i]] = i
    for R in reads:
        for snp in R.keys:
            counts[back[snp]][R.read[snp] % 2] += R.count
    return counts


def der(f: Callable, x: float, dx: float) -> float:
    return (f(x + dx) - f(x - dx)) / (2 * float(dx))


def ld(f, x, dx):
    return f(x + dx) - f(x)


def rd(f, x, dx):
    return f(x) - f(x - dx)


def NW(vec: dict[int, list[int]], rate: float) -> float:
    left = 0.001
    right = 1 - left
    dx = 0.001
    d = 1

    def f(x):
        FD = {}
        return math.log(
            forward(FD, 0, vec, len(vec) - 1, rate, x)
            + forward(FD, 1, vec, len(vec) - 1, rate, x)
        )

    if f(0.501) == f(0.499):
        if f(0.5) >= f(0.499):
            return 0.5
        else:
            left = 0.5 + 0.001

    while abs(d) > 0.001:
        mid = (right + left) / 2.0
        d = der(f, mid, dx)
        if d > 0:
            left = mid
        else:
            right = mid

        if mid > 0.98:
            return 0.98
        if mid < 0.02:
            return 0.02
    return mid


def trans(s1: int, s2: int, rate: float) -> float:
    if s1 == s2:
        return rate
    else:
        return 1 - rate


def forward(
    FD: dict[tuple[int, int], Union[float, int]],
    X: int,
    vec: dict[int, list[int]],
    i: int,
    rate: float,
    p: float,
) -> Union[int, float]:
    if (X, i) in FD:
        return FD[(X, i)]
    e = 0.0001
    p = p * (1 - e) + (1 - p) * (e)
    q = 1 - p
    if i == 0:
        Y = vec[i]
        if X == 1:
            FD[(1, 0)] = 0
            return 0
        else:
            val = choose(Y[0] + Y[1], Y[0]) * (p ** Y[0]) * (q ** Y[1])
            FD[(0, 0)] = 10 ** 10 * val
            return 10 ** 10 * val
    else:
        Y = vec[i]
        f0 = forward(FD, 0, vec, i - 1, rate, p)
        f1 = forward(FD, 1, vec, i - 1, rate, p)
        if (0, i - 1) not in FD:
            FD[(0, i - 1)] = (10 ** 10) * f0
        if (1, i - 1) not in FD:
            FD[(1, i - 1)] = (10 ** 10) * f1

        SUM_log = math.log(trans(0, X, rate) * f0 + trans(1, X, rate) * f1)
        sample_prob_log = (Y[X] * math.log(p)) + (Y[1 - X] * math.log(q))
        binom_log = math.log(choose(sum(Y), Y[0]))
        f_val = math.exp(SUM_log + sample_prob_log + binom_log)

        if (X, i) not in FD:
            FD[(X, i)] = 10 ** 10 * f_val
        return 10 ** 10 * f_val


def choose(n: int, k: int) -> int:
    # n CHOOSE k
    j = min(k, n - k)
    numerator = 1
    for i in range(j):
        numerator = numerator * (n - i)
    return numerator // math.factorial(j)
