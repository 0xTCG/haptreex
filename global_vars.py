import prescoring
import math

from typing import List

global cutoff
cutoff = 0

global quality_cutoff
quality_cutoff = 10

global confidence
confidence = 0.5

global pair_thresh
pair_thresh = 0.8

global V
global vcfChroms
global vcfPositions


def score(pair: List[int]) -> float:
    ##Computing probability that a fair coin generates something as or more
    ##extremely biased than we we are seeing.
    k = min(pair)
    if pair[0] == k:
        var = 0
    else:
        var = 1
    p = 0.5
    n = sum(pair)
    if n > 50:
        return cdf(n, k, var)
    else:
        prob = 0
        for i in range(k + 1):
            if var == 0:
                prob += (p ** i) * ((1 - p) ** (n - i)) * binom(n, i)
            else:
                prob += ((1 - p) ** i) * (p ** (n - i)) * binom(n, i)

        return prob


def rate_tail_score(pair, r):
    r = max(r.values())
    m = min(pair)
    M = max(pair)
    n = sum(pair)
    var = M / float(n) <= r

    prob = 0
    if var:
        if n > 50:
            return cdfp(n, M, r)
        else:
            p = r
            for i in range(M + 1):
                prob += (p ** i) * ((1 - p) ** (n - i)) * binom(n, i)

    else:
        if n > 50:
            return cdfp(n, m, 1 - r)
        else:
            p = 1 - r
            for i in range(m + 1):
                prob += (p ** i) * ((1 - p) ** (n - i)) * binom(n, i)

    return (-1) * prob


def binom(n: int, k: int) -> int:
    # n CHOOSE k
    j = min(k, n - k)
    numerator = 1
    for i in range(j):
        numerator = numerator * (n - i)
    return numerator // math.factorial(j)


def cdf(n: int, k: int, var: int) -> float:
    p = 0.5
    if var == 0:
        mu = n * p
    else:
        mu = n * (1 - p)

    return 0.5 * (1 + math.erf((k + 0.5 - mu) / math.sqrt(2 * n * p * (1 - p))))
