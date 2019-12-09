from math import erf, sqrt, factorial


CUTOFF = 0
QUALITY_CUTOFF = 10
CONFIDENCE = 0.5
PAIR_TRESH = 0.8

DOT = -10 # "."
V = dict[int, dict[int, int]]()


def cdf(n: float, k: int, var: int) -> float:
    p = 0.5
    mu = n * p if var == 0 else n * (1 - p)
    return 0.5 * (1 + erf((k + 0.5 - mu) / sqrt(2 * n * p * (1 - p))))


def binom(n: float, k: float) -> float:
    if k > n:
        raise ValueError('binom k > n')
    if k == 0:
        return 1.0
    if k > n / 2:
        return binom(n, n - k)
    return n * binom(n - 1, k - 1) / k


def score(pair: list[int]) -> float:
    """
    Computing probability that a fair coin generates something as or more
    extremely biased than we we are seeing.
    """
    
    k = min(pair)
    var = 0 if pair[0] == k else 1
    p = 0.5
    n = float(sum(pair))
    if n > 50:
        return cdf(n, k, var)
    else:
        prob = 0.0
        for i in range(k + 1):
            if var == 0:
                prob += (p ** i) * ((1 - p) ** (n - i)) * binom(n, float(i))
            else:
                prob += ((1 - p) ** i) * (p ** (n - i)) * binom(n, float(i))
        return prob
