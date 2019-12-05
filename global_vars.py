from math import erf, sqrt, factorial


CUTOFF = 0
QUALITY_CUTOFF = 10
CONFIDENCE = 0.5
PAIR_TRESH = 0.8


V = dict[int, dict[int, str]]()
vcfChroms = list[str]()
vcfPositions = list[str]()


def cdf(n: int, k: int, var: int) -> float:
    p = 0.5
    mu = n * p if var == 0 else n * (1 - p)
    return 0.5 * (1 + erf((k + 0.5 - mu) / sqrt(2 * n * p * (1 - p))))


def binom(n: int, k: int) -> int:
    # n CHOOSE k
    j = min(k, n - k)
    numerator = 1
    for i in range(j):
        numerator = numerator * (n - i)
    return numerator // factorial(j)


def score(pair: list[int]) -> float:
    ##Computing probability that a fair coin generates something as or more
    ##extremely biased than we we are seeing.
    k = min(pair)
    var = 0 if pair[0] == k else 1
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
