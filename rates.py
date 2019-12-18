from math import log, exp, pow, erf, sqrt
from read import Read, sample_from_reads
from typing import Tuple, Dict, List, Set


def binom(n: float, k: float) -> float:
    if k > n:
        raise ValueError('binom k > n')
    if k == 0:
        return 1.0
    if k > n / 2:
        return binom(n, n - k)
    return n * binom(n - 1, k - 1) / k


def get_counts(
    snps: List[int], reads: List[Read], shuffle: bool = True
) -> Dict[int, List[int]]:
    if shuffle:
        reads = sample_from_reads(reads)
    counts = {i: [0, 0] for i in range(len(snps))}
    back: Dict[int, int] = {}
    for i in range(len(snps)):
        back[snps[i]] = i
    for R in reads:
        for snp in R.snps:
            counts[back[snp]][R.snps[snp] % 2] += R.count  # originally just R.read[snp]
    return counts


def approx_rate(counts: List[List[int]]) -> float:
    return float(sum(max(c) for c in counts)) / float(sum(sum(c) for c in counts))


def trans(s1: int, s2: int, rate: float) -> float:
    if s1 == s2:
        return rate
    else:
        return 1 - rate


def forward(
    FD: Dict[Tuple[int, int], float],
    X: int,
    vec: Dict[int, List[int]],
    i: int,
    rate: float,
    p: float
) -> float:
    if (X, i) in FD:
        return FD[X, i]
    e = 0.0001
    p = p * (1 - e) + (1 - p) * (e)
    q = 1 - p
    Y = [float(y) for y in vec[i]]
    pow10 = pow(10.0, 10.0)
    if i == 0:
        if X == 1:
            FD[1, 0] = 0.0
            return 0.0
        else:
            val = binom(Y[0] + Y[1], Y[0]) * pow(p, Y[0]) * pow(q, Y[1]) * pow10
            FD[0, 0] = val
            return val
    else:
        f0 = forward(FD, 0, vec, i - 1, rate, p)
        f1 = forward(FD, 1, vec, i - 1, rate, p)
        if (0, i - 1) not in FD:
            FD[0, i - 1] = pow10 * f0
        if (1, i - 1) not in FD:
            FD[1, i - 1] = pow10 * f1

        SUM_log = log(trans(0, X, rate) * f0 + trans(1, X, rate) * f1)
        sample_prob_log = (Y[X] * log(p)) + (Y[1 - X] * log(q))
        binom_log = log(binom(sum(Y), Y[0]))
        val = pow10 * exp(SUM_log + sample_prob_log + binom_log)
        if (X, i) not in FD:
            FD[X, i] = val
        return val


def NW(vec: Dict[int, List[int]], rate: float) -> float:
    left = 0.001
    right = 1 - left
    dx = 0.001
    d = 1.0

    def fp(vec: Dict[int, List[int]], rate: float, x: float):
        FD: Dict[Tuple[int, int], float] = {}
        return log(
            forward(FD, 0, vec, len(vec) - 1, rate, x)
            + forward(FD, 1, vec, len(vec) - 1, rate, x)
        )
    f = lambda x: fp(vec, rate, x)

    def der(f, x: float, dx: float) -> float:
        return (f(x + dx) - f(x - dx)) / (2 * float(dx))

    if f(0.501) == f(0.499):
        if f(0.5) >= f(0.499):
            return 0.5
        else:
            left = 0.5 + 0.001
    mid = (right + left) / 2.0
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


def find_rates(snps: List[int], reads: List[Read], rate: float) -> List[float]:
    """
    if sufficient coverage, approximation works well
    snps must be sorted
    """

    lb_for_approx = 200
    counts = get_counts(snps, reads)
    m = max(max(s) for s in counts.values())

    if m == 0:
        raise ValueError(f"[find_rates] m cannot be zero ({snps}; {reads})")
    r = 0.0
    if m > lb_for_approx:
        r = approx_rate(list(counts.values()))
    else:
        r = NW(counts, rate)

    return [r, 1 - r]


def cdf(n: float, k: int, var: int) -> float:
    p = 0.5
    mu = n * p if var == 0 else n * (1 - p)
    return 0.5 * (1 + erf((k + 0.5 - mu) / sqrt(2 * n * p * (1 - p))))


def score(pair: List[int]) -> float:
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
