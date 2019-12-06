from math import log, exp
from common import choose
from read import Read, sample_from_reads


def get_counts(
    snps: list[int], reads: list[Read], shuffle: bool = True
) -> dict[int, list[int]]:
    if shuffle:
        reads = sample_from_reads(reads)
    counts = {i: [0, 0] for i in range(len(snps))}
    back = dict[int, int]()
    for i in range(len(snps)):
        back[snps[i]] = i
    for R in reads:
        for snp in R.keys:
            counts[back[snp]][R.read[snp] % 2] += R.count
    return counts


def approx_rate(counts: list[list[int]]) -> float:
    return float(sum(max(c) for c in counts)) / float(sum(sum(c) for c in counts))


def trans(s1: int, s2: int, rate: float) -> float:
    if s1 == s2:
        return rate
    else:
        return 1 - rate


def forward(
    FD: dict[tuple[int, int], float],
    X: int,
    vec: dict[int, list[int]],
    i: int,
    rate: float,
    p: float,
) -> float:
    if (X, i) in FD:
        return FD[X, i]
    e = 0.0001
    p = p * (1 - e) + (1 - p) * (e)
    q = 1 - p
    if i == 0:
        Y = vec[i]
        if X == 1:
            FD[1, 0] = 0
            return 0
        else:
            val = choose(Y[0] + Y[1], Y[0]) * (p ** Y[0]) * (q ** Y[1])
            FD[0, 0] = 10 ** 10 * val
            return 10 ** 10 * val
    else:
        Y = vec[i]
        f0 = forward(FD, 0, vec, i - 1, rate, p)
        f1 = forward(FD, 1, vec, i - 1, rate, p)
        if (0, i - 1) not in FD:
            FD[0, i - 1] = (10 ** 10) * f0
        if (1, i - 1) not in FD:
            FD[1, i - 1] = (10 ** 10) * f1

        SUM_log = log(trans(0, X, rate) * f0 + trans(1, X, rate) * f1)
        sample_prob_log = (Y[X] * log(p)) + (Y[1 - X] * log(q))
        binom_log = log(choose(sum(Y), Y[0]))
        f_val = exp(SUM_log + sample_prob_log + binom_log)

        if (X, i) not in FD:
            FD[X, i] = 10 ** 10 * f_val
        return 10 ** 10 * f_val


def NW(vec: dict[int, list[int]], rate: float) -> float:
    left = 0.001
    right = 1 - left
    dx = 0.001
    d = 1

    def f(x):
        FD = dict[tuple[int, int], float]()
        return log(
            forward(FD, 0, vec, len(vec) - 1, rate, x)
            + forward(FD, 1, vec, len(vec) - 1, rate, x)
        )

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


def find_rates(snps: list[int], reads: list[Read], rate: float) -> dict[int, float]:
    ##if sufficient coverage, approximation works well
    ##snps must be sorted
    lb_for_approx = 200
    counts = get_counts(snps, reads)
    m = max(max(s) for s in counts.values())

    if m == 0:
        raise ValueError(f"should not happen {snps} {reads}")
    elif m > lb_for_approx:
        r = approx_rate(list(counts.values()))
    else:
        r = NW(counts, rate)

    return {0: r, 1: 1 - r}  ##had swapped