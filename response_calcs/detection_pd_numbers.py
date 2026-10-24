"""Detection numbers quoted in the reply letter (Reviewer 3, Comment 2).

GLRT at the central processor: each of the Nrx receive APs matched-filters the
Ns = Nue + 1 known streams, and the output energies are summed. The normalized
statistic is chi-square with 2*Nrx*Ns degrees of freedom under H0 and
noncentral chi-square with noncentrality 2*N*sum_r gamma_r under H1. The
per-receive-AP SNRs gamma_r are independent and exponentially distributed with
mean gamma_t (random echo phases). gamma_t values are the Fig. 3 numbers.

Also prints the TX-AP-to-target distances, which show why delay separation
cannot justify the non-coherent sum in (2) for the default geometry.

Run: python3 detection_pd_numbers.py  (needs numpy and scipy)
"""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2, ncx2

NRX, NS, PFA = 4, 5, 1e-6
GAMMA_NOAN_DB, GAMMA_PROP_DB = 14.7, 11.1
DF = 2 * NRX * NS
TAU = chi2.isf(PFA, DF)
G = np.random.default_rng(1).exponential(1.0, (400000, NRX)).sum(axis=1)


def pd(gamma_db, n):
    return ncx2.sf(TAU, DF, 2 * n * 10 ** (gamma_db / 10) * G).mean()


def required_snr_db(n, target=0.9):
    return brentq(lambda x: pd(x, n) - target, -60, 40)


tx = [(-200, -200), (200, 200), (-200, 200), (200, -200),
      (0, -200), (0, 200), (-200, 0), (200, 0)]
target = np.array((-75, 75))
dist = sorted(round(float(np.linalg.norm(np.array(p) - target)), 1) for p in tx)
print("TX-AP-to-target distances (m):", dist)

for n in (5, 128):
    req = required_snr_db(n)
    print(f"N = {n:3d}: Pd(no AN, {GAMMA_NOAN_DB} dB) = {pd(GAMMA_NOAN_DB, n):.4f}, "
          f"Pd(proposed, {GAMMA_PROP_DB} dB) = {pd(GAMMA_PROP_DB, n):.4f}, "
          f"Pd = 0.9 needs {req:.1f} dB "
          f"(margins {GAMMA_NOAN_DB - req:.1f} / {GAMMA_PROP_DB - req:.1f} dB)")
print(f"3.6 dB loss -> bistatic distance factor {10 ** (-3.6 / 40):.3f}")
