#!/usr/bin/env python3
"""
Compare spin-correlation / beam-polarisation coefficients against the CMS
measurement.

Reports, for each prediction (MadSpin NLO, Toponium Restricted LO,
Toponium Full LO):
  - signed percentage difference w.r.t. CMS, with propagated uncertainty
  - the PULL (z-score), which is the statistically correct comparison

CAVEATS (read before trusting any percentage):
  * Percentage difference is unreliable when the CMS central value is
    consistent with zero (|x_CMS| < ~2 sigma). All six B coefficients and
    C_rr fall in this regime -> the percentages there are noise.
  * Toponium C coefficients are ~ +-1 by construction (1S0 pseudoscalar
    bound-state spin structure). A percent difference vs the continuum-
    dominated CMS value is a category error UNLESS the toponium is first
    weighted by its near-threshold production fraction and added to the
    continuum (admixture interpretation). Standalone, ignore the C percentages.
"""

import math

# (CMS, MadSpin_NLO, Toponium_Restricted_LO, Toponium_Full_LO)
# each entry: (value, uncertainty)
DATA = {
    "B1k": ((0.005, 0.023), (-0.00258, 0.00245), (0.00190, 0.00260), (-0.00039, 0.00251)),
    "B2k": ((0.007, 0.023), (-0.00099, 0.00245), (0.00108, 0.00259), (-0.00472, 0.00251)),
    "B1r": ((-0.023, 0.017), (0.00155, 0.00245), (0.00267, 0.00259), (0.00052, 0.00251)),
    "B2r": ((-0.010, 0.020), (0.00086, 0.00245), (0.00134, 0.00259), (-0.00074, 0.00251)),
    "B1n": ((0.006, 0.013), (-0.00431, 0.00245), (0.00372, 0.00260), (0.00186, 0.00251)),
    "B2n": ((0.017, 0.013), (0.00280, 0.00245), (-0.00016, 0.00259), (-0.00467, 0.00251)),
    "Ckk": ((0.300, 0.038), (0.36143, 0.00421), (1.00302, 0.00424), (1.01218, 0.00410)),
    "Crr": ((0.081, 0.032), (-0.03414, 0.00425), (-0.99822, 0.00423), (-0.99857, 0.00410)),
    "Cnn": ((0.329, 0.020), (0.32749, 0.00421), (1.00091, 0.00423), (1.01242, 0.00409)),
}

PREDICTIONS = ["MadSpin NLO", "Toponium Restricted LO", "Toponium Full LO"]

# A CMS value is treated as "consistent with zero" (percent diff unreliable)
# when |value| < ZERO_SIGMA * uncertainty.
ZERO_SIGMA = 2.0


def percent_diff(pred, cms):
    """Signed % difference (pred - cms)/cms * 100 with propagated uncertainty."""
    xp, sp = pred
    xc, sc = cms
    if xc == 0:
        return float("nan"), float("nan")
    pct = (xp - xc) / xc * 100.0
    # error propagation: f = (xp - xc)/xc ; df/dxp = 1/xc ; df/dxc = -xp/xc^2
    sigma = 100.0 * math.sqrt((sp / xc) ** 2 + (xp * sc / xc**2) ** 2)
    return pct, sigma


def pull(pred, cms):
    """(pred - cms) / sqrt(sig_pred^2 + sig_cms^2). The correct comparison."""
    xp, sp = pred
    xc, sc = cms
    denom = math.sqrt(sp**2 + sc**2)
    if denom == 0:
        return float("nan")
    return (xp - xc) / denom


def cms_consistent_with_zero(cms):
    xc, sc = cms
    return abs(xc) < ZERO_SIGMA * sc


def main():
    header = (
        f"{'Coeff':<6} {'Prediction':<24} "
        f"{'%diff':>12} {'sig(%)':>10} {'pull[sig]':>10}  flag"
    )
    print(header)
    print("-" * len(header))

    for coeff, vals in DATA.items():
        cms = vals[0]
        flag_zero = cms_consistent_with_zero(cms)
        for name, pred in zip(PREDICTIONS, vals[1:]):
            pct, spct = percent_diff(pred, cms)
            z = pull(pred, cms)

            flags = []
            if flag_zero:
                flags.append("PCT-UNRELIABLE(CMS~0)")
            if coeff.startswith("C") and abs(pred[0]) > 0.9:
                flags.append("PURE-TOPONIUM(|C|~1)")
            flag_str = " ".join(flags)

            print(
                f"{coeff:<6} {name:<24} "
                f"{pct:>11.1f}% {spct:>10.1f} {z:>10.2f}  {flag_str}"
            )
        print()

    print("Notes:")
    print("  - 'pull' is the difference in units of combined sigma. |pull| < ~2")
    print("    means the prediction is statistically compatible with CMS.")
    print("  - Rows flagged PCT-UNRELIABLE: CMS value consistent with 0; ignore %.")
    print("  - Rows flagged PURE-TOPONIUM: |C|~1 is the bound-state spin structure,")
    print("    not a meaningful continuum comparison unless weighted by the")
    print("    toponium production fraction and added to the tt continuum.")


if __name__ == "__main__":
    main()