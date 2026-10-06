"""Reproduce Fig. 1 from the analytic step-wall resolvent.

Requires mpmath and matplotlib. These are high-precision illustrative
evaluations, not the interval certificate (run verify_certificate.py for that).
No waveform search or finite-dimensional Hamiltonian is used.
"""
from pathlib import Path
import csv
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 55
G = mp.mpf(4)
D = G-1


def resolvent(z):
    a, b = mp.sqrt(G+z), mp.sqrt(z)
    return (a+2+b)/(2*(a+b)*(1+a)**2)


def kernel(s):
    return -mp.sqrt(s)*mp.im(resolvent(1j*mp.sqrt(s)))


def decreasing_root(fun, target):
    lo, hi = mp.mpf("1e-40"), mp.mpf("1e40")
    if not (fun(lo) > target > fun(hi)):
        raise ArithmeticError("Multiplier root not bracketed")
    # Geometric bisection resolves the entire range of multiplier scales.
    for _ in range(165):
        mid = mp.sqrt(lo*hi)
        if fun(mid) > target:
            lo = mid
        else:
            hi = mid
    root = mp.sqrt(lo*hi)
    if abs(fun(root)/target-1) > mp.mpf("1e-35"):
        raise ArithmeticError("Multiplier residual too large")
    return root


def envelopes(action):
    radius2 = action**2/2
    sd = decreasing_root(lambda s: -D**2*mp.diff(resolvent, s), radius2)
    sv = decreasing_root(lambda s: D**2*mp.diff(kernel, s), radius2)
    mean = D**2*(resolvent(0)-resolvent(sd))-sd*radius2
    variance = D**2*kernel(sv)-sv*radius2
    if not (0 < mean < 1 and 0 < variance < mp.mpf("4.5")):
        raise ArithmeticError("Invalid positive relaxed value")
    return mean, variance, sd, sv


def main():
    directory = Path(__file__).resolve().parent
    actions = [mp.power(10, -2+4*mp.mpf(i)/120) for i in range(121)]
    rows = [(a, *envelopes(a)) for a in actions]
    for previous, current in zip(rows, rows[1:]):
        if not (current[1] < previous[1] and current[2] < previous[2]):
            raise ArithmeticError("Lower envelopes must decrease")
    with (directory/"control_bounds.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["action_Lambda", "D", "V", "sigma_D", "sigma_V"])
        writer.writerows([[mp.nstr(x, 35) for x in row] for row in rows])

    plt.rcParams.update({"font.family": "serif", "font.size": 9,
                         "mathtext.fontset": "stix", "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(3.4, 2.45), layout="constrained")
    x = [float(row[0]) for row in rows]
    ax.loglog(x, [float(row[1]) for row in rows], color="#155b82", lw=1.7,
              label=r"Dissipation: $D(\Lambda/\sqrt{2})$")
    ax.loglog(x, [float(row[2]/mp.mpf("4.5")) for row in rows],
              color="#b74c20", lw=1.7, ls=(0, (4, 2)),
              label=r"Variance: $\mathcal{V}(\Lambda/\sqrt{2})/(9/2)$")
    ax.set(xlabel=r"Integrated height action $\Lambda$",
           ylabel="Normalized lower bound", xlim=(.01, 100),
           ylim=(1e-17, 2))
    ax.set_yticks([1, 1e-4, 1e-8, 1e-12, 1e-16])
    ax.grid(axis="both", which="major", color="0.9", lw=.5)
    ax.legend(loc="lower left", fontsize=7.2, frameon=False)
    ax.text(.97, .93, r"Step wall, $g=4$", transform=ax.transAxes,
            ha="right", va="top", fontsize=8)
    fig.savefig(directory/"control_bounds.png", dpi=600)
    print("Figure and 121-row data written; positive monotone envelopes and multiplier residuals checked.")


if __name__ == "__main__":
    main()
