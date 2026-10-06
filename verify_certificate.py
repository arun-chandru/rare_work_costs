"""Recompute the continuum final-band and dissipated-work certificate.

Requires Python 3 and the installed package python-flint. No other package,
input data file, workspace path, eigenvalue solver, quadrature or time
propagation is used. Floating arithmetic proposes root brackets only; every
accepted root, spectral membership, overlap and final bound is checked by
FLINT/Arb at 256 bits with exact rational inputs.

Usage:
    python verify_certificate.py
    python verify_certificate.py --output certificate.json

All verification gates raise exceptions and remain active under python -O.
With --output, the complete interval artifact is saved and a short JSON
summary is printed. Without --output, the complete artifact is printed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import time

try:
    import flint
    from flint import arb, ctx, fmpq
except ImportError as exc:
    raise SystemExit("This verifier requires the installed package python-flint.") from exc

ctx.prec = 256
INITIAL_BARRIER = 625
FINAL_BARRIER = 2500
THETA = Fraction(1, 10)
ACTION = Fraction(3, 8)
PROBABILITY_FLOOR = Fraction(293544, 10**10)
EQUILIBRIUM_UPPER = Fraction(65, 10**54)
WORK_FLOOR = Fraction(319, 10000)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def ball(value):
    """Enclose an exact rational; disallow accidental floating inputs."""
    if isinstance(value, arb):
        return value
    if isinstance(value, float):
        raise TypeError("Certified arithmetic requires exact rational inputs.")
    value = Fraction(value)
    return arb(fmpq(value.numerator, value.denominator))


def positive(value, message):
    """Accept a strict inequality only when the whole ball proves it."""
    require(value > 0, f"{message}: {value}")


def hull(lower, upper):
    require(lower <= upper, "Reversed rational root interval.")
    return ball(lower).union(ball(upper))


def subbarrier_secular(k, barrier):
    a = (ball(barrier) - k*k).sqrt()
    return k*k.cos() + a*a.coth()*k.sin()


def propose_subbarrier(n, barrier):
    """Noncertifying candidate search; the returned interval is then checked."""
    lower = (n-0.5)*math.pi + 1e-8
    upper = min(n*math.pi-1e-8, math.sqrt(barrier)-1e-8)

    def value(k):
        a = math.sqrt(barrier-k*k)
        return k*math.cos(k) + a*math.sin(k)/math.tanh(a)

    f_lower = value(lower)
    require(f_lower*value(upper) < 0, "Subbarrier candidate search failed.")
    for _ in range(60):
        middle = (lower+upper)/2
        if f_lower*value(middle) > 0:
            lower, f_lower = middle, value(middle)
        else:
            upper = middle
    center = round((lower+upper)/2 * 10**8)
    return Fraction(center-2, 10**8), Fraction(center+2, 10**8)


def bisect_certified(lower, upper, secular):
    """Exact-rational bisection with interval-certified signs at every step."""
    f_lower, f_upper = secular(ball(lower)), secular(ball(upper))
    positive(-f_lower*f_upper, "Initial opposite secular signs")
    left_sign = 1 if f_lower > 0 else -1
    for _ in range(80):
        middle = (lower+upper)/2
        f_middle = secular(ball(middle))
        if f_middle > 0:
            middle_sign = 1
        elif f_middle < 0:
            middle_sign = -1
        else:
            raise ArithmeticError(f"Unresolved root sign at {middle}: {f_middle}")
        if middle_sign == left_sign:
            lower = middle
        else:
            upper = middle
    positive(-secular(ball(lower))*secular(ball(upper)),
             "Final opposite secular signs")
    require(upper-lower < Fraction(4, 10**32), "Root interval too wide.")
    return lower, upper


def subbarrier_mode(n, barrier):
    lower, upper = propose_subbarrier(n, barrier)
    positive(ball(barrier)-ball(upper)**2, "Candidate below barrier")
    lower, upper = bisect_certified(
        lower, upper, lambda k: subbarrier_secular(k, barrier))
    k = hull(lower, upper)
    positive(k-(ball(n)-ball(Fraction(1, 2)))*arb.pi(), "Left nodal bound")
    positive(ball(n)*arb.pi()-k, "Right nodal bound")
    positive(ball(barrier)-k*k, "Entire root interval below barrier")
    a = (ball(barrier)-k*k).sqrt()
    j_minus = a.coth()/(2*a)-1/(2*a.sinh()**2)
    exterior = k.sin()**2*j_minus
    norm = ball(Fraction(1, 2))-(2*k).sin()/(4*k)+exterior
    positive(norm, "Positive subbarrier normalization")
    return dict(label=n, kind="subbarrier", barrier=barrier,
                lower=lower, upper=upper, k=k, a=a, energy=k*k,
                norm=norm, source=exterior/norm)


def superbarrier_secular(k):
    b = (k*k-FINAL_BARRIER).sqrt()
    return k*k.cos()*b.sin()+b*k.sin()*b.cos()


def propose_superbarrier():
    """Locate a useful subset, not a completeness test of the final spectrum."""
    def value(k):
        b = math.sqrt(k*k-FINAL_BARRIER)
        return k*math.cos(k)*math.sin(b)+b*math.sin(k)*math.cos(b)

    candidates = []
    lower = 50.000001
    f_lower = value(lower)
    end = math.sqrt(INITIAL_BARRIER*7.9)
    while lower < end:
        upper = min(lower+0.001, end)
        f_upper = value(upper)
        if f_lower*f_upper < 0:
            left, right, f_left = lower, upper, f_lower
            for _ in range(50):
                middle = (left+right)/2
                f_middle = value(middle)
                if f_left*f_middle > 0:
                    left, f_left = middle, f_middle
                else:
                    right = middle
            center = round((left+right)/2 * 10**8)
            candidates.append((Fraction(center-2, 10**8),
                               Fraction(center+2, 10**8)))
        lower, f_lower = upper, f_upper
    return candidates


def superbarrier_modes():
    modes = []
    for index, (lower, upper) in enumerate(propose_superbarrier(), 1):
        positive(ball(lower)**2-FINAL_BARRIER, "Candidate above final barrier")
        lower, upper = bisect_certified(lower, upper, superbarrier_secular)
        if modes:
            require(modes[-1]["upper"] < lower, "Final root intervals overlap.")
        k = hull(lower, upper)
        positive(k*k-FINAL_BARRIER, "Entire root interval above barrier")
        b = (k*k-FINAL_BARRIER).sqrt()
        positive(abs(b.sin()), "Nonzero exterior sine denominator")
        j_plus = (ball(Fraction(1, 2))-(2*b).sin()/(4*b))/b.sin()**2
        exterior = k.sin()**2*j_plus
        norm = ball(Fraction(1, 2))-(2*k).sin()/(4*k)+exterior
        positive(norm, "Positive superbarrier normalization")
        modes.append(dict(label=None, subset_index=index, kind="superbarrier",
                          barrier=FINAL_BARRIER, lower=lower, upper=upper,
                          k=k, b=b, energy=k*k, norm=norm, source=exterior/norm))
    require(len(modes) == 22, "Expected 22 certified superbarrier subset roots.")
    return modes


def overlap_square(initial, final):
    a = initial["a"]
    if final["kind"] == "subbarrier":
        b = final["a"]
        positive(abs(a-b), "Nonzero hyperbolic-product denominator")
        integral = ((a+b).sinh()/(2*(a+b))-(a-b).sinh()/(2*(a-b)))
        integral /= a.sinh()*b.sinh()
    else:
        b = final["b"]
        integral = (a*a.coth()-b*b.cot())/(a*a+b*b)
    cross = initial["k"].sin()*final["k"].sin()*integral
    gap = final["energy"]-initial["energy"]
    positive(gap, "Positive endpoint energy gap")
    return ((ball(FINAL_BARRIER-INITIAL_BARRIER)*cross/gap)**2
            /(initial["norm"]*final["norm"]))


def serial_mode(mode):
    secular = ((lambda k: subbarrier_secular(k, mode["barrier"]))
               if mode["kind"] == "subbarrier" else superbarrier_secular)
    return {
        "ordered_label": mode["label"], "subset_index": mode.get("subset_index"),
        "kind": mode["kind"], "barrier": mode["barrier"],
        "k_lower_exact": str(mode["lower"]), "k_upper_exact": str(mode["upper"]),
        "energy_ball": str(mode["energy"]),
        "scaled_energy_ball": str(mode["energy"]/INITIAL_BARRIER),
        "secular_lower_ball": str(secular(ball(mode["lower"]))),
        "secular_upper_ball": str(secular(ball(mode["upper"]))),
        "normalization_integral_ball": str(mode["norm"]),
        "exterior_probability_ball": str(mode["source"]),
    }


def binary_relative_entropy(p, q):
    p, q = ball(p), ball(q)
    return p*(p/q).log()+(1-p)*((1-p)/(1-q)).log()


def verify():
    initial = [subbarrier_mode(n, INITIAL_BARRIER) for n in range(1, 9)]
    final = ([subbarrier_mode(n, FINAL_BARRIER) for n in range(12, 17)]
             + superbarrier_modes())
    require(len(initial) == 8 and len(final) == 27, "Incorrect root counts.")
    for previous, current in zip(final[:-1], final[1:]):
        require(previous["upper"] < current["lower"],
                "Selected final root intervals are not disjoint.")
    for mode in final:
        positive(mode["energy"]/INITIAL_BARRIER-ball(Fraction(21, 10)),
                 "Strict lower interior spectral membership")
        positive(ball(Fraction(79, 10))-mode["energy"]/INITIAL_BARRIER,
                 "Strict upper interior spectral membership")

    theta = ball(THETA)
    weights = [(-theta*mode["energy"]).exp() for mode in initial]
    partial_z = sum(weights, ball(0))
    c = theta*arb.pi()**2/4
    tail_denominator = 1-(-19*c).exp()
    positive(tail_denominator, "Positive geometric-tail denominator")
    remainder = (-81*c).exp()/tail_denominator
    z_lower = partial_z.lower()
    z_upper = (partial_z+remainder).upper()
    positive(z_lower, "Positive Gibbs partition lower bound")
    overlaps = [[overlap_square(i, f) for f in final] for i in initial]
    projected = [sum(row, ball(0)) for row in overlaps]
    static_lower = (sum((w*a for w, a in zip(weights, projected)), ball(0))
                    / z_upper).lower()
    source_upper = ((sum((w*m["source"] for w, m in zip(weights, initial)), ball(0))
                     + remainder)/z_lower).upper()
    positive(static_lower, "Positive static band probability")
    amplitude = static_lower.sqrt()-ball(ACTION)*source_upper.sqrt()
    positive(amplitude, "Positive action-constrained band amplitude")
    probability_floor = amplitude**2
    positive(static_lower-ball(Fraction(22980717, 10**11)),
             "Static band probability exceeds 0.00022980717")
    positive(ball(Fraction(674811, 10**9))-source_upper,
             "Source probability below 0.000674811")
    positive(probability_floor-ball(PROBABILITY_FLOOR),
             "Driven band probability exceeds 0.0000293544")

    threshold = Fraction(1250)
    heat_time = 1/(2*threshold)
    positive(ball(THETA-heat_time), "Heat time strictly below inverse temperature")
    equilibrium_analytic = (2*ball(threshold)/arb.pi()).sqrt() * (
        ball(Fraction(1, 2)-THETA*threshold)+theta*arb.pi()**2).exp()
    positive(equilibrium_analytic, "Positive equilibrium probability bound")
    positive(ball(EQUILIBRIUM_UPPER)-equilibrium_analytic,
             "Final Gibbs band probability below 6.5e-53")
    positive(ball(PROBABILITY_FLOOR-EQUILIBRIUM_UPPER),
             "Probability floor exceeds equilibrium upper bound")
    positive(ball(1-PROBABILITY_FLOOR), "Probability floor below one")
    work_analytic = binary_relative_entropy(PROBABILITY_FLOOR,
                                           equilibrium_analytic)/theta
    work_rational_q = binary_relative_entropy(PROBABILITY_FLOOR,
                                             EQUILIBRIUM_UPPER)/theta
    positive(work_rational_q-ball(WORK_FLOOR), "Dissipated work exceeds 0.0319")
    positive(work_analytic-work_rational_q, "Rational relaxation weakens witness")
    p, q = ball(PROBABILITY_FLOOR), ball(EQUILIBRIUM_UPPER)
    entropy = -p*p.log()-(1-p)*(1-p).log()
    simpler_work = (p*(1/q).log()-entropy)/theta
    positive(simpler_work-ball(WORK_FLOOR), "Entropy-only witness exceeds 0.0319")
    positive((p*(1-q)/(q*(1-p))).log(), "Binary entropy increases in p")
    positive((p-q)/(q*(1-q)), "Binary entropy decreases in q")

    return {
        "status": "all continuum and entropic certificate checks passed",
        "parameters": {
            "interval": [0, 2], "dirichlet_endpoints": [0, 2], "step_location": 1,
            "r_exact": "1/25", "theta_exact": str(THETA),
            "initial_barrier": INITIAL_BARRIER, "final_barrier": FINAL_BARRIER,
            "action_upper_exact": str(ACTION), "scaled_final_band": [2, 8],
            "unscaled_final_band": [1250, 5000],
            "selected_scaled_energies_strictly_inside": ["21/10", "79/10"],
        },
        "method": {
            "preparation": "complete initial Gibbs state; no postselection",
            "continuum_roots_recomputed": True, "input_artifacts": [],
            "selected_final_subset_complete_band_claim": False,
            "matrix_discretization": False, "time_propagation": False,
            "quadrature": False, "classical_comparison": False,
            "detector_or_experimental_claim": False,
            "precision_bits": ctx.prec, "certified_bisections_per_root": 80,
            "initial_root_count": len(initial), "selected_final_root_count": len(final),
            "python_version": platform.python_version(),
            "python_flint_version": getattr(flint, "__version__", "unavailable"),
            "verification_uses_python_assert_statements": False,
        },
        "initial_modes": [serial_mode(mode) for mode in initial],
        "selected_final_modes": [serial_mode(mode) for mode in final],
        "initial_gibbs_weight_balls": [str(w) for w in weights],
        "overlap_square_balls_initial_by_selected_final":
            [[str(value) for value in row] for row in overlaps],
        "initial_projected_probability_balls": [str(value) for value in projected],
        "unnormalized_omitted_gibbs_upper_ball": str(remainder),
        "partition_lower": str(z_lower), "partition_upper": str(z_upper),
        "static_band_probability_lower": str(static_lower),
        "source_probability_upper": str(source_upper),
        "action_band_amplitude_lower_ball": str(amplitude),
        "uniform_band_probability_lower_ball": str(probability_floor),
        "final_gibbs_probability_analytic_upper_ball": str(equilibrium_analytic),
        "heat_time_exact": str(heat_time),
        "dissipated_work_witness_analytic_q_ball": str(work_analytic),
        "dissipated_work_witness_rational_q_ball": str(work_rational_q),
        "entropy_only_work_witness_ball": str(simpler_work),
        "certified_rational_bounds": {
            "static_band_probability_lower_strict": "22980717/100000000000",
            "source_probability_upper_strict": "674811/1000000000",
            "driven_band_probability_lower_strict": str(PROBABILITY_FLOOR),
            "final_gibbs_probability_upper_strict": str(EQUILIBRIUM_UPPER),
            "dissipated_work_lower_strict": str(WORK_FLOOR),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Save the complete JSON artifact.")
    args = parser.parse_args()
    start = time.perf_counter()
    result = verify()
    result["elapsed_seconds"] = time.perf_counter()-start
    result["verifier_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    payload = json.dumps(result, indent=2)+"\n"
    if args.output is not None:
        args.output.write_text(payload, encoding="utf-8")
        print(json.dumps({
            "status": result["status"], "output": str(args.output),
            "elapsed_seconds": result["elapsed_seconds"],
            "certified_rational_bounds": result["certified_rational_bounds"],
            "root_counts": [result["method"]["initial_root_count"],
                            result["method"]["selected_final_root_count"]],
        }, indent=2))
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
