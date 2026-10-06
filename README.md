# Rare work costs under bounded-action quantum confinement control

**R. Arun Chandru | 6 October 2026**

## What the paper says

A particle is confined by a high but finite potential wall. The wall's height is changed, and the particle's initial state is a complete thermal mixture. Almost all particles retain their ordered energy-level label as the wall becomes harder. That familiar measure of successful control nevertheless misses a small population at increasingly high energies. These rare outcomes dominate the leading dissipated work and make the work variance grow.

The paper derives this rare-energy law from the actual continuum Hamiltonian, without discarding the initially occupied high-energy states or replacing the finite wall by a perfectly reflecting boundary. It then proves a limitation on pulse shaping: changing only the wall height with a bounded integrated height action cannot remove the leading energetic cost. This remains true when the waveform and duration change with the wall hardness. The mechanism is an endpoint displacement extending over an entire half-line in the boundary coordinate; the finite-norm correction produced by bounded-action control cannot cancel it.

The leading thermodynamic information is also accessible in principle through final-energy measurements. For each fixed protocol, a finite set of appropriately chosen high-energy bands recovers arbitrarily much of the leading dissipative relative entropy. The energy-basis coherence contribution is smaller order. This does not make the measurement inexpensive: the energy range increases and the counts become rare. A separate finite-height continuum certificate supplies the probability floor and work bound, uniformly over the stated action class.


## Scope and quantifiers

| Statement | Precise scope |
|---|---|
| Leading label-changing probability and positive energetic lower bounds | Uniform over all positive absolutely continuous height paths and finite durations with the common integrated action bound, including hardness-dependent choices |
| Exact rare-event measure and complete leading work coefficients | Fixed path, finite scaled duration, and fixed spatial profile before the hardening limit |
| Finite-bin recovery of leading dissipative entropy | Bins may depend on that fixed protocol, profile, and desired accuracy; chosen before the limit |
| Smooth-wall extension | Each fixed nonnegative smooth transition profile with a constant right plateau; physical transition width shrinks with the penetration length |
| Finite certified numbers | The specified sharp-wall tuple, not a rounded-wall or detector calibration |

The mean dissipated work tends to zero in absolute units at the derived order; its coefficient has a positive lower bound. The work variance diverges at a different order. Dissipated work is not heat exchanged with a bath during the isolated protocol. No general many-body, auxiliary-control, fixed-width-apparatus, or universal speed-limit claim is made. Detector errors and finite-sample confidence intervals are not included.

## Prior art and the contribution

The manuscript cites established two-measurement work statistics, work tails and ideal-wall quenches, weak-driving and fast-driving control, bilinear-control estimates, quadratic norm-ball minimization, quantum entropy monotonicity, and the population/coherence decomposition. These ingredients are not claimed as inventions here.

The candidate contribution is their concrete continuum realization and closure in the stated singular family: the complete energy-weighted boundary-source law; finite-system proofs of action-uniform positive lower bounds; a fixed-profile extension; finite-bin recovery of the leading dissipative entropy; and a reproducible finite-height probability/work certificate. The closest comparisons are:

| Prior work | Distinction retained here |
|---|---|
| Smacchia-Silva local quenches; Jarzynski-Quan-Rahav work tails; Aslangul and Quan-Jarzynski ideal-wall/piston results | These precede the general rare-event lesson. This work retains a finite penetrable wall and derives its thermal source law and control-uniform lower bounds. |
| Guarnieri-Eisert-Miller generalized linear response; Rolandi-Perarnau-Llobet-Miller rapid-driving control; Xia-Nazir-Miller boundary-conformal speed limit | Different perturbative regimes, resource constraints, and models. No claim that those results fail. The endpoint amplitude here is fixed and finite. |
| Boussaid-Caponigro-Chambrion control estimates; More-Sorensen quadratic minimization | General control bounds and the minimization method are established. The source spectrum, hardening limit, finite-system certificates, and constants are model-specific. |
| Deffner-Lutz; Lindblad; Francica-Goold-Plastina; Rubino-Brukner-Manzano | Thermodynamic relative entropy, data processing, coherence decomposition, and coarse-graining distinctions are prior art. The leading-order finite-bin recovery is derived for the actual rare-energy law here. |

## Files

| File | Purpose |
|---|---|
| `r_arun_chandru_rare_work_costs.pdf` | Combined Letter, bibliography, and long Supplemental Material |
| `r_arun_chandru_rare_work_costs.tex` | Complete source with inline bibliography; no private or external TeX inputs |
| `control_bounds.png` | The single figure, at 600 dpi |
| `control_bounds.csv` | Reproducible 121-row figure data |
| `make_figure.py` | Recomputes the analytic relaxation curves and figure |
| `verify_certificate.py` | Fresh continuum root, overlap, Gibbs-tail, and relative-entropy verifier |
| `certificate.json` | Saved successful verification output, including exact rational root brackets and script hash |
| `requirements.txt` | Direct Python dependencies |
| `README.md` | Scope, build instructions, reproducibility, and prior-art distinctions |

The PDF contains all analytic arguments. Code and saved output support the finite numerical certificate; they do not replace missing analytic proofs.

## Compile the single source

Keep the TeX and figure together. With TeX Live or MiKTeX including REVTeX 4.2:

```sh
pdflatex r_arun_chandru_rare_work_costs.tex
pdflatex r_arun_chandru_rare_work_costs.tex
```

Or use Tectonic, which reruns cross-references automatically:

```sh
tectonic r_arun_chandru_rare_work_costs.tex
```

No separate BibTeX command or `.bib` file is needed. Required LaTeX packages are REVTeX 4.2, amsmath, amssymb, bm, graphicx, and hyperref. A first compilation may download standard TeX packages. Python is **not** needed to compile the supplied manuscript or to read its PDF. Tectonic may report missing database entries from REVTeX's automatically generated notes bibliography; the complete inline bibliography supplies all 19 citations. Some underfull-line and float-placement warnings are nonblocking.

## Reproduce the finite certificate

Use a current Python 3 environment; the supplied run used Python 3.12 and python-flint 0.9.0. A virtual environment is recommended.

```sh
python -m pip install -r requirements.txt
python verify_certificate.py --output certificate_recomputed.json
python -O verify_certificate.py --output certificate_optimized.json
```

The verifier needs only `python-flint`, not the plotting dependencies. Each command recomputes its inputs; it never reads the saved `certificate.json`. It encloses 8 initial and 27 selected final continuum eigenvalues, computes all 216 retained squared overlaps analytically, bounds the entire omitted initial thermal tail, and checks the final rational inequalities. All verification gates raise exceptions on failure; none relies on removable Python `assert` statements. Ordinary and optimized-mode executions passed.

For the sharp wall on `(0,2)`, initial height 625, final height 2500, inverse temperature 0.1, and integrated height action at most 3/8, the verified conclusions include

- `Pr(1250 <= final energy <= 5000) > 293544/10^10 = 0.0000293544`;
- the final Gibbs probability of that band is below `65/10^54`;
- the dissipated work exceeds `319/10000 = 0.0319`, in the paper's energy units.

The selected final modes need not exhaust the band. Omitting further final modes can only weaken this lower bound. The initial preparation is not truncated or postselected. The proof combines Sturm oscillation, exact overlap formulas, min-max thermal bounds, Duhamel's formula, and entropy monotonicity. Interval arithmetic uses 256-bit Arb balls and exact rational root brackets. Its trust assumptions include the libraries, interpreter, and hardware; this is not a proof-assistant formalization.

JSON timing and platform fields can differ between runs. The strict rational conclusions are the reproducibility target, not byte-for-byte equality of output files.

## Reproduce the figure

```sh
python make_figure.py
```

This overwrites `control_bounds.csv` and `control_bounds.png` with recomputed versions. The calculation uses the exact step-wall resolvent at endpoint ratio 4, 55-decimal working precision, and scalar multiplier equations. It checks root residuals, positivity, and monotonicity. It is a high-precision illustration, **not** the interval-arithmetic certificate. It does not optimize physical control waveforms or claim the relaxed minimizing vectors are reachable. Recompile the TeX after regenerating the figure.