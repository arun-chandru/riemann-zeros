# Clipped Gram energies and zero counts of the Riemann zeta function

R. Arun Chandru  
arun.chandru@gmail.com  

## Extended abstract

The Riemann zeta function has zeros in a vertical strip of the complex plane. Two basic questions ask which zeros lie on the critical line, where the real part is one half, and which are simple, rather than repeated. This paper proves quantitative lower bounds for several distinguished counts. All statements are asymptotic: they apply with any fixed positive tolerance once the height is sufficiently large. They are not guarantees for every finite height, and no effective starting height is supplied.

The main result is the unconditional lower limiting proportion

\[
\liminf_{T\to\infty}\frac{N_{\mathrm{s0}}(T)}{N(T)}
\ge \frac{66812491}{99194876}
=0.673547805029767868\ldots.
\]

Here the numerator counts simple critical-line zeros, and the denominator counts all nontrivial zeros up to height T with multiplicity. Thus the proportion is at least **67.35478050%**, rounded down. No Riemann Hypothesis, or conjectural higher-moment estimate is assumed.

The mechanism is an energy budget. Repeated zeros consume more quadratic energy than simple zeros, but nearby distinct zeros also interact. An unconditional pair-correlation formula determines the available energy. The difficulty is that zeros away from the critical line lead to a signed Hermitian operator, so the surviving Gram energy is clipped above a spectral threshold. One cannot simply use every positive pair contribution as though clipping had not occurred. An explicit Fourier majorant controls every separated block uniformly in its size. Disjoint close pairs pay for their own removal. This transfers an eight-point inequality to the full clipped energy and gives the displayed fraction. The local coefficients and the general lossless-transfer principle are attributed prior inputs; the numerical application is distinguished from the structural results.

The paper also studies information that scalar proportions discard. It proves that local slack and unused long-range pair energy control how far adjacent gaps are from integral spacing, even when physical spans are unbounded. It then derives a common residual for simultaneous affine deductions using several profiles and clipping thresholds. An explicit two-lattice minimizer and a twelve-fibre analysis describe limitations of that defined deduction method on the full normalized L2 profile space. The resulting method ceilings are not upper bounds on the actual zeta-zero proportions. Finally, exact multiplicity accounting yields further count bounds and explains a sharp scalar limitation of changing the clipping level alone.

## Numerical results and denominators

Every entry below is a lower bound for the liminf of the indicated ratio as T tends to infinity. Terminating percentages are rounded down to eight decimal places. Different rows use the specific profiles and local inequalities stated in the manuscript; their improvements are not added across incompatible profiles.

| Numerator | Denominator | Lower percentage |
|---|---|---:|
| Simple critical-line zeros | All zeros, with multiplicity | 67.35478050% |
| Distinct complex zero locations | All zeros, with multiplicity | 83.90046180% |
| Simple zeros anywhere in the strip | All zeros, with multiplicity | 67.80092360% |
| Simple-zero count plus critical-line multiplicity count | Twice the total multiplicity count | 83.67739025% |
| Multiplicity mass of zeros that are simple or on the critical line | All zeros, with multiplicity | 88.93165002% |
| Mass at locations of multiplicity at most two | All zeros, with multiplicity | 83.79775462% |
| Locations of multiplicity at most two | All distinct complex locations | 93.94524765% |

The average in the fourth row is not a separate lower bound of that size for both of its components. Distinct locations at the same imaginary part are still distinct locations. The fifth row counts the overlap only once.

The exact main fraction is `66812491 / 99194876`. The exact average fraction is `166007367 / 198389752`. The manuscript gives the other rational values and the expression for the union count.

## Contents

- `r_arun_chandru_clipped_gram_energies_zero_counts.pdf`: complete manuscript, including appendices and bibliography.
- `r_arun_chandru_clipped_gram_energies_zero_counts.tex`: the single self-contained LaTeX source. No external figure, bibliography, or custom style file is needed.
- `r_arun_chandru_proof.lean`: the unchanged formal proof source for the main numerical chain, with inherited copyright notices preserved.
- `reproducibility/`: compact exact-arithmetic checks and a runner.
- `BUILD_LEAN.md`, `lean_support/CandidateSpec.lean`, and `lean-toolchain`: pinned formal-build instructions and the small required generated specification.
- `LICENSE`, `THIRD_PARTY_NOTICES.txt`, and `.gitattributes`: software licensing, attribution and line-ending preservation.

## Compile the manuscript

The uploaded PDF has 37 pages and was built with Tectonic 0.17.0.

Use pdfLaTeX twice, or latexmk:

```sh
pdflatex -interaction=nonstopmode -halt-on-error r_arun_chandru_clipped_gram_energies_zero_counts.tex
pdflatex -interaction=nonstopmode -halt-on-error r_arun_chandru_clipped_gram_energies_zero_counts.tex

# Alternative:
latexmk -pdf r_arun_chandru_clipped_gram_energies_zero_counts.tex
```


## Run the small exact checks

Python 3 and Node.js are needed. The exact check modes use only their standard libraries; there is no package installation, network access or numerical sampling used for acceptance.

```sh
python reproducibility/run_all.py

# If Node is not on PATH:
python reproducibility/run_all.py --node /absolute/path/to/node
```

Run with Python assertions enabled: do not use `-O`, `-OO`, or a nonzero `PYTHONOPTIMIZE`. The runner checks this. The twelve-fibre script must be invoked with `--certify`; its default exploratory mode is not part of the exact verification and is not required. The runner selects the correct exact mode automatically.

The eight processes cover:

1. The main cosine profile, rational trigonometric enclosures, majorant and close-pair guards, lag capacities, energy caps, all seven count percentages, and scalar saturation identities: **100 exact checks**.
2. The phase-polynomial certificate and finite-cover constants.
3. The signed Taylor pressure row and strengthened affine coefficients.
4. An independently implemented phase/pressure check.
5. The two-lattice coercivity, positivity and union witness.
6. Independent checks of the printed structural scalar inequalities.
7. The twelve-fibre polynomial, matrix and residual certificate.
8. An independently implemented twelve-fibre check.

These programs check the displayed finite arithmetic. They do not prove Fourier analysis, the analytic pair-correlation theorem, or general variational arguments; those are proved or cited in the manuscript. They also do not rerun the large imported multivariate local-certificate searches. The principal eight-point floor's full finite proof data are in the formal source. The other local inputs have their own public source references and stated computer-assisted verification boundaries.

## Formal coverage

The supplied Lean file is 1,858,713 UTF-8 bytes. Its instantiated internal argument proves the dyadic simple-critical-line bound at the displayed fraction. Its final exported wrappers state distinct-critical-line bounds in dyadic and cumulative height conventions. The manuscript supplies the direct cumulative simple-line argument. The stronger separate counting rows, phase comparison and variational results should not be mistaken for statements all exported by this one file.

The file depends on a pinned Lean/Zeta23 environment; it is not a dependency-free standalone program. Follow `BUILD_LEAN.md`. The exact arithmetic suite is small, but a cold build of the formal dependencies is substantially larger. No new full formal replay is implied by running the Python/Node suite.```
