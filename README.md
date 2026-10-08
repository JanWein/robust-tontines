# Robust Fairness in Open Heterogeneous Tontines

[Read the paper (PDF)](Tontine-Robust-Fairness-Paper-EN.pdf) | [Download the complete source package](Tontine-Paper-Sources-and-R-Code-EN.zip)

Scientific working paper, Version 1.0, English edition, 2 October 2026.

This is the complete English translation of the German working paper.
The mathematical results, proofs, assumptions and numerical results are retained.
Scientific priority and readiness for publication have not been confirmed.
The manuscript is an AI-assisted research draft for expert review; authorship
and institutional affiliation have not yet been assigned.

## Files

- `manuscript.tex`: complete editable English LaTeX manuscript.
- `references.bib`: bibliography with verifiable DOI or links.
- `robust_tontine.R`: independent state checks and event-driven simulation.
- `validate_frontier.py`: independent linear optimization checks.
- `build_results.py`: English result tables and analytical figure.
- `results/`: executed numerical results, R session information and Python versions.
- `frontier.pdf`, `numerical_results.tex`: figure and generated tables.
- `Tontine-Robust-Fairness-Paper-EN.pdf`: typeset English manuscript.

## Reproduction

Requirements: R with Base R; Python with NumPy, SciPy, Pandas and Matplotlib;
LaTeX with pdfLaTeX, BibTeX, amsmath, amsthm, natbib, booktabs, lmodern and Babel.

```sh
Rscript robust_tontine.R 10000 results
python3 validate_frontier.py
python3 build_results.py
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
bibtex manuscript
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

Seed: 20261002. There are five simulation scenarios with 10,000 paths each,
a 20-year horizon and explicit payment of the terminal account balance.
The age parameterization determines an individual intensity that is then held
constant. No empirical calibration or solution of optimal private-information
exit strategies is claimed.

The English edition retains the numerical result files from the German edition.
Simulations were not rerun solely for translation. Tables and the figure were
regenerated in English from the unchanged result files. Numbered equations,
displayed mathematics, statement structure and result-file checksums were checked
against the source edition.

`individual.csv` distinguishes payout rates conditional on survival from rates
conditional on survival and active membership. Individual present values are
measured at entry. Aggregate present values are consistently measured at time zero.

The scope of the theorems is explicitly restricted to nonnegative balanced
death-contingent transfers, no external subsidies, no additional running
compensation premiums and no common death jumps. The closed-form robustness
frontier uses a common relative error half-width.

## Review before publication

Independent mathematical review, a complete priority assessment, authorship
decisions and disclosure according to the target journal's requirements remain
necessary. The mere combination of open accounts, death benefits and exits is
explicitly not presented as a new invention.
