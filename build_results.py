"""Build manuscript tables and the analytic frontier figure from saved results."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ind = pd.read_csv("results/individual.csv")
agg = pd.read_csv("results/aggregate.csv")
states = pd.read_csv("results/state_checks.csv")
lp = pd.read_csv("results/lp_checks.csv")

def number(x, places=4):
    return f"{x:.{places}f}"

lines = [r"\begin{table}[htbp]\centering\small",
         r"\caption{Selected individual present-value ratios. All benefits, including estate, exit and terminal-account payments, are compared with capital at entry. Intervals quantify Monte Carlo uncertainty.}",
         r"\begin{tabular}{llrrr}\toprule",
         r"Scenario & Member & PV / capital & Standard error & 95\% interval\\\midrule"]
selected = [
    ("closed_nominal", [1, 2, 3], "Closed, nominal"),
    ("closed_robust", [1, 2, 3], "Closed, robust"),
    ("open_exit_nominal", [1, 2, 3], "Open with exits"),
    ("open_exit_robust_error", [1, 2, 3], "Open, robust, error"),
    ("singleton", [1], "Singleton")]
for scenario, ids, label in selected:
    for i in ids:
        x = ind[(ind.scenario == scenario) & (ind.id == i)].iloc[0]
        lines.append(f"{label} & {i} & {number(x.epv_ratio)} & {number(x.se)} & "
                     f"$[{number(x.lower)};{number(x.upper)}]$" + r"\\")
lines += [r"\bottomrule\end{tabular}\end{table}"]
lines += [r"\begin{table}[htbp]\centering\small",
          r"\caption{Payout rate after five years relative to its value at entry, conditional only on survival. These closed scenarios have no exits.}",
          r"\begin{tabular}{llrrr}\toprule",
          r"Scenario & Member & Expected rate & Standard error & Survivors\\\midrule"]
for scenario, label in [("closed_nominal", "Nominal"), ("closed_robust", "Robust"), ("singleton", "Singleton")]:
    for _, x in ind[ind.scenario == scenario].iterrows():
        lines.append(f"{label} & {int(x.id)} & {number(x.alive_payment_ratio)} & "
                     f"{number(x.alive_payment_se)} & {int(x.alive_n)}" + r"\\")
lines += [r"\bottomrule\end{tabular}\end{table}"]

err = ind[ind.scenario == "open_exit_robust_error"]
lines += [r"\begin{table}[htbp]\centering\small",
          r"\caption{Misspecification and integrated tolerance for the first entry cohort in the open robust scenario. The relative error half-width is 20\%, and the relative annual error tolerance is 10\%.}",
          r"\begin{tabular}{rrrr}\toprule",
          r"Member & PV / capital minus 1 & Standard error & Expected bound / capital\\\midrule"]
for _, x in err[err.id <= 3].iterrows():
    lines.append(f"{int(x.id)} & {number(x.epv_ratio-1)} & {number(x.se)} & "
                 f"{number(x.expected_bound)}" + r"\\")
lines += [r"\bottomrule\end{tabular}\end{table}"]
maxstate = states[["row_error", "sym_error", "frontier_error", "worst_error"]].max().max()
maxlp = lp.relative_error.max()
maxbal = agg.max_relative_balance.max()
lines += [f"The maximum relative error in the R state checks was ${maxstate:.3g}$, "
          f"that in the 180 independent linear optimization checks was ${maxlp:.3g}$. "
          f"The maximum relative dynamic balance error was ${maxbal:.3g}$. "
          "The immediate change in other accounts at entry and exit was exactly zero. "
          "The aggregate present-value ratio, consistently discounted to time zero, was "
          "one on every path in all scenarios, up to rounding error; its Monte Carlo standard error "
          "is therefore merely a numerical rounding measure. Individual members, in contrast, have "
          "stochastic individual present values."]
# Replace computer exponent syntax in TeX prose by valid scientific notation.
import re
text = "\n".join(lines)
text = re.sub(r"\$([0-9.]+)e(-?\d+)\$", lambda m: "$"+m[1]+r"\cdot10^{"+str(int(m[2]))+"}$", text)
Path("numerical_results.tex").write_text(text+"\n", encoding="utf-8")

plt.rcParams.update({"font.family": "DejaVu Serif", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})
x = np.linspace(0, 1.2, 241); rho = np.minimum(1, x)
fig, ax = plt.subplots(figsize=(6.8, 3.5))
for s, color in [(0.4, "#163b57"), (0.65, "#387b73"), (0.9, "#b26b38")]:
    ax.plot(x, 1-rho+rho*max(0, 2*s-1), color=color, lw=2,
            label=rf"$\max_i k_i/K={s:g}$")
ax.set(xlabel=r"Relative tolerance / twice the error half-width: $\eta/(2\varepsilon)$",
       ylabel=r"Minimum death-benefit rate $\mathcal{H}_{\min}/K$",
       xlim=(0, 1.2), ylim=(-0.02, 1.04))
ax.grid(alpha=0.15); ax.legend(frameon=False, loc="lower left")
fig.tight_layout(); fig.savefig("frontier.pdf"); plt.close(fig)
print("Generated numerical_results.tex and frontier.pdf")
