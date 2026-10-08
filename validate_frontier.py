"""Independent LP verification of the sharp frontier, not a proof."""
from pathlib import Path
import csv
import numpy as np
from scipy.optimize import linprog

rng = np.random.default_rng(20261002)
rows = []
for n in [1, 2, 3, 5, 10, 20]:
    for rep in range(30):
        k = np.exp(rng.normal(0, 3, n))
        delta = np.exp(rng.normal(0, 3, n))
        b = np.minimum(k, delta / 0.4)
        B, K = b.sum(), k.sum()
        if n == 1:
            optimum = 0.0
        else:
            edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
            incidence = np.zeros((n, len(edges)))
            for col, (i, j) in enumerate(edges):
                incidence[i, col] = incidence[j, col] = 1
            # Scale exposures to avoid conditioning problems in raw units.
            solution = linprog(-2 * np.ones(len(edges)), A_ub=incidence,
                               b_ub=b/B, bounds=(0, None), method="highs",
                               options={"primal_feasibility_tolerance": 1e-9,
                                        "dual_feasibility_tolerance": 1e-9})
            if not solution.success:
                raise RuntimeError(solution.message)
            optimum = -solution.fun * B
        predicted = K-B+max(0, 2*b.max()-B)
        error = abs(K-optimum-predicted)/max(1, K)
        assert error < 1e-8, (n, rep, error)
        rows.append(dict(n=n, rep=rep+1, lp_estate_rate=K-optimum,
                         formula_estate_rate=predicted, relative_error=error))
Path("results").mkdir(exist_ok=True)
with open("results/lp_checks.csv", "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
print("LP states:", len(rows), "maximum relative error:", max(x["relative_error"] for x in rows))
