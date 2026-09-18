"""Regenerate the cross-language parity fixtures under inst/extdata.

Takes a seeded 300-gene subset of GSE11923 (mouse liver, hourly for 48 h,
Hughes et al. 2009) and records what the canonical Python implementation
(par2-circadian) returns for each gene. The R test suite refits the same series
and fails if any gene diverges by more than a few multiples of machine epsilon.

Usage:
    PYTHONPATH=/path/to/par2discovery python3 data-raw/make_reference.py \
        /path/to/GSE11923_Liver_1h_48h_genes.csv
"""

import csv
import sys

import numpy as np

from par2 import __version__, fit_ar2

N_GENES = 300
SEED = 20260918
OUT_SERIES = "inst/extdata/gse11923_subset.csv"
OUT_REFERENCE = "inst/extdata/python_reference.csv"


def number(value) -> str:
    """Full-precision repr that R's read.csv can parse, with NA for missing."""
    if value is None:
        return "NA"
    value = float(value)
    if value != value:
        return "NA"
    return repr(value)


def main(source: str) -> None:
    with open(source) as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = [r for r in reader if len(r) == len(header)]

    rng = np.random.default_rng(SEED)
    idx = sorted(rng.choice(len(rows), size=N_GENES, replace=False).tolist())
    subset = [rows[i] for i in idx]

    with open(OUT_SERIES, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(subset)

    with open(OUT_REFERENCE, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "gene",
                "eigenvalue",
                "phi1",
                "phi2",
                "r2",
                "root_type",
                "half_life",
                "eigenperiod",
            ]
        )
        for row in subset:
            values = [float(v) for v in row[1:]]
            fit = fit_ar2(values)
            writer.writerow(
                [
                    row[0],
                    number(fit["eigenvalue"]),
                    number(fit["phi1"]),
                    number(fit["phi2"]),
                    number(fit["r2"]),
                    fit["root_type"],
                    number(fit["half_life"]),
                    number(fit["eigenperiod"]),
                ]
            )

    print(f"wrote {N_GENES} genes; par2-circadian {__version__}")


if __name__ == "__main__":
    main(sys.argv[1])
