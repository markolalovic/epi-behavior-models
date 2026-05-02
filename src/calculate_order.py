#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/calculate_order.py

Calculates location order for fit figures.

Locations are ordered by decreasing median posterior NSSE distance `Baseline - Behavioral model`

Locations with the largest improvement from behavioral feedback appear first in the figure.

Inputs:
- results/model_comparison/median_nssr_distances_{RESULT}.csv

Outputs:
- data/plotting/order_mixed.csv
- data/plotting/order_exp.csv
- data/plotting/order_rational.csv
"""

import os
import pandas as pd

from config import RESULT


def save_order(df, behavior_col, out_path):
    cols = ["Abbr", "Baseline", behavior_col]
    order = df[cols].copy()

    order["Difference"] = order["Baseline"] - order[behavior_col]

    order = order.rename(columns={"Abbr": "Location"})
    order = order[["Location", "Difference", "Baseline", behavior_col]]

    order = order.sort_values(
        by=["Difference", "Location"],
        ascending=[False, True]
    ).reset_index(drop=True)

    order["Difference"] = order["Difference"].round(4)
    order["Baseline"] = order["Baseline"].round(4)
    order[behavior_col] = order[behavior_col].round(4)

    order.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")
    print(order.to_string(index=False))
    print()


if __name__ == "__main__":
    input_path = f"../results/model_comparison/median_nssr_distances_{RESULT}.csv"
    out_dir = "../data/plotting"

    df = pd.read_csv(input_path)

    save_order(
        df,
        behavior_col="Behavioral (Mixed)",
        out_path=os.path.join(out_dir, "order_mixed.csv")
    )

    save_order(
        df,
        behavior_col="Behavioral (Exponential)",
        out_path=os.path.join(out_dir, "order_exp.csv")
    )

    save_order(
        df,
        behavior_col="Behavioral (Rational)",
        out_path=os.path.join(out_dir, "order_rational.csv")
    )