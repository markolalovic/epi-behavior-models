#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/state_selection.py

Objective state selection for the first-wave mortality analysis.

A location is excluded based on 3 rules:

Rule 1: sufficient mortality signal:
   - total smoothed mortality in the analysis window is too small
   - peak 7-day averaged daily mortality is too small

Rule 2: completed mortality wave:
   - the peak occurs too close to the end of the observation window
   - the final tail remains too high relative to the peak

Rule 3: coherent mortality wave:
   - the high-mortality part of the trajectory is not a single dominant block
   - the smoothed trajectory contains too many large local deviations


"""

import numpy as np
import pandas as pd

from config import LOCATIONS, LOCATIONS_BAD, LOCATION_NAME


# -------------------- selection thresholds --------------------

# Rule 1: sufficient mortality signal
MIN_TOTAL_DEATHS = 50.0
MIN_PEAK_DEATHS = 10.0

# Rule 2: completed first mortality wave
MIN_DAYS_AFTER_PEAK = 21
TAIL_DAYS = 14
MAX_TAIL_TO_PEAK = 0.40

# Rule 3: single dominant high-mortality block
HIGH_MORTALITY_FRACTION = 0.50
MIN_HIGH_BLOCK_LENGTH = 5
MAX_GAP_TO_FILL = 3
MAX_HIGH_BLOCKS = 1

# Rule 3: no large local reporting irregularities
OUTLIER_WINDOW = 7
OUTLIER_ABS_THRESHOLD = 5.0
OUTLIER_REL_THRESHOLD = 0.35
MAX_OUTLIER_DAYS = 3


def fill_short_gaps(mask, max_gap):
    mask = np.asarray(mask, dtype=bool).copy()
    n = len(mask)

    i = 0
    while i < n:
        if mask[i]:
            i += 1
            continue

        start = i
        while i < n and not mask[i]:
            i += 1
        end = i

        left_true = start > 0 and mask[start - 1]
        right_true = end < n and mask[end]

        if left_true and right_true and (end - start) <= max_gap:
            mask[start:end] = True

    return mask


def count_true_blocks(mask, min_length):
    """Counts blocks with length at least min_length"""
    mask = np.asarray(mask, dtype=bool)

    blocks = []
    i = 0
    n = len(mask)

    while i < n:
        if not mask[i]:
            i += 1
            continue

        start = i
        while i < n and mask[i]:
            i += 1
        end = i

        length = end - start
        if length >= min_length:
            blocks.append((start, end, length))

    return blocks


def count_local_outliers(y, peak):
    """Counts large deviations from a centered local median"""
    local_median = pd.Series(y).rolling(
        OUTLIER_WINDOW,
        center=True,
        min_periods=1
    ).median().to_numpy()

    residual = np.abs(y - local_median)
    threshold = max(OUTLIER_ABS_THRESHOLD, OUTLIER_REL_THRESHOLD * peak)

    n_outlier_days = int(np.sum(residual > threshold))
    max_local_residual = float(np.max(residual))

    return n_outlier_days, max_local_residual


def format_list(name, values):
    print(f"{name} = [")
    for i in range(0, len(values), 6):
        chunk = values[i:i + 6]
        quoted = [f"'{x}'" for x in chunk]
        print("    " + ", ".join(quoted) + ",")
    print("]")


if __name__ == "__main__":
    data_path = "../data/processed/smoothed_mortality.csv"
    df = pd.read_csv(data_path, parse_dates=["date"])

    dates = df["date"]
    all_locations = sorted(LOCATION_NAME.keys())

    rows = []

    for loc in all_locations:
        y = df[loc].to_numpy(dtype=float)

        total_deaths = float(np.sum(y))
        peak = float(np.max(y))
        peak_idx = int(np.argmax(y))
        peak_date = dates.iloc[peak_idx].date()

        days_after_peak = len(y) - peak_idx - 1
        tail_mean = float(np.mean(y[-TAIL_DAYS:]))

        if peak > 0.0:
            tail_to_peak = tail_mean / peak
            decline_from_peak = 1.0 - tail_to_peak
        else:
            tail_to_peak = np.inf
            decline_from_peak = -np.inf

        high_mask = y >= HIGH_MORTALITY_FRACTION * peak
        high_mask = fill_short_gaps(high_mask, MAX_GAP_TO_FILL)
        high_blocks = count_true_blocks(high_mask, MIN_HIGH_BLOCK_LENGTH)

        n_outlier_days, max_local_residual = count_local_outliers(y, peak)

        # Rule 1
        fail_low_deaths = (
            total_deaths < MIN_TOTAL_DEATHS
            or peak < MIN_PEAK_DEATHS
        )

        # Rule 2
        fail_no_completed_wave = (
            days_after_peak < MIN_DAYS_AFTER_PEAK
            or tail_to_peak > MAX_TAIL_TO_PEAK
        )

        # Rule 3
        fail_not_single_wave = len(high_blocks) > MAX_HIGH_BLOCKS
        fail_noisy_trajectory = n_outlier_days > MAX_OUTLIER_DAYS

        fail_noncoherent_wave = (
            fail_not_single_wave
            or fail_noisy_trajectory
        )

        excluded = (
            fail_low_deaths
            or fail_no_completed_wave
            or fail_noncoherent_wave
        )

        if loc in LOCATIONS:
            original_status = "selected"
        elif loc in LOCATIONS_BAD:
            original_status = "excluded"
        else:
            original_status = "not listed"


        if fail_low_deaths:
            reason = "insufficient mortality signal"
        elif fail_no_completed_wave:
            reason = "no completed mortality wave"
        elif fail_noncoherent_wave:
            reason = "non-coherent mortality wave"
        else:
            reason = "included"

        rows.append({
            "location": loc,
            "name": LOCATION_NAME[loc],
            "original_status": original_status,
            "objective_status": "excluded" if excluded else "selected",
            "reason": reason,
            "total_deaths": total_deaths,
            "peak": peak,
            "peak_date": str(peak_date),
            "days_after_peak": days_after_peak,
            "tail_mean": tail_mean,
            "tail_to_peak": tail_to_peak,
            "decline_from_peak": decline_from_peak,
            "high_blocks": len(high_blocks),
            "n_outlier_days": n_outlier_days,
            "max_local_residual": max_local_residual,
            "fail_low_deaths": fail_low_deaths,
            "fail_no_completed_wave": fail_no_completed_wave,
            "fail_not_single_wave": fail_not_single_wave,
            "fail_noisy_trajectory": fail_noisy_trajectory,
            "fail_noncoherent_wave": fail_noncoherent_wave,
        })

    summary = pd.DataFrame(rows)

    selected = summary.loc[
        summary["objective_status"] == "selected",
        "location"
    ].tolist()

    excluded = summary.loc[
        summary["objective_status"] == "excluded",
        "location"
    ].tolist()

    columns = [
        "location",
        "name",
        "original_status",
        "objective_status",
        "reason",
        "total_deaths",
        "peak",
        "peak_date",
        "days_after_peak",
        "tail_mean",
        "tail_to_peak",
        "decline_from_peak",
        "high_blocks",
        "n_outlier_days",
        "max_local_residual",
    ]

    display = summary[columns].copy()
    display["total_deaths"] = display["total_deaths"].round(1)
    display["peak"] = display["peak"].round(2)
    display["tail_mean"] = display["tail_mean"].round(2)
    display["tail_to_peak"] = display["tail_to_peak"].round(2)
    display["decline_from_peak"] = display["decline_from_peak"].round(2)
    display["max_local_residual"] = display["max_local_residual"].round(2)

    print("\nObjective state-selection summary")
    print("---------------------------------\n")
    print(display.to_string(index=False))

    print("\nSelected locations")
    print("------------------")
    print(selected)
    print(f"n_selected = {len(selected)}")

    print("\nExcluded locations")
    print("------------------")
    print(excluded)
    print(f"n_excluded = {len(excluded)}")

    print("\nPreviously selected, now excluded")
    print("---------------------------------")
    print(summary.loc[
        (summary["original_status"] == "selected")
        & (summary["objective_status"] == "excluded"),
        [
            "location",
            "name",
            "reason",
            "total_deaths",
            "peak",
            "peak_date",
            "tail_to_peak",
            "high_blocks",
            "n_outlier_days",
            "max_local_residual",
        ]
    ].to_string(index=False))

    print("\nPreviously excluded, now selected")
    print("---------------------------------")
    print(summary.loc[
        (summary["original_status"] == "excluded")
        & (summary["objective_status"] == "selected"),
        [
            "location",
            "name",
            "reason",
            "total_deaths",
            "peak",
            "peak_date",
            "tail_to_peak",
            "high_blocks",
            "n_outlier_days",
            "max_local_residual",
        ]
    ].to_string(index=False))

    print("\nExclusion counts by reason")
    print("--------------------------")
    print(
        summary.loc[
            summary["objective_status"] == "excluded",
            "reason"
        ].value_counts().to_string()
    )

    print("\nConfig-style lists")
    print("------------------\n")
    format_list("LOCATIONS_OBJECTIVE", selected)
    print()
    format_list("LOCATIONS_BAD_OBJECTIVE", excluded)

    out_path = "../results/state_selection_summary.csv"
    summary.to_csv(out_path, index=False)
    print(f"\nSaved state-selection summary to: {out_path}")

'''
Selected locations
------------------
['CO', 'CT', 'DC', 'DE', 'IA', 'IL', 'IN', 'LA', 'MA', 'MD', 'MI', 'MN', 'NJ', 'NM', 'NY', 'OH', 'PA', 'RI', 'VA', 'WA']
n_selected = 20

Excluded locations
------------------
['AK', 'AL', 'AR', 'AZ', 'CA', 'FL', 'GA', 'HI', 'ID', 'KS', 'KY', 'ME', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NV', 'OK', 'OR', 'SC', 'SD', 'TN', 'TX', 'UT', 'VT', 'WI', 'WV', 'WY']
n_excluded = 31

Previously selected, now excluded
---------------------------------
location          name                        reason  total_deaths     peak  peak_date  tail_to_peak  high_blocks  n_outlier_days  max_local_residual
      ID         Idaho insufficient mortality signal     91.000000 3.428571 2020-04-15      0.083333            1               0            0.857143
      ND  North Dakota insufficient mortality signal     78.714286 1.857143 2020-05-12      0.181319            1               0            0.571429
      NH New Hampshire insufficient mortality signal    366.714286 8.285714 2020-05-13      0.406404            1               0            1.571429
      NV        Nevada insufficient mortality signal    502.142857 9.714286 2020-04-13      0.278361            1               0            2.000000
      OK      Oklahoma insufficient mortality signal    383.142857 8.285714 2020-04-25      0.197044            1               0            1.285714

Previously excluded, now selected
---------------------------------
Empty DataFrame
Columns: [location, name, reason, total_deaths, peak, peak_date, tail_to_peak, high_blocks, n_outlier_days, max_local_residual]
Index: []

Exclusion counts by reason
--------------------------
reason
insufficient mortality signal    18
no completed mortality wave      11
non-coherent mortality wave       2

Config-style lists
------------------

LOCATIONS_OBJECTIVE = [
    'CO', 'CT', 'DC', 'DE', 'IA', 'IL',
    'IN', 'LA', 'MA', 'MD', 'MI', 'MN',
    'NJ', 'NM', 'NY', 'OH', 'PA', 'RI',
    'VA', 'WA',
]

LOCATIONS_BAD_OBJECTIVE = [
    'AK', 'AL', 'AR', 'AZ', 'CA', 'FL',
    'GA', 'HI', 'ID', 'KS', 'KY', 'ME',
    'MO', 'MS', 'MT', 'NC', 'ND', 'NE',
    'NH', 'NV', 'OK', 'OR', 'SC', 'SD',
    'TN', 'TX', 'UT', 'VT', 'WI', 'WV',
    'WY',
]

Saved state-selection summary to: ../results/state_selection_summary.csv
'''