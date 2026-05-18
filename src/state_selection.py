#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/state_selection.py

Location selection for the first-wave mortality analysis.

Locations are excluded using 2 criteria applied to the 7-day moving average of daily reported deaths:

- Rule 1: insufficient mortality signal
  * total smoothed mortality is below 50 deaths
  * peak smoothed daily mortality does not exceed 10 deaths per day

- Rule 2: incomplete first mortality wave
  * the peak occurs fewer than 21 days before the end of the window
  * the final 14-day mean exceeds 40% of the peak

- Colorado is additionally excluded 
  * because of anomaly / discontinuity in the cumulative mortality data / reporting

"""

import numpy as np
import pandas as pd

from config import LOCATIONS_ALL, LOCATION_NAME

MIN_TOTAL_DEATHS = 50.0
MIN_PEAK_DEATHS = 10.0

MIN_DAYS_AFTER_PEAK = 21
TAIL_DAYS = 14
MAX_TAIL_TO_PEAK = 0.40

ADDITIONAL_EXCLUDED = ["CO"]

def format_list(name, values):
    print(f"{name} = [")
    for i in range(0, len(values), 5):
        chunk = values[i:i + 5]
        quoted = [f"'{x}'" for x in chunk]
        suffix = "," if i + 5 < len(values) else ""
        print("    " + ", ".join(quoted) + suffix)
    print("]")

if __name__ == "__main__":
    data_path = "../data/processed/smoothed_mortality.csv"
    df = pd.read_csv(data_path, parse_dates=["date"])

    dates = df["date"]

    rows = []

    for loc in LOCATIONS_ALL:
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

        # Rule 1
        fail_low_deaths = (
            total_deaths < MIN_TOTAL_DEATHS or peak <= MIN_PEAK_DEATHS
        )

        # Rule 2
        fail_no_completed_wave = (
            days_after_peak < MIN_DAYS_AFTER_PEAK
            or tail_to_peak > MAX_TAIL_TO_PEAK
        )

        # data reporting anomaly exclusion
        fail_reporting_discontinuity = loc in ADDITIONAL_EXCLUDED

        excluded = (
            fail_low_deaths
            or fail_no_completed_wave
            or fail_reporting_discontinuity
        )

        if fail_low_deaths:
            reason = "insufficient mortality signal"
        elif fail_no_completed_wave:
            reason = "incomplete mortality wave"
        elif fail_reporting_discontinuity:
            reason = "reporting discontinuity"
        else:
            reason = "included"

        rows.append({
            "location": loc,
            "name": LOCATION_NAME[loc],
            "status": "excluded" if excluded else "selected",
            "reason": reason,
            "total_deaths": total_deaths,
            "peak": peak,
            "peak_date": str(peak_date),
            "days_after_peak": days_after_peak,
            "tail_mean": tail_mean,
            "tail_to_peak": tail_to_peak,
            "decline_from_peak": decline_from_peak,
            "fail_low_deaths": fail_low_deaths,
            "fail_no_completed_wave": fail_no_completed_wave,
            "fail_reporting_discontinuity": fail_reporting_discontinuity,
        })

    summary = pd.DataFrame(rows)

    selected = summary.loc[
        summary["status"] == "selected",
        "location"
    ].tolist()

    excluded = summary.loc[
        summary["status"] == "excluded",
        "location"
    ].tolist()

    display_cols = [
        "location",
        "name",
        "status",
        "reason",
        "total_deaths",
        "peak",
        "peak_date",
        "days_after_peak",
        "tail_mean",
        "tail_to_peak",
        "decline_from_peak",
    ]

    display = summary[display_cols].copy()
    display["total_deaths"] = display["total_deaths"].round(1)
    display["peak"] = display["peak"].round(2)
    display["tail_mean"] = display["tail_mean"].round(2)
    display["tail_to_peak"] = display["tail_to_peak"].round(2)
    display["decline_from_peak"] = display["decline_from_peak"].round(2)

    print("\nLocation-selection summary: ")
    print(display.to_string(index=False))

    print("\nCounts: ")
    print(f"n_total    = {len(LOCATIONS_ALL)}")
    print(f"n_selected = {len(selected)}")
    print(f"n_excluded = {len(excluded)}")


    print(
        summary.loc[
            summary["status"] == "excluded",
            "reason"
        ].value_counts().to_string()
    )

    format_list("LOCATIONS", selected)
    print()
    format_list("LOCATIONS_BAD", excluded)

    out_path = "../results/state_selection_summary.csv"
    summary.to_csv(out_path, index=False)
    print(f"\nSaved location-selection summary to: {out_path}")



'''
Location-selection summary: 
location                 name   status                        reason  total_deaths    peak  peak_date  days_after_peak  tail_mean  tail_to_peak  decline_from_peak
      AK               Alaska excluded insufficient mortality signal          13.7    0.57 2020-04-05               87       0.13          0.23               0.77
      AL              Alabama excluded     incomplete mortality wave         927.4   17.14 2020-05-12               50      10.95          0.64               0.36
      AR             Arkansas excluded insufficient mortality signal         260.6    7.00 2020-06-23                8       5.65          0.81               0.19
      AZ              Arizona excluded     incomplete mortality wave        1596.1   36.86 2020-07-01                0      28.80          0.78               0.22
      CA           California excluded     incomplete mortality wave        5879.6   80.57 2020-04-24               68      61.92          0.77               0.23
      CO             Colorado excluded       reporting discontinuity        1717.6   85.14 2020-04-07               85       3.03          0.04               0.96
      CT          Connecticut selected                      included        4314.0  113.86 2020-04-26               66       8.91          0.08               0.92
      DC District of Columbia selected                      included         548.9   12.14 2020-04-30               62       2.55          0.21               0.79
      DE             Delaware selected                      included         563.4   12.57 2020-04-28               64       1.60          0.13               0.87
      FL              Florida excluded     incomplete mortality wave        3429.1   50.71 2020-05-08               54      35.45          0.70               0.30
      GA              Georgia excluded     incomplete mortality wave        2783.6   44.29 2020-04-20               72      22.41          0.51               0.49
      HI               Hawaii excluded insufficient mortality signal          17.9    0.86 2020-04-27               65       0.05          0.06               0.94
      IA                 Iowa selected                      included         707.7   15.71 2020-05-26               36       3.78          0.24               0.76
      ID                Idaho excluded insufficient mortality signal          91.0    3.43 2020-04-15               77       0.29          0.08               0.92
      IL             Illinois selected                      included        6884.9  116.86 2020-05-13               49      40.23          0.34               0.66
      IN              Indiana selected                      included        2549.7   42.29 2020-04-27               65      12.86          0.30               0.70
      KS               Kansas excluded insufficient mortality signal         269.6    6.57 2020-04-16               76       1.76          0.27               0.73
      KY             Kentucky excluded insufficient mortality signal         558.3   10.00 2020-04-22               70       4.00          0.40               0.60
      LA            Louisiana selected                      included        3200.1   65.86 2020-04-18               74      13.04          0.20               0.80
      MA        Massachusetts selected                      included        8043.9  189.43 2020-04-25               67      30.93          0.16               0.84
      MD             Maryland selected                      included        3166.6   68.57 2020-05-08               54      16.34          0.24               0.76
      ME                Maine excluded insufficient mortality signal         104.1    2.57 2020-04-24               68       0.26          0.10               0.90
      MI             Michigan selected                      included        6161.3  145.86 2020-04-23               69      10.58          0.07               0.93
      MN            Minnesota selected                      included        1461.0   25.00 2020-05-30               32       9.83          0.39               0.61
      MO             Missouri excluded     incomplete mortality wave        1025.0   18.86 2020-04-26               66       9.11          0.48               0.52
      MS          Mississippi excluded     incomplete mortality wave        1046.6   19.29 2020-05-07               55      10.70          0.56               0.44
      MT              Montana excluded insufficient mortality signal          21.9    1.00 2020-04-22               70       0.22          0.22               0.78
      NC       North Carolina excluded     incomplete mortality wave        1360.6   22.86 2020-06-02               29      15.80          0.69               0.31
      ND         North Dakota excluded insufficient mortality signal          78.7    1.86 2020-05-12               50       0.34          0.18               0.82
      NE             Nebraska excluded insufficient mortality signal         268.4    6.14 2020-06-17               14       3.44          0.56               0.44
      NH        New Hampshire excluded insufficient mortality signal         366.7    8.29 2020-05-13               49       3.37          0.41               0.59
      NJ           New Jersey selected                      included       14973.4  345.00 2020-04-20               72      33.48          0.10               0.90
      NM           New Mexico selected                      included         492.4   10.29 2020-05-15               47       4.10          0.40               0.60
      NV               Nevada excluded insufficient mortality signal         502.1    9.71 2020-04-13               79       2.70          0.28               0.72
      NY             New York selected                      included       31900.7 1013.57 2020-04-12               80      33.55          0.03               0.97
      OH                 Ohio selected                      included        2575.7   44.71 2020-04-28               64      12.59          0.28               0.72
      OK             Oklahoma excluded insufficient mortality signal         383.1    8.29 2020-04-25               67       1.63          0.20               0.80
      OR               Oregon excluded insufficient mortality signal         203.1    3.71 2020-04-13               79       1.88          0.51               0.49
      PA         Pennsylvania selected                      included        6623.6  158.71 2020-05-05               57      28.07          0.18               0.82
      RI         Rhode Island selected                      included         968.9   20.43 2020-05-08               54       4.07          0.20               0.80
      SC       South Carolina excluded     incomplete mortality wave         719.9   15.14 2020-05-04               58       8.50          0.56               0.44
      SD         South Dakota excluded insufficient mortality signal          90.3    2.29 2020-05-06               56       1.07          0.47               0.53
      TN            Tennessee excluded     incomplete mortality wave         588.6   11.29 2020-04-08               84       8.04          0.71               0.29
      TX                Texas excluded     incomplete mortality wave        2778.6   49.14 2020-07-01                0      37.18          0.76               0.24
      UT                 Utah excluded insufficient mortality signal         168.1    3.00 2020-06-17               14       1.96          0.65               0.35
      VA             Virginia selected                      included        1731.4   34.00 2020-05-28               34      13.00          0.38               0.62
      VT              Vermont excluded insufficient mortality signal          56.0    1.86 2020-04-18               74       0.07          0.04               0.96
      WA           Washington selected                      included        1291.3   35.57 2020-04-07               85       5.80          0.16               0.84
      WI            Wisconsin selected                      included         809.4   13.00 2020-05-29               33       4.87          0.37               0.63
      WV        West Virginia excluded insufficient mortality signal          92.7    2.71 2020-04-23               69       0.36          0.13               0.87
      WY              Wyoming excluded insufficient mortality signal          20.0    0.71 2020-04-22               70       0.14          0.20               0.80

Counts: 
n_total    = 51
n_selected = 20
n_excluded = 31

Reason:
insufficient mortality signal    19
incomplete mortality wave        11
reporting discontinuity           1

LOCATIONS = [
    'CT', 'DC', 'DE', 'IA', 'IL',
    'IN', 'LA', 'MA', 'MD', 'MI',
    'MN', 'NJ', 'NM', 'NY', 'OH',
    'PA', 'RI', 'VA', 'WA', 'WI'
]

LOCATIONS_BAD = [
    'AK', 'AL', 'AR', 'AZ', 'CA',
    'CO', 'FL', 'GA', 'HI', 'ID',
    'KS', 'KY', 'ME', 'MO', 'MS',
    'MT', 'NC', 'ND', 'NE', 'NH',
    'NV', 'OK', 'OR', 'SC', 'SD',
    'TN', 'TX', 'UT', 'VT', 'WV',
    'WY'
]

Saved location-selection summary to: ../results/state_selection_summary.csv

'''