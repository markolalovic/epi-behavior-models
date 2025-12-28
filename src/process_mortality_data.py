#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""src/process_mortality_data.py
Processes raw mortality data into `data/processed/smoothed_mortality.csv`
Columns:
  date, CT, IL, MA, MI, NY, PA, ...
Dates:
  2020-03-01 ... 2020-07-01 (inclusive)

Outputs: 
- data/processed/smoothed_mortality.csv
- data/processed/population_sizes.csv
"""

import os
import numpy as np
import pandas as pd
from utils import trailing_ma7

from config import LOCATIONS, LOCATIONS_BAD, LOCATION_NAME, START_DATE, END_DATE

K = 7        # trailing MA window
PAD = K - 1  # 6 leading zeros from trailing MA(7)

if __name__ == "__main__":
    # -------------------- config --------------------
    INPUT_CSV = "../data/raw/time_series_covid19_deaths_US.csv"
    OUTPUT_DIR = "../data/processed"

    OUTPUT_SMOOTHED = os.path.join(OUTPUT_DIR, "smoothed_mortality.csv")
    OUTPUT_POPS = os.path.join(OUTPUT_DIR, "population_sizes.csv")

    df_all = pd.read_csv(INPUT_CSV)

    # date range bounds
    start_date = pd.to_datetime(START_DATE)
    end_date   = pd.to_datetime(END_DATE)
    extended_start = start_date - pd.Timedelta(days=PAD)

    # date index for output
    dates = pd.date_range(start_date, end_date, freq="D")
    df_out = pd.DataFrame({"date": dates})

    # collect population sizes
    population_size = {}

    all_locations = sorted(LOCATIONS + LOCATIONS_BAD)
    for abbr in all_locations:
        state_name = LOCATION_NAME[abbr]
        print(f"Processing {abbr} ...")

        # extract state rows
        mask = df_all["Province_State"] == state_name
        df_state_all = df_all.loc[mask]

        # sum population across state rows
        population_size[abbr] = int(df_state_all["Population"].sum())

        # get cumulative deaths time series columns 
        # JHU format: first 12 cols are metadata
        ts_cols = df_state_all.iloc[:, 12:]
        cumulative = ts_cols.sum(axis=0)

        # convert to dataframe with datetime
        df_state = cumulative.reset_index()
        df_state.columns = ["date", "cumulative_daily_deaths"]
        df_state["date"] = pd.to_datetime(df_state["date"], format="%m/%d/%y", errors="coerce")
        df_state = df_state.sort_values("date").dropna(subset=["date"]).reset_index(drop=True)

        # filter extended window for smoothing
        mask_ext = (df_state["date"] >= extended_start) & (df_state["date"] <= end_date)
        df_ext = df_state.loc[mask_ext].copy().reset_index(drop=True)

        # daily deaths as float 
        # allow for small negatives due to JHU revisions
        df_ext["daily_deaths"] = (
            df_ext["cumulative_daily_deaths"].diff().fillna(0.0).astype(float)
        )

        # apply trailing MA(7), produces leading zeros
        avg7_ext = trailing_ma7(df_ext["daily_deaths"].to_numpy())

        # trim the zeros to align to [start_date, end_date]
        avg7_trim = avg7_ext[PAD:] # length should match target length
        dates_trim = df_ext["date"].iloc[PAD:].to_numpy()

        # keep exactly the target window
        mask_target = (dates_trim >= start_date) & (dates_trim <= end_date)
        obs_smoothed = np.asarray(avg7_trim[mask_target], dtype=float)

        # clip negative values caused by revisions
        obs_smoothed = np.maximum(obs_smoothed, 0.0)

        # insert column into the wide output frame
        df_out[abbr] = obs_smoothed

    # save wide single csv
    df_out.to_csv(OUTPUT_SMOOTHED, index=False, date_format="%Y-%m-%d")
    print(f"\nSaved smoothed mortality to: {OUTPUT_SMOOTHED}")
    print(f"Columns: {list(df_out.columns)}  |  rows: {len(df_out)}")

    # save population sizes
    pd.DataFrame(
        {"state": list(population_size.keys()),
         "population": list(population_size.values())}
    ).to_csv(OUTPUT_POPS, index=False)
    print(f"Saved population sizes to: {OUTPUT_POPS}")

