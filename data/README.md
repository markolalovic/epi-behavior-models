# Data directory
This directory contains the input data, and processed data:

- [`./raw/`](./raw/): contains raw time series data for COVID-19 mortality in the US
    * See [`./raw/README.md`](./raw/README.md) for details.

- [`./processed/`](./processed/): contains processed time series data for COVID-19 mortality in the US
    * See [`./processed/README.md`](./processed/README.md) for details.

- [`./synthetic/`](./synthetic/): contains generated synthetic time series data
    * See [`src/run_synthetic_experiments.py`](../src/run_synthetic_experiments.py) for details.

- [`./local-db/`](./local-db/): temporary local database SQLite files for logging and storing results from ABC-SMC runs
    * See [`./local-db/README.md`](./local-db/README.md) for details on usage and configuration.

- [`./plotting/`](./plotting/): temporary files used for plotting results from ABC-SMC runs
    * See [./plotting/README.md`](./plotting/README.md) for details. 

