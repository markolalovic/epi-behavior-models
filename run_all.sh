#!/bin/bash

RUN_DATA=false

# time-consuming:
RUN_INFERENCE=false
RUN_SELECTION=false
RUN_SYNTHETIC=false

RUN_POSTPROCESS=false
RUN_PLOTS=true

# to set RUN_INFERENCE=true: ./run_all.sh --inference 
# to set RUN_SELECTION=true: ./run_all.sh --selection 
for arg in "$@"; do
  case $arg in
    --data)      RUN_DATA=true ;;
    --inference) RUN_INFERENCE=true ;;
    --selection) RUN_SELECTION=true ;;
    --synthetic) RUN_SYNTHETIC=true ;;
    --post)      RUN_POSTPROCESS=true ;;
    --all)       RUN_DATA=true; RUN_INFERENCE=true; RUN_SELECTION=true; RUN_SYNTHETIC=true; RUN_POSTPROCESS=true; RUN_PLOTS=true ;;
  esac
done

if [ "$RUN_DATA" = true ]; then
    echo "Processing mortality data..."
    cd src
    python3 process_mortality_data.py
    cd ..
    echo "Done"
fi

if [ "$RUN_INFERENCE" = true ]; then
    echo "Running ABC-SMC parameter inference (30 locations)..."
    cd src
    python3 run_inference.py
    cd ..
    echo "Done"
fi

if [ "$RUN_SELECTION" = true ]; then
    echo "Running ABC-SMC model selection..."
    cd src
    python3 run_selection.py
    echo "Done"
fi

if [ "$RUN_SYNTHETIC" = true ]; then
    echo "Running synthetic experiments..."
    cd src
    python3 run_synthetic_experiments.py
    cd ..
    echo "Done"
fi

if [ "$RUN_POSTPROCESS" = true ]; then
    cd src 
    echo "Creating fits plotting data..."
    python3 create_fits_data.py
    python3 create_r0_data.py
    python3 create_final_epi_size_data.py
    echo "Done"

    echo "Extracting statistics for tables..."
    python3 median_nssr_distances.py
    echo "Done"

    echo "Extracting statistics for synthetic-data figures..."
    python3 calculate_synthetic_r0_and_final_sizes.py
    python3 create_synthetic_fits_data.py
    python3 create_synthetic_r0_data.py
    python3 create_synthetic_final_sizes_data.py
    echo "Done"

    echo "Creating effective transmission rate beta(t) plotting data..."
    python3 create_beta_eff_data.py
    echo "Done"

    echo "Creating effective reproduction numbers R_e plotting data..."
    python3 create_re_data.py
    echo "Done"
fi

if [ "$RUN_PLOTS" = true ]; then
    echo "Generating figures..."
    
    echo "Plotting S1: Excluded locations..."
    Rscript R/plot_bad_locations.R

    echo "Plotting Figures 1, S4, S5: Model fits..."
    Rscript R/plot_fits_data.R

    echo "Plotting Figures 2, S6, S7: R_0 and final sizes..."
    Rscript R/plot_r0_and_final_size_boxplots.R
    
    echo "Generating Tables 1, S1, S2..."
    Rscript R/generate_tables.R

    echo "Plotting Figures 3, S2, S3: Synthetic experiments..."
    Rscript R/plot_synthetic_results.R

    echo "Plotting Figures S9, S10, S11: Effective transmission rates..."
    Rscript R/plot_beta_eff_all.R

    echo "Plotting Figures S12, S13, S14: Effective reproduction numbers..."
    Rscript R/plot_re_eff_all.R
    echo "Done"    
fi
