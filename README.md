# epi-behavior-models-dev
Parameter inference of epidemic models that incorporate behavior using ABC method based on sequential Monte Carlo (ABC-SMC).

It includes:
- Parameter inference using `pyABC` library.
- Analysis and diagnostics of the inference results.

### Models:
1. Basic SEIRD model called `model_baseline`
2. Extended SEIRD model with behavior feedback loop, using adjusted `beta` with parameter `zeta`, this is called `model_behavior`
3. SEIRPD model, which is delayed version of behavior model, with additional parameter `tau`, called  `model_behavior_delayed`

## Examples
In the examples below we generate synthetic data from empirically plausible parameters and recover them using the ABC-SMC. 

### 1. Baseline model (SEIRD)
Recovery of parameters under the assumption of constant transmission.
```bash
python3 src/run_example.py --model baseline
```
| Estimated posterior (KDE) | Model fit to synthetic data |
| :---: | :---: |
| ![Baseline Posteriors](./figures/examples/posteriors_baseline.svg) | ![Baseline Fit](./figures/examples/fit_baseline.svg) |

*The red dashed lines and markers indicate the ground-truth parameters used to generate the synthetic observations.*

### 2. Behavioral model (Mixed form)
Recovery of parameters including the behavioral sensitivity $\zeta$.
```bash
python3 src/run_example.py --model behavior
```
| Estimated posterior (KDE) | Model fit to synthetic data |
| :---: | :---: |
| ![Behavior Posteriors](./figures/examples/posteriors_behavior.svg) | ![Behavior Fit](./figures/examples/fit_behavior.svg) |

*The red dashed lines and markers indicate the ground-truth parameters used to generate the synthetic observations.*

### Installation
Clone this repository:
```bash
git clone https://github.com/markolalovic/epi-behavior-models.git
cd epi-behavior-models
```

Install the necessary Python packages in a virtual environment:
```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install numpy
python3 -m pip install scipy
python3 -m pip install pandas
python3 -m pip install matplotlib
python3 -m pip install pyabc
```

## TODO:
Later, create a requirements file to simplify installation, i.e.:
```bash
python3 -m pip install -r requirements.txt
```