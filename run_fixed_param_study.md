# Parametric Study of Turbulence Intensity and Eddy Viscosity Ratio

**File:** `run_fixed_param_study.py`

## Purpose

Conduct CFD simulations for laminar and turbulence models (Menter SST, Langtry Menter SST) with fixed turbulence intensity and eddy viscosity ratio parameters. Supports both single analysis and polar sweep tasks. Saves convergence histories and aerodynamic function evaluations for post-processing.

---

## Features

- Supports multiple turbulence models including laminar, Menter SST, and Langtry Menter SST.
- Runs simulations on specified grid files with fixed Mach, Reynolds number, and angle of attack.
- Saves convergence history in pickle format.
- Two task modes:
  - **analysis:** Single simulation run per parameter set.
  - **polar:** Sweep angle of attack from 0° to 5° and extract lift and drag.
- MPI parallel support for distributed runs.
- Handles directory creation and output organization per simulation.

---

## Input Parameters

- `--output` : Output directory root (default `"sim"`).
- `--task` : `"analysis"` or `"polar"` (default `"analysis"`).
- `--turbulenceModels` : List of turbulence models (default set to `["Langtry Menter SST"]` if not specified).

---

## Simulation Setup

- Fixed grids: `sd7003_L2.cgns`, `sd7003_L1.cgns`.
- Fixed flow conditions: Mach 0.1, Reynolds 60,000, angle of attack 4° (for analysis task).
- Fixed turbulence intensity and eddy viscosity ratio (single values).
- AeroProblem configured with reference length and flow parameters.
- ADflow solver configured with solver and turbulence options.

---

## Workflow

1. Parse command-line arguments.
2. Setup MPI environment and output folders.
3. Loop over turbulence intensities, eddy viscosity ratios, turbulence models, flow conditions, and grids.
4. Setup solver and problem configuration.
5. Run CFD solver with ADflow.
6. Save convergence history and aerodynamic functions.
7. For polar task, sweep alpha, run simulations, collect lift and drag, and print results.

---

## Output

- Convergence histories saved in `conv_hist/` as pickle files.
- Console output summarizes aerodynamic function values.
- Organized output directory structure by grid, model, and flow parameters.

---

## Dependencies

- Python packages: `numpy`, `argparse`, `os`, `pickle`, `shutil`, `mpi4py`
- CFD solver: `adflow`
- Aerodynamics base: `baseclasses`

---

## Usage Examples

Run a single analysis:

```bash
python run_fixed_param_study.py --task analysis --output sim --turbulenceModels "Langtry Menter SST"