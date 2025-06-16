Here is the cleaned and finalized version of the Markdown documentation for run_param_study.py, formatted for inclusion as a single .md file:

# Parametric Study of Turbulence Intensity and Eddy Viscosity Ratio

**File:** `run_param_study.py`

## Purpose

Conduct CFD simulations for laminar and turbulence models (Menter SST, Langtry Menter SST) while varying free-stream turbulence intensity (`turbIntensityInf`) and eddy viscosity ratio (`eddyVisInfRatio`). Supports both steady-state analysis and angle-of-attack sweep (polar) tasks. Stores convergence histories and aerodynamic coefficients for post-processing.

---

## Features

- Supports multiple turbulence models:
  - Laminar (internally mapped to "Menter SST")
  - Menter SST
  - Langtry Menter SST
- Runs simulations on multiple grids and flow conditions.
- Two task modes:
  - **analysis**: Single run at fixed AoA.
  - **polar**: AoA sweep from 0° to 5°.
- Saves convergence history as `.pkl`.
- Organized output per simulation case.
- Supports MPI parallel execution.

---

## Input Parameters

- `--output`: Output directory root (default: `"sim"`).
- `--task`: `"analysis"` or `"polar"` (default: `"analysis"`).
- `--turbulenceModels`: List of turbulence models (defaults to `"Langtry Menter SST"` if `"laminar"` is passed).

---

## Simulation Setup

- **Grids**:
  - `sd7003_L2.cgns`
  - `sd7003_L1.cgns`

- **Flow Conditions**:
  - Example tuple: Mach 0.1, Reynolds 60,000, AoA 4°
  - Extendable list of `(mach, reynolds, alpha)` tuples

- **Parameter Ranges**:
  - Turbulence Intensity: `[0.0008]`
  - Eddy Viscosity Ratio: `[0.001]`

- **ADflow Solver Settings**:
  - RANS or laminar NS equations
  - Steady-state solver
  - Turbulence modeling and discretization settings
  - Output: surface/volume data, residuals, function evaluations

---

## Workflow

1. Parse command-line arguments.
2. Set up MPI communication and output directories.
3. Loop over:
   - Turbulence intensities
   - Eddy viscosity ratios
   - Turbulence models
   - Flow conditions
   - Grid files
4. For each simulation:
   - Configure solver options.
   - Initialize `ADFLOW` and define an `AeroProblem`.
   - Run simulation.
   - Evaluate functions and save convergence history.
5. For `polar` task:
   - Sweep AoA from 0° to 5° in 1° steps.
   - Store and print lift and drag coefficients.

---

## Output

- Convergence histories saved in:

conv_hist/conv_hist_.pkl

- Output directories structured as:

///M_Re_alpha_TI_EVR/

- Terminal output includes:
- Residual convergence
- Function values: `cl`, `cd`, `cmz`
- AoA sweep summaries for polar runs

---

## Dependencies

- `numpy`
- `argparse`
- `os`
- `pickle`
- `shutil`
- `stat`
- `mpi4py`
- `adflow`
- `baseclasses`

---

## Example Usage

Run a single analysis case:
```bash
python run_param_study.py --task analysis --output sim --turbulenceModels "Langtry Menter SST"

---

Notes
	•	If "laminar" is passed, it is internally remapped to "Langtry Menter SST".
	•	Restart file functionality is included in comments and can be activated.
	•	Additional models, grid levels, and parameter combinations can be added by extending the lists at the top of the script.
