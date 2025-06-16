# Turbulence and Laminar Flow Simulations with ADflow

**File:** `run_turb_laminar_sims.py`

## Purpose

Run CFD simulations for laminar and turbulence models (Menter SST, Langtry Menter SST) using the ADflow solver. Supports steady-state analysis and polar angle-of-attack sweeps. Saves convergence histories for post-processing and comparison.

---

## Features

- Supports turbulence models: laminar (mapped to Menter SST internally), Menter SST, Langtry Menter SST.
- Runs simulations on specified grid files and flow conditions.
- Saves convergence history in pickle format.
- Supports two task modes:
  - **analysis:** Run single simulation per parameter set.
  - **polar:** Sweep angle of attack from 0° to 5° and record lift and drag coefficients.
- MPI parallel support for distributed execution.
- Output directory organized by grid, turbulence model, Mach, Reynolds number, and angle of attack.

---

## Inputs

- Grid files to run (e.g., `"naca0015_L0.cgns"`).
- Tuples of flow conditions: `(Mach, Reynolds number, angle of attack)`.
- Turbulence models specified via command-line or defaulted internally.

---

## Workflow

1. Parse command-line arguments.
2. Set turbulence models to `["Langtry Menter SST"]` if default laminar specified.
3. Initialize MPI communicator and output directories.
4. Loop over turbulence models, flow conditions, and grids:
   - Setup solver options based on turbulence model.
   - Initialize ADflow solver.
   - Define AeroProblem with flow parameters.
   - Run ADflow in analysis or polar mode.
   - Save convergence history.
   - Print aerodynamic coefficients and convergence info.

---

## Key Functions

- `save_conv_history(Solver, AP, model)`: Saves convergence history as pickle file.
- `handle_remove_readonly`: Utility for removing read-only files when cleaning directories.

---

## Output

- Convergence histories saved under `conv_hist/` directory.
- Console output shows aerodynamic function values per run.
- Organized output folders for results per run.

---

## Dependencies

- Python packages: `numpy`, `argparse`, `os`, `stat`, `pickle`, `shutil`, `mpi4py`
- CFD solver: `adflow`
- Aerodynamics base: `baseclasses`

---

## Usage Examples

Run analysis:

```bash
python run_turb_laminar_sims.py --task analysis --output sim --turbulenceModels "Langtry Menter SST"