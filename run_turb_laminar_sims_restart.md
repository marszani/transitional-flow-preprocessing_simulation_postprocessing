# Laminar and Turbulence Model Simulations with Restart Support

**File:** `run_turb_laminar_sims_restart.py`

## Purpose

Run CFD simulations for laminar and turbulence models (Menter SST, Langtry Menter SST) using the ADflow solver. Supports restart files for continuing previous simulations and saves convergence histories. Supports steady-state analysis and polar angle-of-attack sweeps.

---

## Features

- Supports turbulence models: laminar (runs as Menter SST), Menter SST, Langtry Menter SST.
- Uses restart files if available to continue simulations.
- Saves convergence history to pickle files for post-processing.
- Supports two tasks:
  - **analysis:** run single case.
  - **polar:** sweep angle of attack from 0° to 5°.
- Uses MPI for parallel execution.
- Organizes output in structured directories by grid, turbulence model, Mach, Reynolds, and angle of attack.

---

## Inputs

- Grid files: e.g., `e387_L0.cgns`.
- Flow conditions: Mach, Reynolds number, angle of attack tuples.
- Turbulence models specified by command-line argument.

---

## Workflow

1. Parse command-line arguments.
2. Setup MPI and output directories.
3. For each turbulence model, Mach-Reynolds-alpha tuple, and grid:
   - Create output directory.
   - Construct restart file path.
   - Check restart file existence.
   - If restart file exists:
     - Setup solver options and initialize ADflow solver.
     - Define AeroProblem with flow parameters.
     - Run simulation with ADflow.
     - Save convergence history.
     - Print summary of aerodynamic functions.
4. For polar task:
   - Sweep alpha in increments.
   - Run simulation for each alpha.
   - Save and print aerodynamic coefficients (CL, CD).

---

## Key Functions

- `save_conv_history(Solver, AP, model)`: Saves convergence history pickle file.
- `check_restart_file(filepath)`: Checks if restart file exists.

---

## Output

- Convergence histories saved in `conv_hist/` folder.
- Summary printouts for each case.

---

## Dependencies

- Python packages: `numpy`, `argparse`, `os`, `pickle`, `shutil`, `mpi4py`.
- CFD solver: `adflow`.
- Aerodynamics base: `baseclasses`.

---

## Usage Example

```bash
python run_turb_laminar_sims_restart.py --task analysis --output sim --turbulenceModels "Langtry Menter SST"