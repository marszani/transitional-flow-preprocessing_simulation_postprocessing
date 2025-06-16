# Turbulence and Laminar Flow Simulation with Restart Support

**File:** `run_param_turb_laminar_restart.py`

## Purpose

Automate CFD simulations for laminar and turbulence models (Menter SST, Langtry Menter SST) using the ADflow solver. Supports restart files for continuing simulations and saves convergence histories for post-processing. Allows running steady-state analyses or polar sweeps over angle of attack.

---

## Features

- Runs multiple turbulence models including laminar, Menter SST, and Langtry Menter SST.
- Supports reading and using restart files for faster convergence.
- Saves convergence histories as pickle files.
- Supports two main tasks:
  - **analysis:** single simulation run per parameter set.
  - **polar:** sweep over angle of attack and extract lift and drag coefficients.
- Saves aerodynamic coefficient summaries in text files (`analysis_summary.txt`, `polar_summary.txt`).

---

## Input Parameters

- **Grid files:** List of grids to simulate (e.g., `sd7003_L2.cgns`, `sd7003_L1.cgns`).
- **Flow conditions:** Mach number, Reynolds number, angle of attack.
- **Turbulence intensities** and **eddy viscosity ratios** as simulation parameters.
- **Turbulence models:** Specified via command-line or defaults to Langtry Menter SST.

---

## Workflow

1. **Argument Parsing:**  
   Handles output directory, task type (`analysis` or `polar`), and turbulence models.

2. **Setup:**  
   Creates output directories and sets MPI synchronization.

3. **For each turbulence model, grid, and flow condition:**
   - Construct output directories.
   - Generate restart file path and check its existence.
   - If restart file is present:
     - Setup solver options.
     - Initialize ADflow solver.
     - Define AeroProblem with flow and geometry parameters.
     - Run ADflow in specified mode.
     - Save convergence history.
     - Save aerodynamic results summary.

4. **Polar Task:**  
   - Sweep alpha from 0° to 5° in 1° increments.
   - Run simulations, save convergence histories, and print lift and drag results.

---

## Important Functions

- `check_restart_file(filepath)`: Checks if restart file exists.
- `save_conv_history(Solver, AP, model)`: Saves convergence history pickle file.
- ADflow solver setup with detailed aerodynamic and solver options.

---

## Output

- Convergence histories saved in `conv_hist/`.
- Analysis results appended to `analysis_summary.txt`.
- Polar sweep results appended to `polar_summary.txt`.
- Console prints detailed simulation status.

---

## Dependencies

- Python packages: `numpy`, `argparse`, `os`, `shutil`, `pickle`, `mpi4py`
- CFD solver: `adflow`
- Aerodynamics base: `baseclasses`

---

## Usage Example

```bash
python run_param_turb_laminar_restart.py --task analysis --output sim --turbulenceModels "Langtry Menter SST"