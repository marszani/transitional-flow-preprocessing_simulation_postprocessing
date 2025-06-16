# Convergence History Plotting

**File:** `plot_convergence.py`

## Purpose

Load, normalize, and visualize CFD convergence histories from saved pickle (`.pkl`) files. Generate plots of normalized residuals, turbulence quantities, and aerodynamic coefficients against iterations to assist in solution monitoring and analysis.

---

## Input Data

- Loads all `.pkl` files found in the `conv_hist` subdirectory relative to the script location.
- Each `.pkl` file contains a dictionary with keys such as:
  - Residual norms (e.g., `RSDMassRMS`, `RSDTurbulentEnergyKineticRMS`, `RSDTurbulentDissRateRMS`, `RSDTransitionGammaRMS`, `RSDTransitionReThetat`)
  - Aerodynamic coefficients (`CoefLift`, `CoefDrag`, `CoefMomentZ`)
  - Iteration counts (`total minor iters`)
  - Additional metadata in `map` dictionary (optional)

---

## Processing Steps

1. **Data Loading:**  
   Reads each `.pkl` file and inspects its contents, printing key summary info.

2. **Normalization:**  
   For residuals and turbulence quantities, normalize each variable by the mean of its first 5 values to enable comparison across runs.

3. **Directory Setup:**  
   Creates a subdirectory in `convergence_history/<filename_without_extension>/` to store plots.

4. **Plot Generation:**

   - **Normalized Residuals Plot:**  
     Log-scale plot of normalized residuals and turbulence RMS values vs iterations.

   - **Normalized Variables Plot:**  
     Linear-scale plot of normalized residuals vs iterations.

   - **Aerodynamic Coefficients Plots:**  
     Individual plots of lift (`Cl`), drag (`Cd`), and moment (`CMz`) coefficients vs iterations, annotated with final values.

   - **Additional Variables:**  
     Plots for CFL number and iteration types if available.

   - **Map Variables:**  
     Plots for any variables under the `map` dictionary key vs iterations.

---

## Plot Styling

- Uses `matplotlib` with LaTeX rendering for axis labels and legends.
- Residuals plots use logarithmic y-axis with limits set between `1e-15` and `1`.
- Grid lines, legends, and titles are configured for clarity and print quality.
- Plots saved as high-resolution PNG files (`dpi=600`).

---

## Usage

- Place `.pkl` convergence history files in the `conv_hist` directory.
- Run the script; it processes all `.pkl` files automatically.
- Generated plots are saved in respective folders under `convergence_history/`.

---

## Notes

- Skips normalization if mean of first 5 values is zero to avoid division by zero.
- Handles missing expected data gracefully with warnings.
- Designed for monitoring iterative CFD runs with turbulence modeling.