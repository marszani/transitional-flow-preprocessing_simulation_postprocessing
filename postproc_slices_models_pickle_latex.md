# Post-Processing and Comparison of Slice Data for Turbulence Models

**File:** `postproc_slices_models_pickle_latex.py`

## Purpose

Compare aerodynamic flow slice data between two turbulence models, **Langtry Menter SST** and **Menter SST**, for various cases defined by Mach number, Reynolds number, and angle of attack. Generate and save plots of pressure coefficient (`Cp`), skin friction coefficient (`Cf`), and other flow variables along the chord.

---

## Input Data

- CFD slice data files named `*_slices.dat` located in:

sim/e387_L0/Langtry Menter SST//
sim/e387_L0/Menter SST//

- Each slice file contains multiple spanwise slices indexed by `z` coordinate.
- Each slice holds tabular flow data with columns for spatial coordinates and flow variables.

---

## Variables Extracted

| Variable Name    | Column Index | Description               |
|------------------|--------------|---------------------------|
| `XoC`            | 3            | Chordwise normalized coordinate (x/c) |
| `YoC`            | 4            | Vertical normalized coordinate (y/c)  |
| `ZoC`            | 5            | Spanwise coordinate (z)              |
| `CoefPressure`   | 9            | Pressure coefficient (`Cp`)          |
| `SkinFrictionX`  | 11           | Skin friction coefficient (`Cf`)    |
| Others (Velocities, Mach, etc.) | Varies | Additional flow quantities          |

---

## Workflow

1. **Case Discovery:**  
 Detect cases present in both models’ directories.

2. **Data Loading:**  
 Parse slice files into numpy arrays grouped by spanwise location (`z`).

3. **Surface Filtering:**  
 Separate upper and lower surfaces using vertical coordinate thresholds:
 - Upper surface: `YoC` in [0, 1]
 - Lower surface: `YoC` in [-1, 0]
 Also limit chordwise coordinate (`XoC`) to [0, 1].

4. **Plotting:**  
 - For each case and spanwise slice:
   - Plot `Cp` vs. `x/c` for both models on the same plot.
   - Plot other variables (e.g., `Cf`) similarly.
   - Markers and colors differentiate models and surfaces.
   - Titles include Mach, Reynolds number, and angle of attack extracted from case name.
   - Invert y-axis for `Cp` plots to follow aerodynamic convention.

5. **Save Outputs:**  
 - Save plots as PDF files.
 - Serialize matplotlib figure objects as pickle files (`.fig.pickle`) alongside PDFs.

6. **Combined Plots:**  
 Generate aggregate plots of all cases and models per spanwise location and variable for overview comparison.

---

## Output

- Saved in:

slices_plots/e387_L0//

- Combined plots saved under:

slices_plots/e387_L0/all_models_all_cases/

---

## Dependencies

- Python standard libraries: `os`, `re`, `glob`, `pickle`, `itertools`
- Third-party: `numpy`, `matplotlib`

---

## Notes

- The filtering thresholds for upper and lower surfaces and chordwise limits are configurable.
- The script automatically skips cases or data files if missing or incomplete.
- Designed for the Eppler 387 airfoil dataset but easily adaptable to other cases or turbulence models.

---

## Usage

Run the script in a directory where the simulation results are organized as above. The script processes all common cases and produces detailed comparative plots for aerodynamic slice data.

