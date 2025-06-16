# Comparison of Turbulence Models via Slice Data Plots

**File:** `postproc_slices_models.py`

## Purpose

Compare and visualize aerodynamic slice data between two turbulence models, **Langtry Menter SST** and **Menter SST**, across multiple simulation cases varying in Mach number, Reynolds number, and angle of attack. The script plots pressure coefficient (`Cp`), skin friction coefficient (`Cf`), and other flow variables along the chord at various spanwise locations.

---

## Input Data

- Slice data files named `*_slices.dat` located in the directory structure:

sim/e387_L0/Langtry Menter SST//
sim/e387_L0/Menter SST//

- Each slice file contains data organized in multiple slices at discrete spanwise locations (`z`).
- Data columns include spatial coordinates (`XoC`, `YoC`, `ZoC`) and flow variables (e.g., `CoefPressure`, `SkinFrictionX`).

---

## Workflow

1. **Case Discovery:**  
 - Detect simulation cases present in both turbulence models’ directories to ensure fair comparison.

2. **Data Loading:**  
 - Read slice files for each case and model, parsing data into numpy arrays keyed by spanwise `z` values.

3. **Data Filtering:**  
 - Separate upper and lower surface data based on vertical coordinate thresholds (`YoC`):
   - Upper surface: `YoC` in [0, 1]
   - Lower surface: `YoC` in [-1, 0]
 - Limit chordwise coordinate (`XoC`) to [0, 1] for plotting.

4. **Plot Generation:**  
 - For each case and spanwise slice:
   - Plot `Cp` versus `XoC` for both models on the same figure.
   - Similarly plot other variables such as `SkinFrictionX`, velocities, and Mach number.
   - Use different markers and colors to distinguish models and surfaces.
   - Titles include Mach, Reynolds number, and angle of attack parsed from the case name.
   - Invert y-axis for `Cp` plots following aerodynamic convention.

5. **Combined Plots:**  
 - Aggregate plots comparing all cases and both models at each spanwise location for selected variables.

6. **Saving Outputs:**  
 - Save all plots as high-resolution PNG files in the folder:
   ```
   slices_plots/e387_L0/<case>/
   ```
 - Combined plots are saved in:
   ```
   slices_plots/e387_L0/all_models_all_cases/
   ```

---

## Variables and Column Indices

| Variable        | Column Index | Description                  |
|-----------------|--------------|------------------------------|
| `YoC`           | 4            | Vertical coordinate (y/c)    |
| `ZoC`           | 5            | Spanwise coordinate (z)      |
| `VelocityX`     | 6            | Velocity in x-direction      |
| `VelocityY`     | 7            | Velocity in y-direction      |
| `VelocityZ`     | 8            | Velocity in z-direction      |
| `CoefPressure`  | 9            | Pressure coefficient (`Cp`)  |
| `Mach`          | 10           | Local Mach number            |
| `SkinFrictionX` | 11           | Skin friction coefficient (`Cf`) |

---

## Dependencies

- Python standard libraries: `os`, `re`, `glob`, `itertools`
- Third-party libraries: `numpy`, `matplotlib`

---

## Notes

- The script uses LaTeX rendering for plot labels and titles.
- Handles missing data gracefully by skipping incomplete cases.
- User-configurable thresholds control surface and chordwise filtering.
- Suitable for detailed comparison of transitional and baseline turbulence model results.

---

## Usage

Run the script in the directory containing the `sim/e387_L0` folder structured by turbulence model and case. The script automatically processes all valid cases and saves comparison plots.
