# Plotting Pressure and Skin Friction Distributions at Different Grid Levels

**File:** `plot_L_comparison.py`

## Purpose

Load CFD slice data for selected aerodynamic cases at multiple grid refinement levels, plot pressure coefficient (`Cp`) and skin friction coefficient (`Cf`) distributions along the chord, and save both plots and serialized data for later analysis.

---

## Input Data

- Reads slice data files (`*_slices.dat`) organized in subdirectories by grid level:
  - `sim/sd7003_L0`
  - `sim/sd7003_L1`
  - `sim/sd7003_L2`
- Data files contain multiple slices at different spanwise positions (`z`).
- Each slice has tabular data with columns including spatial coordinates and flow quantities.

---

## Configuration

- **Models:** Currently only supports `"Langtry Menter SST"`.
- **Variables plotted:**
  - `CoefPressure` (pressure coefficient, `Cp`) at column index 9
  - `SkinFrictionX` (skin friction coefficient, `Cf`) at column index 11
- **Output Directory:** `grouped_plots` (created if it does not exist).

---

## Processing Steps

1. **Data Loading:**  
   For each subdirectory and case folder, load slice data grouped by spanwise location `z`.

2. **Data Filtering and Sorting:**  
   Separate data points above (`y >= 0`) and below (`y < 0`) the chord line to plot upper and lower surfaces. Data filtered for \( 0 \leq x/c \leq 1 \).

3. **Plotting:**
   - Plot `Cp` or `Cf` vs. normalized chord position (`x/c`).
   - Use different markers for each grid level (`L0`, `L1`, `L2`).
   - Invert y-axis for pressure coefficient plots to follow aerodynamic convention.
   - Add descriptive titles with flow parameters extracted from case names (Mach number, Reynolds number, angle of attack, turbulence intensity, eddy viscosity ratio).
   - Add legend specifying grid level.

4. **Data Export:**
   - Save plotted data (`x`, `y` values) along with labels and titles into `.pkl` files inside `grouped_plots`.
   - Save plot figures as PNG files.

---

## Utility Functions

- `read_slices(file_path)`: Parses slice data grouped by spanwise position `z` from slice data files.
- `format_case_title(case_str, var)`: Extracts flow parameters from case name string and formats a LaTeX title string for plots.

---

## Plot Styling

- Uses LaTeX for text rendering.
- Marker styles cycle through a predefined list for visual distinction.
- Grid and axis labels styled for clarity.

---

## Usage

- Place CFD slice data in the prescribed folder structure.
- Run the script; it will automatically find and process common cases across all grid levels.
- Resulting plots and `.pkl` files will be saved under the `grouped_plots` directory.

---

## Notes

- Only cases present in all three grid level folders (`L0`, `L1`, `L2`) are plotted.
- Handles missing or malformed data gracefully.
- Designed for analysis of transitional turbulence model results on SD7003 airfoil.