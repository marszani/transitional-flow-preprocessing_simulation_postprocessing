# Grid Convergence and Order of Convergence Analysis

**File:** `compute_order_of_convergence.py`

## Purpose

Calculate the achieved order of convergence and Grid Convergence Index (GCI) for aerodynamic coefficients (drag `Cd` and lift `Cl`) from CFD results at multiple mesh refinement levels. Produce plots showing the convergence behavior.

---

## Input Data

- The script reads from text files (e.g., `grid_refinement_sc20714.txt`) containing aerodynamic coefficients and grid info.
- Example input file for Eppler 387 airfoil: `grid_refinement_data.txt`
- Expected file contents include:
  - Lines with Python dict-like data specifying `cd` and `cl` values at levels L0, L1, L2.
  - A line specifying total number of cells for each grid level, e.g.:
    ```
    Total number of cells: 10000 40000 160000
    ```
  - A line specifying refinement factors, e.g.:
    ```
    refinement = [1.0, 2.0, 4.0]
    ```

---

## Main Steps

1. **Data Reading:**
   - Parse coefficient values, cell counts, and refinement from the input file.

2. **Grid Spacing Calculation:**
   - Calculate grid spacing \( h = N_{\text{cells}}^{-1/d} \) for \( d=2 \).
   - Compute refinement ratio \( r = h_{L2} / h_{L0} \).

3. **Order of Convergence:**
   - Calculate achieved order \( p \) using three grid levels.
   - Calculate Richardson extrapolated values.

4. **GCI Computation:**
   - Compute GCI percentages between grid levels.

5. **Plotting:**
   - Generate plots of \( C_d \), \( C_l \) vs. grid spacing.
   - Plot lift coefficient vs. refinement level.
   - Plot bar chart comparing GCI values.

---

## Functions

- `read_results_from_file(file_path)`: Parses input file for data.
- `achieved_order_of_convergence(f_L0, f_L1, f_L2, r)`: Computes order \( p \).
- `compute_gci(fine, coarse, r, p)`: Computes GCI percentage.

---

## Usage

- Specify input files in `file_names` list.
- Run script; outputs are printed convergence metrics and saved plots.
- Plots saved as PNG files in the script directory.

---

## Notes

- Assumes 2D grids.
- Input files must follow the expected data formatting.
- Requires `numpy` and `matplotlib`.