# Post-Processing of Slice Data for Parametric Study of Langtry Menter SST Model

**File:** `postproc_param_slices_models.py`

## Purpose

Analyze and compare aerodynamic surface data (e.g., pressure coefficient `Cp`, skin friction `Cf`) from CFD simulations using the **Langtry Menter SST** turbulence model. The script processes slice data for multiple parameter cases (varying Mach, Reynolds, AoA, turbulence intensity, and eddy viscosity ratio), and generates detailed and grouped plots of flow variables along the airfoil chord.

---

## Inputs

- Slice data files in the directory:  

sim/sd7003_L2/Langtry Menter SST//*_slices.dat

where `<case>` encodes Mach, Reynolds, angle of attack, TI, and EVR (e.g., `M0.10_Re60000_alpha4_TI0.0100_EVR1.0000`).

- Required columns in the data:
- `XoC`: chordwise position (index 3)
- `YoC`: vertical position (index 4)
- `CoefPressure`: pressure coefficient, `Cp` (index 9)
- `SkinFrictionX`: skin friction coefficient, `Cf` (index 11)

---

## Processing Steps

1. **Parse All Cases:**
 - Identify common cases in all models (only `Langtry Menter SST` is used).
 - Load slice files and group data by spanwise position `z`.

2. **Surface Filtering:**
 - Separate upper and lower surfaces based on `YoC`.
 - Filter points to ensure `0 ≤ X/C ≤ 1`.

3. **Plotting per Case and Slice:**
 - Generate `Cp` vs `X/C` plots for each case and spanwise slice.
 - Generate `Cf` vs `X/C` plots and additional variable plots if defined.

4. **Grouped Plots Across All Cases:**
 - Overlay multiple cases on a single plot per variable for fixed slice (`z`).
 - Save grouped plots and serialized figures as `.pkl`.

5. **Metadata Extraction:**
 - Automatically extracts Mach, Reynolds, AoA, TI, and EVR from case folder names.

---

## Outputs

- Plots saved in:

slices_plots/sd7003_L2//           # individual case plots
slices_plots/sd7003_L2/grouped_plots/    # grouped across cases
slices_plots/sd7003_L2/all_models_all_cases/  # combined by z-slice

- Filenames include variable names and `z` position (e.g., `Cp_vs_XoC_z0.5_all_cases.png`).

- Pickled matplotlib figures for reproducibility.

---

## Dependencies

- `numpy`
- `matplotlib`
- `pickle`
- `os`, `glob`, `re`, `itertools`

---

## Notes

- Only `Langtry Menter SST` is currently processed.
- Assumes `.dat` files contain formatted zone headers like `z = 0.5`.
- Handles multiple spanwise slices per case.
- Marker styles and LaTeX rendering are used for clear, publication-quality plots.

---

## Usage

Ensure that simulation data is located in `sim/sd7003_L2/Langtry Menter SST/<case>/`, with each case containing a `*_slices.dat` file. Then run:

```bash
python postproc_param_slices_models.py

Plots will be saved in structured subfolders under slices_plots/.

