import numpy as np
import matplotlib.pyplot as plt
import os

def read_results_from_file(file_path):
    cd_values = {}
    cl_values = {}
    cell_counts = {}
    refinement = []

    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                # Read Python dict-like lines safely
                if line.startswith("{") and line.endswith("}"):
                    try:
                        data_dict = eval(line)  # Assumes file is trusted
                        for key, value in data_dict.items():
                            if 'cd' in key:
                                if 'L0' in key:
                                    cd_values["L0"] = float(value)
                                elif 'L1' in key:
                                    cd_values["L1"] = float(value)
                                elif 'L2' in key:
                                    cd_values["L2"] = float(value)
                            elif 'cl' in key:
                                if 'L0' in key:
                                    cl_values["L0"] = float(value)
                                elif 'L1' in key:
                                    cl_values["L1"] = float(value)
                                elif 'L2' in key:
                                    cl_values["L2"] = float(value)
                    except Exception as e:
                        print(f"Error parsing line as dict: {e}")
                elif 'Total number of cells' in line:
                    cell_counts_list = line.split(":")[1].strip().split()
                    cell_counts["L0"], cell_counts["L1"], cell_counts["L2"] = map(int, cell_counts_list)
                elif 'refinement' in line:
                    refinement = list(map(float, line.split('=')[1].strip().strip('[]').split(',')))
    except Exception as e:
        print(f"Error reading file: {e}")
    
    return cd_values, cl_values, cell_counts, refinement

# File names
file_names = [
    "grid_refinement_sc20714.txt"
]

# Loop over each case
for file_name in file_names:
    print(f"\nProcessing: {file_name}")

    current_directory = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(current_directory, file_name)
    case_id = r"SD7003, M=0.1, Re=60000, $\alpha=4^\circ$"

    cd_values, cl_values, cell_counts, refinement = read_results_from_file(file_path)

    print("Refinement levels:", refinement)
    print("CD values:", cd_values)
    print("CL values:", cl_values)
    print("Number of cells:", cell_counts)

    d = 2
    h = {key: count**(-1/d) for key, count in cell_counts.items()}
    r = h["L2"] / h["L0"]
    print(f"Grid refinement ratio (r): {r:.4f}")

    def achieved_order_of_convergence(f_L0, f_L1, f_L2, r):
        numerator = (f_L2 - f_L1)
        denominator = (f_L1 - f_L0)
        if numerator * denominator <= 0:
            return np.nan
        else:
            return np.log(abs(numerator / denominator)) / np.log(r)

    def compute_gci(fine, coarse, r, p):
        Fs = 1.0
        if fine < coarse:
            fine, coarse = coarse, fine
        return (Fs * abs(fine - coarse) / (r**p - 1) / fine) * 100

    # ---------- Cd analysis ----------
    try:
        p_cd = achieved_order_of_convergence(cd_values["L0"], cd_values["L1"], cd_values["L2"], r)
        cd_h0 = cd_values["L0"] + (cd_values["L0"] - cd_values["L1"]) / (r**p_cd - 1)
        print(f"Achieved order (p) for Cd: {p_cd:.4f}")
        print(f"Richardson extrapolated Cd: {cd_h0:.8f}")

        gci_L0_L1 = compute_gci(cd_values["L1"], cd_values["L0"], r, p_cd)
        gci_L1_L2 = compute_gci(cd_values["L2"], cd_values["L1"], r, p_cd)
        print(f"GCI Cd L0-L1: {gci_L0_L1:.4f}%")
        print(f"GCI Cd L1-L2: {gci_L1_L2:.4f}%")
    except Exception as e:
        print(f"Error in Cd convergence calculation: {e}")

    # ---------- Cl analysis ----------
    try:
        p_cl = achieved_order_of_convergence(cl_values["L0"], cl_values["L1"], cl_values["L2"], r)
        cl_h0 = cl_values["L0"] + (cl_values["L0"] - cl_values["L1"]) / (r**p_cl - 1)
        print(f"Achieved order (p) for Cl: {p_cl:.4f}")
        print(f"Richardson extrapolated Cl: {cl_h0:.8f}")

        gci_cl_L0_L1 = compute_gci(cl_values["L1"], cl_values["L0"], r, p_cl)
        gci_cl_L1_L2 = compute_gci(cl_values["L2"], cl_values["L1"], r, p_cl)
        print(f"GCI Cl L0-L1: {gci_cl_L0_L1:.4f}%")
        print(f"GCI Cl L1-L2: {gci_cl_L1_L2:.4f}%")
    except Exception as e:
        print(f"Error in Cl convergence calculation: {e}")

    # ---------- Plotting ----------
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')

    levels = ["L0", "L1", "L2"]
    colors = {"L0": "#1f77b4", "L1": "#ff7f0e", "L2": "#2ca02c", "Richardson": "#d62728"}
    h_values = [h[level] for level in levels]

    # Cd vs h
    cd_values_list = [cd_values[level] for level in levels]
    plt.figure(figsize=(8, 6))
    plt.plot(h_values, cd_values_list, 'o-', label="Computed $C_d$", color=colors["L0"])
    plt.plot(0, cd_h0, 's', markersize=8, color=colors["Richardson"], label="Richardson Extrapolated $C_d$")
    plt.xlabel(r"$h$", fontsize=12)
    plt.ylabel(r"$C_d$", fontsize=12)
    plt.title(fr"Grid Convergence for $C_d$: {case_id}", fontsize=14)
    plt.grid(True, which="both", linestyle="--")
    plt.legend(fontsize=12)
    plt.savefig(os.path.join(current_directory, f'grid_convergence_Cd_{case_id}.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Cl vs h
    cl_values_list = [cl_values[level] for level in levels]
    plt.figure(figsize=(8, 6))
    plt.plot(h_values, cl_values_list, 'o-', label="Computed $C_l$", color=colors["L1"])
    plt.plot(0, cl_h0, 's', markersize=8, color=colors["Richardson"], label="Richardson Extrapolated $C_l$")
    plt.xlabel(r"$h$", fontsize=12)
    plt.ylabel(r"$C_l$", fontsize=12)
    plt.title(fr"Grid Convergence for $C_l$: {case_id}", fontsize=14)
    plt.grid(True, which="both", linestyle="--")
    plt.legend(fontsize=12)
    plt.savefig(os.path.join(current_directory, f'grid_convergence_Cl_{case_id}.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Cl vs refinement level
    levels_numeric = [0, 1, 2]
    plt.plot(levels_numeric, cl_values_list, 'o-', color=colors["L2"], label="Computed $C_l$")
    plt.scatter([-1], [cl_h0], s=100, marker='s', color=colors["Richardson"], label="Richardson Extrapolated $C_l$")
    plt.xticks(levels_numeric, levels)
    plt.xlabel(r"Grid Refinement Level", fontsize=12)
    plt.ylabel(r"$C_l$", fontsize=12)
    plt.title(fr"Lift Coefficient vs. Refinement: {case_id}", fontsize=14)
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.savefig(os.path.join(current_directory, f'lift_vs_refinement_{case_id}.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # GCI bar plot for Cl
    plt.figure(figsize=(8, 6))
    plt.bar(["L0 - L1", "L1 - L2"], [gci_cl_L0_L1, gci_cl_L1_L2], color=[colors["L0"], colors["L2"]])
    plt.xlabel(r"Grid Level Comparison", fontsize=12)
    plt.ylabel(r"GCI ($\%$)", fontsize=12)
    plt.title(fr"GCI Comparison (Cl): {case_id}", fontsize=14)
    plt.grid(axis='y', linestyle='--')
    plt.savefig(os.path.join(current_directory, f'gci_comparison_Cl_{case_id}.png'), dpi=300, bbox_inches='tight')
    plt.close()