# comparison of LM-SST, SST and laminar sim for different cases of alpha, mach and reynolds

import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pickle

# Get the directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define subdirectories containing simulation data
subdirs = ['sim/sd7003_L2']

# Directory to save plots
slices_plots_dir = os.path.join(script_dir, 'slices_plots')
os.makedirs(slices_plots_dir, exist_ok=True)

# Variable column indices
variables_index = {
   # 'YoC': 4,
   # 'ZoC': 5,
   # 'VelocityX': 6,
   # 'VelocityY': 7,
   # 'VelocityZ': 8,
    'CoefPressure': 9,
   # 'Mach': 10,
    'SkinFrictionX': 11
}

plt.rcParams.update({
    "text.usetex": True,
    "font.size": 14,          # Default font size
    "axes.labelsize": 12,     # Axis label size
    "axes.titlesize": 14,     # Title size
    "xtick.labelsize": 12,    # X tick label
    "ytick.labelsize": 12     # Y tick label
})

import glob
import itertools
line_styles = ['-', '-', '-.', ':']
markers = ['o', 's', '^', 'd', 'x', '*', '+']
style_cycle = itertools.cycle([(ls, m) for ls in line_styles for m in markers])

# Generalized model/case loop, using intersection of cases for all models
models = ["Langtry Menter SST"]
base_case_dir = os.path.join(script_dir, subdirs[0])

# Use intersection of case names available in all models
case_sets = [
    set(
        entry for entry in os.listdir(os.path.join(base_case_dir, model))
        if os.path.isdir(os.path.join(base_case_dir, model, entry))
    )
    for model in models
]
case_names = sorted(set.intersection(*case_sets))

all_case_data = {}

for case in case_names:
    sim_folder_name = os.path.basename(subdirs[0])
    sim_save_dir = os.path.join(slices_plots_dir, sim_folder_name)
    os.makedirs(sim_save_dir, exist_ok=True)
    case_save_dir = os.path.join(sim_save_dir, case)
    os.makedirs(case_save_dir, exist_ok=True)

    model_paths = {
        "Langtry Menter SST": os.path.join(base_case_dir, "Langtry Menter SST", case)
    }

    # Look for *_slices.dat file in each model/case directory
    data_files = {}
    for model, path in model_paths.items():
        match = glob.glob(os.path.join(path, '*_slices.dat'))
        if match:
            data_files[model] = match[0]
        else:
            data_files[model] = None

    if not all(os.path.exists(path) for path in model_paths.values()):
        print(f"Missing data files for {case}, skipping...")
        continue

    langtry_path = data_files["Langtry Menter SST"]

    if not langtry_path:
        print(f"Missing data files for {case}, skipping...")
        continue

    print(f"Processing {langtry_path}")

    def read_slices(file_path):
        slices = {}
        current_z = None
        data = []
        
        with open(file_path, 'r') as f:
            for line in f:
                zone_match = re.search(r'z\s*=\s*([\d.]+)', line)
                if zone_match:
                    if current_z is not None and data:
                        slices[current_z] = np.array(data)
                    current_z = float(zone_match.group(1))
                    data = []
                    continue

                if current_z is not None:
                    try:
                        values = list(map(float, line.split()))
                        if len(values) >= max(variables_index.values()) + 1:
                            data.append(values)
                    except ValueError:
                        continue
        
        if current_z is not None and data:
            slices[current_z] = np.array(data)
        
        return slices

    slices_langtry = read_slices(langtry_path)

    all_case_data[case] = {
        "langtry": slices_langtry,
        "case_save_dir": case_save_dir
    }

for z_value in next(iter(all_case_data.values()))["langtry"].keys():
    for case, data in all_case_data.items():
        slices_langtry = data["langtry"]

        # Define Y/C thresholds (user-controlled)
        upper_min_yoc = 0
        upper_max_yoc = 1
        lower_min_yoc = -1
        lower_max_yoc = 0

        # Define X/C thresholds (user-controlled)
        upper_min_xoc = 0
        upper_max_xoc = 1-0
        lower_min_xoc = 0
        lower_max_xoc = 1-0

        # Filter surfaces based on thresholds
        upper_surface_langtry = slices_langtry[z_value][
            (slices_langtry[z_value][:, 4] >= upper_min_yoc) & (slices_langtry[z_value][:, 4] <= upper_max_yoc) &
            (slices_langtry[z_value][:, 3] >= upper_min_xoc) & (slices_langtry[z_value][:, 3] <= upper_max_xoc)
        ]
        lower_surface_langtry = slices_langtry[z_value][
            (slices_langtry[z_value][:, 4] >= lower_min_yoc) & (slices_langtry[z_value][:, 4] <= lower_max_yoc) &
            (slices_langtry[z_value][:, 3] >= lower_min_xoc) & (slices_langtry[z_value][:, 3] <= lower_max_xoc)
        ]

        var_name_label = '$C_p$'

        try:
            parts = case.split('_')
            mach = float(parts[0][1:])
            reynolds = float(parts[1][2:])
            alpha = float(parts[2][5:])
            turb_intensity = float(parts[3][2:])
            eddy_vis_ratio = float(parts[4][3:])
        except Exception:
            mach = reynolds = alpha = turb_intensity = eddy_vis_ratio = 0.0

        plt.figure(figsize=(8, 6))
        combined_surface = np.vstack([upper_surface_langtry, lower_surface_langtry])
        plt.scatter(combined_surface[:, 3], combined_surface[:, variables_index["CoefPressure"]], color='r', marker='o', s=5, label=r'$k\mathrm{-}\omega\mathrm{-}\gamma\mathrm{-}Re_{\theta_t}$')

        plt.title(rf'{sim_folder_name}: {var_name_label} vs $X/C$ | $M={mach:.2f}$, $Re={reynolds:.2e}$, $\alpha={alpha:.0f}^\circ$, TI={turb_intensity:.4f}, EVR={eddy_vis_ratio:.4f}')
        plt.xlabel(r'$X/C$')
        plt.ylabel(var_name_label)
        plt.legend()
        plt.grid(True)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        filename = f'{sim_folder_name}_Cp_vs_XoC.png'
        plt.savefig(os.path.join(data["case_save_dir"], filename), dpi=600)
        plt.close()

for z_value in next(iter(all_case_data.values()))["langtry"].keys():
    for case, data in all_case_data.items():
        slices_langtry = data["langtry"]

        # Define Y/C thresholds (user-controlled)
        upper_min_yoc = 0
        upper_max_yoc = 1
        lower_min_yoc = -1
        lower_max_yoc = 0

        # Define X/C thresholds (user-controlled)
        upper_min_xoc = 0
        upper_max_xoc = 1-0
        lower_min_xoc = 0
        lower_max_xoc = 1-0

        # Filter surfaces based on thresholds
        upper_surface_langtry = slices_langtry[z_value][
            (slices_langtry[z_value][:, 4] >= upper_min_yoc) & (slices_langtry[z_value][:, 4] <= upper_max_yoc) &
            (slices_langtry[z_value][:, 3] >= upper_min_xoc) & (slices_langtry[z_value][:, 3] <= upper_max_xoc)
        ]
        lower_surface_langtry = slices_langtry[z_value][
            (slices_langtry[z_value][:, 4] >= lower_min_yoc) & (slices_langtry[z_value][:, 4] <= lower_max_yoc) &
            (slices_langtry[z_value][:, 3] >= lower_min_xoc) & (slices_langtry[z_value][:, 3] <= lower_max_xoc)
        ]

        for var_name, col_idx in variables_index.items():
            if var_name == "CoefPressure":
                continue

            var_data_langtry_upper = upper_surface_langtry[:, col_idx]
            var_data_langtry_lower = lower_surface_langtry[:, col_idx]

            plt.figure(figsize=(8, 6))
            combined_surface = np.vstack([upper_surface_langtry, lower_surface_langtry])
            plt.scatter(combined_surface[:, 3], combined_surface[:, col_idx], color='r', marker='o', s=5, label=r'$k\mathrm{-}\omega\mathrm{-}\gamma\mathrm{-}Re_{\theta_t}$')

            try:
                parts = case.split('_')
                mach = float(parts[0][1:])
                reynolds = float(parts[1][2:])
                alpha = float(parts[2][5:])
                turb_intensity = float(parts[3][2:])
                eddy_vis_ratio = float(parts[4][3:])
            except Exception:
                mach = reynolds = alpha = turb_intensity = eddy_vis_ratio = 0.0
            case_label = f'M{mach:.2f}_Re{int(reynolds)}_alpha{alpha:.0f}_TI{turb_intensity:.4f}_EVR{eddy_vis_ratio:.4f}'

            var_name_label = var_name

            if var_name in ['YoC', 'ZoC']:  # Adjust for Y/C vs X/C plot
                plt.axis('equal')

            plt.title(rf'{sim_folder_name}: {var_name_label} vs $X/C$  | $M={mach:.2f}$, $Re={reynolds:.2e}$, $\alpha={alpha:.0f}^\circ$, TI={turb_intensity:.4f}, EVR={eddy_vis_ratio:.4f}')
            plt.xlabel(r'$X/C$')
            plt.ylabel(var_name_label)
            plt.legend()
            plt.grid(True)
            #plt.axis('equal')
            plt.tight_layout()
            filename = f'{sim_folder_name}_{var_name}_vs_XoC.png'
            plt.savefig(os.path.join(data["case_save_dir"], filename), dpi=600)
            plt.close()

# Grouped plots per variable across all cases
grouped_plot_dir = os.path.join(slices_plots_dir, sim_folder_name, 'grouped_plots')
os.makedirs(grouped_plot_dir, exist_ok=True)

unique_markers = itertools.cycle(['o', 's', '^', 'd', 'x', '*', '+', 'v', '<', '>', 'p', 'h'])

for var_name, col_idx in variables_index.items():
    plt.figure(figsize=(8, 6))
    var_name_label = var_name
    if var_name == "CoefPressure":
        var_name_label = r"$C_p$"
    elif var_name == "SkinFrictionX":
        var_name_label = r"$C_f$"
    for case, data in all_case_data.items():
        marker = next(unique_markers)
        slices_langtry = data["langtry"]
        try:
            parts = case.split('_')
            mach = float(parts[0][1:])
            reynolds = float(parts[1][2:])
            alpha = float(parts[2][5:])
            turb_intensity = float(parts[3][2:])
            eddy_vis_ratio = float(parts[4][3:]) # label=r'$k\mathrm{-}\omega\mathrm{-}\gamma\mathrm{-}Re_{\theta_t}$'
        except Exception:
            mach = reynolds = alpha = turb_intensity = eddy_vis_ratio = 0.0
        label = (
    rf"$k\mathrm{{-}}\omega\mathrm{{-}}\gamma\mathrm{{-}}Re_{{\theta_t}},\ "
    rf"TI_\infty={turb_intensity:.4f},\ "
    rf"\left(\mu_t/\mu\right)_\infty={eddy_vis_ratio:.4f}$"
)
        for z_value in slices_langtry:
            upper_surface = slices_langtry[z_value][
                (slices_langtry[z_value][:, 4] >= 0) & (slices_langtry[z_value][:, 4] <= 1 ) &
                (slices_langtry[z_value][:, 3] >= 0) & (slices_langtry[z_value][:, 3] <= 1 )
            ]
            lower_surface = slices_langtry[z_value][
                (slices_langtry[z_value][:, 4] >= -1) & (slices_langtry[z_value][:, 4] < 0) &
                (slices_langtry[z_value][:, 3] >= 0) & (slices_langtry[z_value][:, 3] <= 1 )
            ]
            if var_name == "SkinFrictionX":
                plt.scatter(upper_surface[:, 3], upper_surface[:, col_idx], marker=marker, s=5 ,label=f'{label}')
                #plt.scatter(lower_surface[:, 3], lower_surface[:, col_idx], marker=marker, s=5 ,label=f'{label} (Lower Surface)')
                break  # Only use first z for SkinFrictionX
            else:
                combined_surface = np.vstack([upper_surface, lower_surface])
                plt.scatter(combined_surface[:, 3], combined_surface[:, col_idx], marker=marker, s=5, label=label)
    if var_name in ['YoC', 'ZoC']:
        plt.axis('equal')
    plt.title(f'{var_name_label} vs X/C (All Cases)')
    plt.xlabel('X/C')
    plt.ylabel(var_name_label)
    plt.legend( loc='best')
    plt.grid(True)
    if var_name == "CoefPressure":
        plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(grouped_plot_dir, f'{var_name}_vs_XoC_all_cases_grouped.png'), dpi=600)
    # Save figure as pickle file
    with open(os.path.join(grouped_plot_dir, f'{var_name}_vs_XoC_all_cases_grouped.pkl'), 'wb') as f:
        pickle.dump(plt.gcf(), f)
    plt.close()

# New combined plots for all cases and models per z_value and variable
all_models_all_cases_dir = os.path.join(slices_plots_dir, 'sd7003_L2', 'all_models_all_cases')
os.makedirs(all_models_all_cases_dir, exist_ok=True)

unique_markers = itertools.cycle(['o', 's', '^', 'd', 'x', '*', '+', 'v', '<', '>', 'p', 'h'])

for z_value in next(iter(all_case_data.values()))["langtry"].keys():
    plt.figure(figsize=(8, 6))
    for case, data in all_case_data.items():
        marker = next(unique_markers)
        try:
            parts = case.split('_')
            mach = float(parts[0][1:])
            reynolds = float(parts[1][2:])
            alpha = float(parts[2][5:])
            turb_intensity = float(parts[3][2:])
            eddy_vis_ratio = float(parts[4][3:])
        except Exception:
            mach = reynolds = alpha = turb_intensity = eddy_vis_ratio = 0.0
        case_label = f'M{mach:.2f}_Re{int(reynolds)}_alpha{alpha:.0f}_TI{turb_intensity:.4f}_EVR{eddy_vis_ratio:.4f}'
        for model_key, model_label in zip(["langtry"], ["Langtry Menter SST"]):
            model_label_fixed = r"$k\mathrm{-}\omega\mathrm{-}\gamma\mathrm{-}Re_{\theta_t}$"
            slices = data[model_key]
            if z_value not in slices:
                continue
            XoC = slices[z_value][:, 3]
            Cp = slices[z_value][:, variables_index["CoefPressure"]]
            ls, _ = next(style_cycle)
            plt.scatter(XoC, Cp, marker=marker, label=rf'{model_label_fixed} | {case_label}', s=5)
    plt.title(rf'{sim_folder_name}: $C_p$ vs $X/C$ at $z = {z_value}$ (All Cases \& Models)')
    plt.xlabel(r'$X/C$')
    plt.ylabel(r'$C_p$')
    plt.legend()
    plt.grid(True)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(all_models_all_cases_dir, f'{sim_folder_name}_Cp_vs_XoC_z{z_value}_all_cases.png'), dpi=600)
    plt.close()

    for var_name, col_idx in variables_index.items():
        if var_name == "CoefPressure":
            continue
        plt.figure(figsize=(8, 6))
        for case, data in all_case_data.items():
            marker = next(unique_markers)
            try:
                parts = case.split('_')
                mach = float(parts[0][1:])
                reynolds = float(parts[1][2:])
                alpha = float(parts[2][5:])
                turb_intensity = float(parts[3][2:])
                eddy_vis_ratio = float(parts[4][3:])
            except Exception:
                mach = reynolds = alpha = turb_intensity = eddy_vis_ratio = 0.0
            case_label = f'M{mach:.2f}_Re{int(reynolds)}_alpha{alpha:.0f}_TI{turb_intensity:.4f}_EVR{eddy_vis_ratio:.4f}'

            var_name_label = var_name

            for model_key, model_label in zip(["langtry"], ["Langtry Menter SST"]):
                model_label_fixed = r"$k\mathrm{-}\omega\mathrm{-}\gamma\mathrm{-}Re_{\theta_t}$"
                slices = data[model_key]
                if z_value not in slices:
                    continue
                XoC = slices[z_value][:, 3]
                Ydata = slices[z_value][:, col_idx]
                ls, _ = next(style_cycle)
                plt.scatter(XoC, Ydata, marker=marker, label=rf'{model_label_fixed} | {case_label}', s=5)
        if var_name in ['YoC', 'ZoC']:  # Adjust for Y/C vs X/C plot
            plt.axis('equal')
        plt.title(rf'{sim_folder_name}: {var_name_label} vs $X/C$ at $z = {z_value}$ (All Cases \& Models)')
        plt.xlabel(r'$X/C$')
        plt.ylabel(var_name_label)
        plt.legend()
        plt.grid(True)
        # plt.axis('equal')
        plt.tight_layout()
        filename = f'{sim_folder_name}_{var_name}_vs_XoC_z{z_value}_all_cases.png'
        plt.savefig(os.path.join(all_models_all_cases_dir, filename), dpi=600)
        plt.close()

print(f'Plots saved in: {slices_plots_dir}')