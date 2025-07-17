# comparison of LM-SST, SST and laminar sim for different cases of alpha, mach and reynolds


import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pickle

# Get the directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define subdirectories containing simulation data
subdirs = ['sim/e387_L0']

plot_title_prefix = os.path.basename(subdirs[0]).replace('_L0', '').upper()

# Directory to save plots
slices_plots_dir = os.path.join(script_dir, 'slices_plots')
os.makedirs(slices_plots_dir, exist_ok=True)

# Variable column indices
variables_index = {
    'YoC': 4,
    'ZoC': 5,
    'VelocityX': 6,
    'VelocityY': 7,
    'VelocityZ': 8,
    'CoefPressure': 9,
    'Mach': 10,
    'SkinFrictionX': 11
}

# Enable LaTeX for plotting
plt.rcParams['text.usetex'] = True

import glob
import itertools
line_styles = ['-', '-', '-.', ':']
markers = ['o', 's', '^', 'd', 'x', '*', '+']
style_cycle = itertools.cycle([(ls, m) for ls in line_styles for m in markers])

models = ["Langtry Menter SST", "Menter SST"]
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
        model: os.path.join(base_case_dir, model, case)
        for model in models
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
    menter_path = data_files["Menter SST"]

    if not (langtry_path and menter_path):
        print(f"Missing data files for {case}, skipping...")
        continue

    print(f"Processing {langtry_path}, {menter_path}")

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
    slices_menter = read_slices(menter_path)

    all_case_data[case] = {
        "langtry": slices_langtry,
        "menter": slices_menter,
        "case_save_dir": case_save_dir
    }

for z_value in next(iter(all_case_data.values()))["menter"].keys():
    for case, data in all_case_data.items():
        slices_menter = data["menter"]
        slices_langtry = data["langtry"]

        # Define Y/C thresholds (user-controlled)
        upper_min_yoc = 0
        upper_max_yoc = 1
        lower_min_yoc = -1
        lower_max_yoc = 0

        # Define X/C thresholds (user-controlled)
        upper_min_xoc = 1e-9
        upper_max_xoc = 1-1e-9
        lower_min_xoc = 1e-9
        lower_max_xoc = 1-1e-9

        # Filter surfaces based on thresholds
        upper_surface_menter = slices_menter[z_value][
            (slices_menter[z_value][:, 4] >= upper_min_yoc) & (slices_menter[z_value][:, 4] <= upper_max_yoc) &
            (slices_menter[z_value][:, 3] >= upper_min_xoc) & (slices_menter[z_value][:, 3] <= upper_max_xoc)
        ]
        lower_surface_menter = slices_menter[z_value][
            (slices_menter[z_value][:, 4] >= lower_min_yoc) & (slices_menter[z_value][:, 4] <= lower_max_yoc) &
            (slices_menter[z_value][:, 3] >= lower_min_xoc) & (slices_menter[z_value][:, 3] <= lower_max_xoc)
        ]
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
            mach = float(case.split('_')[0][1:])
            reynolds = float(case.split('_')[1][2:])
            alpha = float(case.split('_')[2][5:])
        except Exception:
            mach = 0.0
            reynolds = 0.0
            alpha = 0.0

        plt.figure(figsize=(6, 6))

        # Combine upper and lower surfaces for each model
        for model_key, model_label in zip(["menter", "langtry"], [r'$k$-$\omega$ SST', r'$k$-$\omega$-$\gamma$-$Re_{\theta_t}$']):
            if model_key == "menter":
                combined_surface = np.vstack((upper_surface_menter, lower_surface_menter))
            else:
                combined_surface = np.vstack((upper_surface_langtry, lower_surface_langtry))
            
            # Sort combined_surface by x-coordinate (column 3)
            sorted_indices = np.argsort(combined_surface[:, 3])
            combined_surface_sorted = combined_surface[sorted_indices]

            x = combined_surface_sorted[:, 3]
            y = combined_surface_sorted[:, variables_index["CoefPressure"]]

            plt.plot(
                x,
                y,
                label=rf'{model_label}',
                color='b' if model_key == "menter" else 'r',
                marker='o',
                linestyle='None',
                markersize=2,
                linewidth=1
            )

        plt.title(rf'{plot_title_prefix}: {var_name_label} vs $X/C$ | $M={mach:.2f}$, $Re={reynolds:.2e}$, $\alpha={alpha:.0f}^\circ$')
        plt.xlabel(r'$X/C$')
        plt.ylabel(var_name_label)
        plt.legend()
        plt.grid(True)
        plt.gca().invert_yaxis()
        plt.tight_layout()
        filename = f'{sim_folder_name}_Cp_vs_XoC.pdf'
        with open(os.path.join(data["case_save_dir"], filename.replace('.pdf', '.fig.pickle')), 'wb') as f:
            pickle.dump(plt.gcf(), f)
        plt.close()

for z_value in next(iter(all_case_data.values()))["menter"].keys():
    for case, data in all_case_data.items():
        slices_menter = data["menter"]
        slices_langtry = data["langtry"]

        # Define Y/C thresholds (user-controlled)
        upper_min_yoc = 0
        upper_max_yoc = 1
        lower_min_yoc = -1
        lower_max_yoc = 0

        # Define X/C thresholds (user-controlled)
        upper_min_xoc = 1e-9
        upper_max_xoc = 1-1e-9
        lower_min_xoc = 1e-9
        lower_max_xoc = 1-1e-9

        # Filter surfaces based on thresholds
        upper_surface_menter = slices_menter[z_value][
            (slices_menter[z_value][:, 4] >= upper_min_yoc) & (slices_menter[z_value][:, 4] <= upper_max_yoc) &
            (slices_menter[z_value][:, 3] >= upper_min_xoc) & (slices_menter[z_value][:, 3] <= upper_max_xoc)
        ]
        lower_surface_menter = slices_menter[z_value][
            (slices_menter[z_value][:, 4] >= lower_min_yoc) & (slices_menter[z_value][:, 4] <= lower_max_yoc) &
            (slices_menter[z_value][:, 3] >= lower_min_xoc) & (slices_menter[z_value][:, 3] <= lower_max_xoc)
        ]
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

            var_data_menter_upper = upper_surface_menter[:, col_idx]
            var_data_langtry_upper = upper_surface_langtry[:, col_idx]
            var_data_menter_lower = lower_surface_menter[:, col_idx]
            var_data_langtry_lower = lower_surface_langtry[:, col_idx]

            plt.figure(figsize=(6, 6))
            for model_key, model_label in zip(
                ["menter", "langtry"],
                [r'$k$-$\omega$ SST', r'$k$-$\omega$-$\gamma$-$Re_{\theta_t}$']
            ):
                plt.plot(
                    upper_surface_menter[:, 3] if model_key == "menter" else upper_surface_langtry[:, 3],
                    var_data_menter_upper if model_key == "menter" else var_data_langtry_upper,
                    label=rf'{model_label} (Upper Surface)',
                    color='b' if model_key == "menter" else 'r',
                    marker='o',
                    linestyle='None',
                    markersize=2,
                    linewidth=1
                )
                plt.plot(
                    lower_surface_menter[:, 3] if model_key == "menter" else lower_surface_langtry[:, 3],
                    var_data_menter_lower if model_key == "menter" else var_data_langtry_lower,
                    label=rf'{model_label} (Lower Surface)',
                    color='b' if model_key == "menter" else 'r',
                    marker='x',
                    linestyle='None',
                    markersize=2,
                    linewidth=1
                )

            try:
                mach = float(case.split('_')[0][1:])
                reynolds = float(case.split('_')[1][2:])
                alpha = float(case.split('_')[2][5:])
            except Exception:
                mach = 0.0
                reynolds = 0.0
                alpha = 0.0
            case_label = f'M{mach:.2f}_Re{int(reynolds)}_alpha{alpha:.0f}'

            var_name_label = r'$C_f$' if var_name == 'SkinFrictionX' else var_name

            if var_name in ['YoC', 'ZoC']:  # Adjust for Y/C vs X/C plot
                plt.axis('equal')

            plt.title(rf'{plot_title_prefix}: {var_name_label} vs $X/C$  | $M={mach:.2f}$, $Re={reynolds:.2e}$, $\alpha={alpha:.0f}^\circ$')
            plt.xlabel(r'$X/C$')
            plt.ylabel(var_name_label)
            plt.legend()
            plt.grid(True)
            #plt.axis('equal')
            plt.tight_layout()
            filename = f'{sim_folder_name}_{var_name}_vs_XoC.pdf'
            with open(os.path.join(data["case_save_dir"], filename.replace('.pdf', '.fig.pickle')), 'wb') as f:
                pickle.dump(plt.gcf(), f)
            plt.close()

# New combined plots for all cases and models per z_value and variable
all_models_all_cases_dir = os.path.join(slices_plots_dir, 'e387_L0', 'all_models_all_cases')
os.makedirs(all_models_all_cases_dir, exist_ok=True)

for z_value in next(iter(all_case_data.values()))["menter"].keys():
    plt.figure(figsize=(6, 6))
    for case, data in all_case_data.items():
        mach = float(case.split('_')[0][1:])
        reynolds = float(case.split('_')[1][2:])
        alpha = float(case.split('_')[2][5:])
        case_label = f'M{mach:.2f}_Re{int(reynolds)}_alpha{alpha:.0f}'
        for model_key, model_label in zip(
            ["menter", "langtry"],
            [r'$k$-$\omega$ SST', r'$k$-$\omega$-$\gamma$-$Re_{\theta_t}$']
        ):
            slices = data[model_key]
            if z_value not in slices:
                continue
            XoC = slices[z_value][:, 3]
            Cp = slices[z_value][:, variables_index["CoefPressure"]]
            ls, marker = next(style_cycle)
            plt.plot(XoC, Cp, linestyle=ls, marker=marker, label=rf'{model_label} | {case_label}')
    plt.title(rf'{plot_title_prefix}: $C_p$ vs $X/C$ at $z = {z_value}$ (All Cases \& Models)')
    plt.xlabel(r'$X/C$')
    plt.ylabel(r'$C_p$')
    plt.legend()
    plt.grid(True)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    filename = f'{sim_folder_name}_Cp_vs_XoC_z{z_value}_all_cases.pdf'
    with open(os.path.join(all_models_all_cases_dir, filename.replace('.pdf', '.fig.pickle')), 'wb') as f:
        pickle.dump(plt.gcf(), f)
    plt.close()

    for var_name, col_idx in variables_index.items():
        if var_name == "CoefPressure":
            continue
        plt.figure(figsize=(6, 6))
        for case, data in all_case_data.items():
            mach = float(case.split('_')[0][1:])
            reynolds = float(case.split('_')[1][2:])
            alpha = float(case.split('_')[2][5:])
            case_label = f'M{mach:.2f}_Re{int(reynolds)}_alpha{alpha:.0f}'

            var_name_label = r'$C_f$' if var_name == 'SkinFrictionX' else var_name

            for model_key, model_label in zip(["menter", "langtry"], [r'$k$-$\omega$ SST', r'$k$-$\omega$-$\gamma$-$Re_{\theta_t}$']):
                slices = data[model_key]
                if z_value not in slices:
                    continue
                XoC = slices[z_value][:, 3]
                Ydata = slices[z_value][:, col_idx]
                ls, marker = next(style_cycle)
                plt.plot(XoC, Ydata, linestyle=ls, marker=marker, label=rf'{model_label} | {case_label}')
        if var_name in ['YoC', 'ZoC']:  # Adjust for Y/C vs X/C plot
            plt.axis('equal')
        plt.title(rf'{plot_title_prefix}: {var_name_label} vs $X/C$ at $z = {z_value}$ (All Cases \& Models)')
        plt.xlabel(r'$X/C$')
        plt.ylabel(var_name_label)
        plt.legend()
        plt.grid(True)
        # plt.axis('equal')
        plt.tight_layout()
        filename = f'{sim_folder_name}_{var_name}_vs_XoC_z{z_value}_all_cases.pdf'
        with open(os.path.join(all_models_all_cases_dir, filename.replace('.pdf', '.fig.pickle')), 'wb') as f:
            pickle.dump(plt.gcf(), f)
        plt.close()

print(f'Plots saved in: {slices_plots_dir}')