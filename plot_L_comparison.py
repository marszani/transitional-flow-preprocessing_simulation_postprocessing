import os
import re
import numpy as np
import matplotlib.pyplot as plt
import itertools
import glob
import pickle

script_dir = os.path.dirname(os.path.abspath(__file__))
subdirs = ['sim/sd7003_L0', 'sim/sd7003_L1', 'sim/sd7003_L2']
models = ['Langtry Menter SST']
variables = {'CoefPressure': 9, 'SkinFrictionX': 11}
output_dir = os.path.join(script_dir, 'grouped_plots')
os.makedirs(output_dir, exist_ok=True)

def read_slices(file_path):
    slices, current_z, data = {}, None, []
    with open(file_path, 'r') as f:
        for line in f:
            if 'z =' in line:
                if current_z is not None and data:
                    slices[current_z] = np.array(data, dtype=float)
                current_z = float(re.search(r'z\s*=\s*([\d.]+)', line).group(1))
                data = []
            else:
                parts = line.split()
                if len(parts) < 12:
                    continue
                try:
                    data.append([float(x) for x in parts])
                except ValueError:
                    continue
        if current_z is not None and data:
            slices[current_z] = np.array(data, dtype=float)
    return slices

def format_case_title(case_str, var):
    M = re.search(r'M(\d*\.?\d+)', case_str)
    Re = re.search(r'Re(\d+)', case_str)
    alpha = re.search(r'alpha(\d*\.?\d+)', case_str)
    TI = re.search(r'TI(\d*\.?\d+)', case_str)
    EVR = re.search(r'EVR(\d*\.?\d+)', case_str)

    M_val = M.group(1) if M else '?'
    Re_val = f"{int(Re.group(1)):.1e}" if Re else '?'
    alpha_val = alpha.group(1) if alpha else '?'
    TI_val = TI.group(1) if TI else '?'
    EVR_val = EVR.group(1) if EVR else '?'

    if var == "SkinFrictionX":
        title = (r"$C_f$ at "
                 rf"$M={M_val}$, "
                 rf"$Re={Re_val}$, "
                 rf"$\alpha={alpha_val}^\circ$, "
                 rf"$TI_{{\infty}}={TI_val}$, "
                 rf"$\mu_t/\mu_{{\infty}}={EVR_val}$")
    elif var == "CoefPressure":
        title = (r"$C_p$ at "
                 rf"$M={M_val}$, "
                 rf"$Re={Re_val}$, "
                 rf"$\alpha={alpha_val}^\circ$, "
                 rf"$TI_{{\infty}}={TI_val}$, "
                 rf"$\mu_t/\mu_{{\infty}}={EVR_val}$")
    else:
        title = case_str

    return title

plt.rcParams.update({
    "text.usetex": True,
    "font.size": 12
})

marker_cycle = itertools.cycle(['o', 's', '^'])

# --- Load all case data
case_data = {}
for subdir in subdirs:
    model_dir = os.path.join(script_dir, subdir, models[0])
    if not os.path.isdir(model_dir):
        print(f"Skipping: Missing model path {model_dir}")
        continue
    for case in os.listdir(model_dir):
        case_path = os.path.join(model_dir, case)
        match = glob.glob(os.path.join(case_path, '*_slices.dat'))
        if not match:
            continue
        slices = read_slices(match[0])
        case_data.setdefault(case, {})[os.path.basename(subdir)] = slices
        print(f"Loaded {case} from {subdir}")

# --- Plot only common cases across all subdirs
common_cases = [
    case for case, subdir_data in case_data.items()
    if all(f'sd7003_L{i}' in subdir_data for i in range(3))
]

for var, col in variables.items():
    for case in common_cases:
        subdir_data = case_data[case]
        plt.figure(figsize=(8, 6))
        for subdir, slices in subdir_data.items():
            if not slices:
                continue
            z_iter = iter(slices)
            try:
                z0 = next(z_iter)
            except StopIteration:
                continue
            data = slices[z0]
            upper = data[(data[:,4] >= 0) & (data[:,3] >= 0) & (data[:,3] <= 1)]
            lower = data[(data[:,4] < 0) & (data[:,3] >= 0) & (data[:,3] <= 1)]
            combined = np.vstack([upper, lower])
            marker = next(marker_cycle)
            label_raw = os.path.basename(subdir).split('_')[-1]  # L0, L1, L2
            label = fr"$L_{{{label_raw[-1]}}}$"
            plt.plot(combined[:,3], combined[:,col], marker=marker, linestyle='None', markersize=2, label=label)

            plot_data = {
                'case': case,
                'variable': var,
                'x': combined[:,3],
                'y': combined[:,col],
                'label': label,
                'title': format_case_title(case, var)
            }
            pkl_filename = os.path.join(output_dir, f"{case}_{var}_{label_raw}.pkl")
            with open(pkl_filename, 'wb') as pkl_file:
                pickle.dump(plot_data, pkl_file)

        plt.xlabel(r"$x/c$", fontsize=14)
        plt.ylabel(r"$C_p$" if var == 'CoefPressure' else r"$C_f$", fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        if var == 'CoefPressure':
            plt.gca().invert_yaxis()

        plt.title(format_case_title(case, var), fontsize=15)
        plt.legend(title=r"Grid Level", fontsize=12, title_fontsize=13, loc='best', frameon=True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{case}_{var}.png'), dpi=300)