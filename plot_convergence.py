import pickle
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 16,
})
import os
import numpy as np

plt.rcParams['axes.formatter.limits'] = [-15, 15]  # Force scientific notation with high precision

# Automatically find all .pkl files in the folder
script_dir = os.path.dirname(os.path.abspath(__file__))
conv_hist_dir = os.path.join(script_dir, 'conv_hist')
files = sorted([f for f in os.listdir(conv_hist_dir) if f.endswith('.pkl')])
if not files:
    raise FileNotFoundError("No convergence history .pkl file found in script directory.")

for filename in files:
    filepath = os.path.join(conv_hist_dir, filename)
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
        print(f"\nLoaded {filename}")
        for key, value in data.items():
            print(f"{key} -> {type(value)}")
            if key == 'map' and isinstance(value, dict):
                for subkey, subvalue in value.items():
                    print(f"  map[{subkey}] -> {type(subvalue)}")

    print("\nSummary of data contents:")
    for key, value in data.items():
        if isinstance(value, (list, np.ndarray)):
            shape = np.shape(value)
            print(f"{key}: type={type(value).__name__}, shape={shape}")

    var_label_map = {
        'RSDMassRMS': r'$\rho$',
        'RSDTurbulentEnergyKineticRMS': r'$k$',
        'RSDTurbulentDissRateRMS': r'$\omega$',
        'RSDTransitionGammaRMS': r' $\gamma$',
        'RSDTransitionReThetat': r'$Re_{{\theta}_t}$',
    }

    # Include all variables for normalization and plotting
    vars_to_normalize = list(var_label_map.keys())
    vars_to_plot = list(var_label_map.keys())

    for var in vars_to_normalize:
        if var in data:
            arr = np.array(data[var])
            mean_first5 = np.mean(arr[:5])
            if mean_first5 != 0:
                data[var] = arr / mean_first5
                print(f"Normalized {var} by mean of first 5 values: {mean_first5}")
            else:
                print(f"Skipped normalization of {var} (mean of first 5 values is zero)")


    # Extract a prefix name from the filename (remove directory and extension, replace spaces with underscores)
    prefix_name = os.path.splitext(filename)[0].replace(' ', '_')
    plot_dir = os.path.join(script_dir, 'convergence_history', prefix_name)
    os.makedirs(plot_dir, exist_ok=True)

    # Plot all normalized variables in one figure
    if 'total minor iters' in data:
        steps = data['total minor iters']
        plt.figure(figsize=(10, 6))
        for var in vars_to_plot:
            if var in data and len(data[var]) == len(steps):
                plt.plot(steps, data[var], label=var_label_map[var])
        plt.xlabel('Iterations')
        plt.ylabel('Normalized Value')
        #plt.title(f'{prefix_name}: Normalized Variables vs Iterations')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        #plt.savefig(os.path.join(plot_dir, f'{prefix_name}_normalized_variables.png'), dpi=600)
        plt.close()

        # Plot residuals only in a standalone figure (fancy version for thesis)
        plt.figure(figsize=(10, 6))
        for var in vars_to_normalize:
            if var in data and len(data[var]) == len(steps):
                plt.plot(steps, data[var], label=var_label_map[var])
        plt.yscale('log')
        plt.ylim(1e-15, 1e0)
        plt.xlabel('Iterations', fontsize=12)
        plt.ylabel('Normalized Residual', fontsize=12)
        #plt.title(f'{prefix_name}: Residuals Convergence', fontsize=14)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.legend(title='Residual Type', fontsize=10, title_fontsize=11, loc='upper right', frameon=True, fancybox=True, shadow=True)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f'{prefix_name}_residuals_only.png'), dpi=600)
        plt.close()


        # Plot zoomed-in normalized variables together with CL, CD, CMz as a 3x2 subplot
        fig, axes = plt.subplots(3, 2, figsize=(14, 12))
        #fig.suptitle(f'{prefix_name}: Normalized Variables and Coefficients', fontsize=14)
        axes = axes.flatten()

        # Plot normalized variables in log scale
        axes[0].set_title('Normalized Residuals')
        for var in vars_to_plot:
            if var in data and len(data[var]) == len(steps):
                axes[0].plot(steps, data[var], label=var_label_map[var])
        axes[0].set_yscale('log')
        axes[0].set_ylim(1e-15, 1e0)
        axes[0].set_xlabel('Iterations')
        axes[0].set_ylabel('Normalized Value')
        axes[0].grid(True, which='both', linestyle='--', linewidth=0.5)
        axes[0].legend(fontsize=8)

        # Plot CL (CoefLift)
        if 'CoefLift' in data and len(data['CoefLift']) == len(steps):
            axes[1].plot(steps, data['CoefLift'], label=r'$C_l$')
            axes[1].set_title(r'$C_l$ vs Iterations')
            axes[1].set_xlabel('Iterations')
            axes[1].set_ylabel(r'$C_l$')
            axes[1].set_ylim(0, 1)
            axes[1].grid(True)
            axes[1].legend()
            final_cl = data['CoefLift'][-1]
            axes[1].text(steps[-1], final_cl, f'Final: {final_cl:.4f}', fontsize=8, verticalalignment='bottom', horizontalalignment='right')
            

        # Plot CD (CoefDrag)
        if 'CoefDrag' in data and len(data['CoefDrag']) == len(steps):
            axes[2].plot(steps, data['CoefDrag'], label=r'$C_d$')
            axes[2].set_title(r'$C_d$ vs Iterations')
            axes[2].set_xlabel('Iterations')
            axes[2].set_ylabel(r'$C_d$')
            axes[2].set_ylim(-1, 1)
            axes[2].grid(True)
            axes[2].legend()
            final_cd = data['CoefDrag'][-1]
            axes[2].text(steps[-1], final_cd, f'Final: {final_cd:.4f}', fontsize=8, verticalalignment='bottom', horizontalalignment='right')
            

        # Plot CMz (CoefMomentZ)
        if 'CoefMomentZ' in data and len(data['CoefMomentZ']) == len(steps):
            axes[3].plot(steps, data['CoefMomentZ'], label=r'$C_{Mz}$')
            axes[3].set_title(r'$C_{Mz}$ vs Iterations')
            axes[3].set_xlabel('Iterations')
            axes[3].set_ylabel(r'$C_{Mz}$')
            axes[3].set_ylim(-1, 1)
            axes[3].grid(True)
            axes[3].legend()
            final_cmz = data['CoefMomentZ'][-1]
            axes[3].text(steps[-1], final_cmz, f'Final: {final_cmz:.4f}', fontsize=8, verticalalignment='bottom', horizontalalignment='right')

        # Plot CFL
        if 'CFL' in data and len(data['CFL']) == len(steps):
            axes[4].plot(steps, data['CFL'], label='CFL')
            axes[4].set_title('CFL vs Iterations')
            axes[4].set_xlabel('Iterations')
            axes[4].set_ylabel('CFL')
            axes[4].grid(True)
            axes[4].legend()
        
        # Plot iter type
        if 'iter type' in data and len(data['iter type']) == len(steps):
            axes[5].plot(steps, data['iter type'], label='Iteration type')
            axes[5].set_title('Iteration type vs Iterations')
            axes[5].set_xlabel('Iterations')
            axes[5].set_ylabel('Iteration type')
            axes[5].grid(True)
            axes[5].legend()

        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f'{prefix_name}_normalized_and_coefficients_subplots.png'), dpi=600)
        plt.close()

    if 'total minor iters' in data:
        steps = data['total minor iters']
        for key, value in data.items():
            if key not in ['map', 'total minor iters'] and isinstance(value, (list, np.ndarray)) and len(value) == len(steps):
                plt.figure(figsize=(8, 5))
                plt.plot(steps, value, label=key)
                plt.xlabel('Iterations')
                plt.ylabel(key)
                #plt.title(f'{prefix_name}: {key} vs Iterations')
                plt.grid(True)
                plt.legend()
                plt.tight_layout()
                #plt.savefig(os.path.join(plot_dir, f'{prefix_name}_{key}.png'), dpi=600)
                plt.close()
    else:
        print("Warning: 'total minor iters' not found in data. Cannot plot against Iterations.")

    map_data = data.get('map', None)
    if isinstance(map_data, dict):
        for key, value in map_data.items():
            if isinstance(value, (list, np.ndarray)) and len(value) == len(steps):
                plt.figure(figsize=(8, 5))
                plt.plot(steps, value, label=f'map[{key}]')
                plt.xlabel('Iterations')
                plt.ylabel(key)
                #plt.title(f'{prefix_name}: map[{key}] vs Iterations')
                plt.grid(True)
                plt.legend()
                plt.tight_layout()
                #plt.savefig(os.path.join(plot_dir, f'{prefix_name}_map_{key}.png'), dpi=600)
                plt.close()