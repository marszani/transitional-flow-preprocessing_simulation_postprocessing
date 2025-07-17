from pyhyp import pyHyp
from prefoil import sampling, Airfoil
import os
import numpy as np
import matplotlib.pyplot as plt

def read_dat_file(filepath):
    """Reads a .dat airfoil file and returns numpy array of coordinates."""
    coords = []
    with open(filepath, 'r') as f:
        for line in f:
            try:
                parts = line.strip().split()
                if len(parts) == 2:
                    x, y = map(float, parts)
                    coords.append([x, y])
            except ValueError:
                continue
    return np.array(coords)

# Input airfoil .dat path
airfoil_dat_path = './e387.dat'
initCoords = read_dat_file(airfoil_dat_path)
airfoil = Airfoil(initCoords)

# Output directory
output_dir = './input'
os.makedirs(output_dir, exist_ok=True)

# Base L2 mesh parameters
nTE_cells_L2 = 1
nSurfPts_L2 = 300
nLayers_L2 = 150
s0_L2 = 1e-6

# Refinement setup
refinements = [0.2, 1.0, 2.0]
levels = ["L2", "L1", "L0"]

for r, level in zip(refinements, levels):
    nSurfPts = int(r * nSurfPts_L2)
    nTEPts = int(r * nTE_cells_L2)
    nExtPts = int(r * nLayers_L2)
    s0 = s0_L2 / r

    coords = airfoil.getSampledPts(
        nSurfPts,
        spacingFunc=sampling.polynomial,
        func_args={"order": 8},
        nTEPts=nTEPts,
    )

    base_filename = f"{os.path.splitext(os.path.basename(airfoil_dat_path))[0]}_{level}"
    
    # Save as .dat file
    dat_path = os.path.join(output_dir, f"{base_filename}.dat")
    np.savetxt(dat_path, coords, fmt="%.10f")

    # Save as plot3d (.xyz)
    airfoil.writeCoords(os.path.join(output_dir, base_filename), file_format="plot3d")

    # Plot and save figure
    plt.figure()
    plt.plot(coords[:, 0], coords[:, 1])
    plt.axis("equal")
    plt.title(f"{base_filename}")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"{base_filename}.png"), dpi=300, bbox_inches='tight')
    plt.close()

    # pyHyp mesh generation
    options = {
        "inputFile": os.path.join(output_dir, f"{base_filename}.xyz"),
        "unattachedEdgesAreSymmetry": False,
        "outerFaceBC": "farfield",
        "autoConnect": True,
        "BC": {1: {"jLow": "zSymm", "jHigh": "zSymm"}},
        "families": "wall",
        "N": nExtPts,
        "s0": s0,
        "marchDist": 20.0,
    }

    try:
        hyp = pyHyp(options=options)
        hyp.run()
        hyp.writeCGNS(os.path.join(output_dir, f"{base_filename}_20_march_dist.cgns"))
    except Exception as e:
        print(f"Error processing {level}: {e}")