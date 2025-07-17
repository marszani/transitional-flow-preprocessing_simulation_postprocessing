# runs of cases laminar, SST, LM-SST and saves convergence history, uses also a restart file if present 

import numpy as np
import argparse
import os
import stat
from adflow import ADFLOW
from baseclasses import AeroProblem
from mpi4py import MPI
import shutil
import pickle

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, default="sim")
parser.add_argument("--task", choices=["analysis", "polar"], default="analysis")
parser.add_argument("--turbulenceModels", type=str, nargs="*", default=["laminar"], 
                    help="Specify turbulence models. Defaults to ['laminar'].")
args = parser.parse_args()

# Comment the next three lines for only laminar simulation

# Fallback if turbulenceModels is unchanged from default
if args.turbulenceModels == ["laminar"]:  
    args.turbulenceModels = ["Menter SST" ,"Langtry Menter SST"]
    args.turbulenceModels = ["Langtry Menter SST"]

# MPI setup
comm = MPI.COMM_WORLD

# Ensure main output directory exists before creating subdirectories
if comm.rank == 0:
    os.makedirs(args.output, exist_ok=True)

comm.barrier()  # Synchronize MPI ranks before proceeding

# Function to handle read-only files when removing directories
def handle_remove_readonly(func, path, exc):
    os.chmod(path, stat.S_IWRITE)
    func(path)

def save_conv_history(Solver, AP, model):
    hist = Solver.getConvergenceHistory()
    if MPI.COMM_WORLD.rank == 0:
        os.makedirs('conv_hist', exist_ok=True)
        filepath = os.path.join('conv_hist', f'conv_hist_{AP.name}_{model}.pkl')
        with open(filepath, "wb") as f:
            pickle.dump(hist, f)

def check_restart_file(restartFile):
    if os.path.isfile(restartFile):
        print(f"[Info] Restart file {restartFile} found.")
        return True
    else:
        print(f"[Warning] Restart file {restartFile} not found.")
        return False

# List of grid files to run
gridFiles = [
    "e387_L0.cgns"
]

# List of (mach, reynolds, alpha) tuples
mach_reynolds_alpha_tuples = [
    (0.1, 2e5, 0),
    (0.1, 2e5, 4),
]

# Loop through specified turbulence models
for turbulenceModel in args.turbulenceModels:
    aeroOptions = {
        "gridFile": None,
        "outputDirectory": None,
        'equationMode': 'steady',

        "monitorvariables": ['resrho', 'resturb','yplus',"cl","cd","cmz",],
        'surfaceVariables': ['cp', 'vx', 'vy', 'vz', 'mach', 'cfx'],
        "writeTecplotSurfaceSolution": False,
        'writeSurfaceSolution': False,
        'writeVolumeSolution': True,


        #'ANKJacobianLag': 3,
        #"ANKUnsteadyLSTol": 1.5,
        #"ANKPhysicalLSTol": 0.5,   
        
        #'MGCycle': '3w',

        
        "turbulenceModel": turbulenceModel if turbulenceModel != "laminar" else "Menter SST",
        'turbResScale': [1e0,1e0,1e0,1e0], # [1e3, 1e-14, 1e-14, 1e-25],
        #"turbulenceOrder": "second order",
        #"turbulenceProduction": "Kato-Launder",
        #'coarseDiscretization': 'central plus matrix dissipation' ,  
        'discretization': 'upwind',
        
        #'limiter': 'minmod',

        #'ANKGlobalPreconditioner': 'multigrid',
        #'ANKUnsteadyLSTol': 1.5, # increase till 1.9 if solver stalls 
        #'ANKPhysicalLSTol': 1.5,

        

        #'vis4': 0.1,  #  default is 0.0156 . Increased to avoid solver stall
        #'vis2': 0.25, # default is 0.25 . Set to zero only for subsonic cases

        #'loadBalanceIter':20,
        #'ANKSwitchTol': 1e-2,

        "useblockettes": False,
        #'smoother': 'Runge-Kutta',
        #'ANKASMOverlap': 3,
        #'ANKPCILUFill': 3,
        "ANKUseTurbDADI": False,
        #'ANKNSubiterTurb': 22,
        #'ANKTurbCFLScale': 1e-1,
        "useANKSolver": True,
        "ANKSecondOrdSwitchTol": 1e-6,
        "ANKCoupledSwitchTol": 1e4,
        #"nSubiter":12,
        #"useNKSolver": True,
        #'ANKTurbKSPDebug': True,

        'solutionPrecision':'double',

        "ANKADPC": True,
        "NKADPC": True,
        #"NKSwitchTol": 1e-1,
        'ANKCFLLimit': 1e2,
        'ANKCFL0': 5e-1,
        "L2Convergence": 1e-14,  
        "nCycles": 5000,
        'turbIntensityInf' : 0.0008, # range  0 : 3
        'eddyVisInfRatio' : 0.001, # range 0.1 : 100'

    }
    if turbulenceModel == "laminar":
        aeroOptions["equationType"] = "laminar NS"
    else:
        aeroOptions["equationType"] = "RANS"

    for mach, reynolds, alpha in mach_reynolds_alpha_tuples:
        for gridFile in gridFiles:
            # Create a unique output directory for each grid file, turbulence model, Mach, Reynolds, and alpha (only rank 0)
            outputDir = os.path.join(args.output, gridFile.split(".")[0], turbulenceModel, f"M{mach}_Re{int(reynolds)}_alpha{alpha}")

            if comm.rank == 0:
                # if os.path.exists(outputDir):
                #     shutil.rmtree(outputDir, onerror=handle_remove_readonly)
                os.makedirs(outputDir, exist_ok=True)
                print(f"[Setup] Overwritten output directory: {outputDir} | Model: {turbulenceModel} | Mach: {mach} | Re: {reynolds:.1e} | Alpha: {alpha}")

            # Dynamically generate the restart file path for each case
            restart_dir = os.path.join(
                args.output, gridFile.split(".")[0], turbulenceModel,
                f'M{mach}_Re{int(reynolds)}_alpha{alpha}'
            )
            restart_filename = f'{gridFile.split(".")[0]}_M{mach}_Re{int(reynolds)}_alpha{alpha}_000_vol.cgns'
            restartFile = os.path.join(restart_dir, restart_filename)
            print(f"Generated restart file path: {restartFile}")  # Add this print statement for debugging
            aeroOptions["restartFile"] = restartFile

            if check_restart_file(restartFile):

                aeroOptions["gridFile"] = gridFile
                aeroOptions["outputDirectory"] = outputDir

                # Initialize ADflow solver
                CFDSolver = ADFLOW(options=aeroOptions)
                CFDSolver.addLiftDistribution(150, "z")
                CFDSolver.addSlices("z", np.linspace(0.1, 1, 1))

                # Define AeroProblem
                ap = AeroProblem(
                    name=f"{gridFile.split('.')[0]}_M{mach}_Re{int(reynolds)}_alpha{alpha}",
                    mach=mach,
                    T=288,
                    reynolds=reynolds,
                    reynoldsLength = 1.0,
                    alpha=3,
                    areaRef=1.0,
                    chordRef=1.0,
                    evalFuncs=["cl", "cd","cmz"],
                )
                ap.alpha = alpha

                # Run ADflow in analysis mode
                if args.task == "analysis":
                    CFDSolver(ap)
                    funcs = {}
                    CFDSolver.evalFunctions(ap, funcs)
                    save_conv_history(CFDSolver, ap, turbulenceModel)

                    if comm.rank == 0:
                            print(f"[Analysis] {gridFile} | Model: {turbulenceModel} | Mach: {mach} | Re: {reynolds:.1e} => {funcs}")

                # Run ADflow in polar mode (sweeping AoA)
                elif args.task == "polar":
                    alphaList = np.linspace(0, 5, 6)
                    CLList, CDList = [], []

                    for alpha in alphaList:
                        ap.name = f"{gridFile.split('.')[0]}_M{mach}_Re{int(reynolds)}_{alpha:.2f}"
                        ap.alpha = alpha

                        if comm.rank == 0:
                            print(f"Running alpha = {ap.alpha} for {gridFile} with {turbulenceModel}")

                        CFDSolver(ap)
                        funcs = {}
                        CFDSolver.evalFunctions(ap, funcs)
                        save_conv_history(CFDSolver, ap, turbulenceModel)

                        if comm.rank == 0:
                                print(f"[Polar] {gridFile} | Model: {turbulenceModel} | Mach: {mach} | Re: {reynolds:.1e} | Alpha: {alpha:.2f} => {funcs}")

                        CLList.append(funcs.get(f"{ap.name}_cl", np.nan))
                        CDList.append(funcs.get(f"{ap.name}_cd", np.nan))

                    if comm.rank == 0:
                        print(f"\nPolar results for {gridFile} with {turbulenceModel}:")
                        print("{:>6} {:>8} {:>8}".format("Alpha", "CL", "CD"))
                        print("=" * 24)
                        for alpha, cl, cd in zip(alphaList, CLList, CDList):
                            print(f"{alpha:6.1f} {cl:8.4f} {cd:8.4f}")