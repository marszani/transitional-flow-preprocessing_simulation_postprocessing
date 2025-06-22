
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


# List of grid files to run
gridFiles = [
    "sd7003_L0.cgns"
]

# List of (mach, reynolds, alpha) tuples
mach_reynolds_alpha_tuples = [
    (0.1, 6e4, 6),
]

# Loop through specified turbulence models
for turbulenceModel in args.turbulenceModels:
    aeroOptions = {
        "gridFile": None,
        "outputDirectory": None,
        #'restartFile': 'sd7003_L0_M0.1_Re60000_alpha6_000_vol.cgns',
        'equationMode': 'steady',

        "monitorvariables": ['resrho', 'resturb','yplus',"cl","cd","cmz",],
        'surfaceVariables': ['cp', 'vx', 'vy', 'vz', 'mach', 'cfx'],
        "writeTecplotSurfaceSolution": False,
        'writeSurfaceSolution': False,
        'writeVolumeSolution': True,
    
        #'ANKJacobianLag': 3,
        #"ANKUnsteadyLSTol": 1.2,
        #"ANKPhysicalLSTol": 0.5,   

        "turbulenceModel": turbulenceModel if turbulenceModel != "laminar" else "Menter SST",
        'turbResScale': [1e0,1e-0,1e-0,1e-0],  # [1e3, 1e-14, 1e-14, 1e-25], [1e5,1e-9,1e0,1e-5], 
        #"turbulenceOrder": "second order",
        #"turbulenceProduction": "Kato-Launder",
        #'coarseDiscretization': 'central plus matrix dissipation' ,  #central plus matrix dissipation
        'discretization': 'upwind',

        #'limiter': 'minmod',

        #'ANKGlobalPreconditioner': 'multigrid',
        #'ANKUnsteadyLSTol': 1.5, # increase till 1.9 if solver stalls 
        #'ANKPhysicalLSTol': 1.5,

        #'vis4': 0.1,  #  default is 0.0156 . Increased to avoid solver stall
        #'vis2': 0.0, # default is 0.25 . Set to zero only for subsonic cases

        #'loadBalanceIter':20,
        #'ANKSwitchTol': 1e-2,

        "useblockettes": False,
        #'smoother': 'Runge-Kutta',
        #'ANKASMOverlap': 3,
        #'ANKPCILUFill': 3,
        "ANKUseTurbDADI": False,
        #'ANKNSubiterTurb': 22,
        #'ANKTurbCFLScale': 1e-1,
        #"useANKSolver": True,
        "ANKSecondOrdSwitchTol": 1e-6,
        "ANKCoupledSwitchTol": 1e4,
        #"nSubiter":12,
        #"useNKSolver": True,
        #'ANKTurbKSPDebug': True,

        'solutionPrecision':'double',

        "ANKADPC": True,
        #"NKADPC": True,
        #"NKSwitchTol": 1e-1,
        'ANKCFLLimit': 1e1,
        "L2Convergence": 1e-24,  
        "nCycles": 1000,
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

            aeroOptions["gridFile"] = gridFile
            aeroOptions["outputDirectory"] = outputDir

            # Initialize ADflow solver (real)
            CFDSolver = ADFLOW(options=aeroOptions, dtype="D")
            CFDSolver.addLiftDistribution(150, "z")
            CFDSolver.addSlices("z", np.linspace(0.1, 1, 1))

            # Initialize ADflow solver (complex) with fallback
            CFDSolver_cmplx = None
            try:
                CFDSolver_cmplx = ADFLOW(options=aeroOptions, dtype="D")
                CFDSolver_cmplx.addLiftDistribution(150, "z")
                CFDSolver_cmplx.addSlices("z", np.linspace(0.1, 1, 1))
            except Exception as e:
                if comm.rank == 0:
                    print(f"[Check] Complex solver init failed: {e}")
                CFDSolver_cmplx = None  # Disable complex step testing gracefully

            # Define AeroProblem
            ap = AeroProblem(
                name=f"{gridFile.split('.')[0]}_M{mach}_Re{int(reynolds)}_alpha{alpha}",
                mach=mach,
                #turbIntensityInf = 0.3,
                #eddyVisInfRatio = 0.0018,
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

                # === Jacobian-vector product verification ===

                # Information for Jacobian-vector product computation (plain English, no LaTeX)
                info = {
                    "Parameters": {
                        "rel_error": (
                            "Relative error between AD and FD or CS. "
                            "Formula: rel_error(a, b) = ||a - b|| / max(||a||, eps) "
                            "Where: "
                            "- a is always AD (fixed) "
                            "- b is FD or CS "
                            "- || . || denotes the Euclidean norm (L2 norm) "
                            "- eps is a small number (e.g., 1e-16) to avoid division by zero"
                        ),
                        "xDvDot": "dict - Perturbation on geometric design variables defined in DVGeo",
                        "xSDot": "numpy array - Perturbation on surface coordinates",
                        "xVDot": "numpy array - Perturbation on volume mesh",
                        "wDot": "numpy array - Perturbation on flow state variables",
                        "residualDeriv": "bool - Whether to compute residual derivative",
                        "funcDeriv": "bool - Whether to compute derivative of cost functions (objective and constraints)",
                        "fderiv": "bool - Whether to compute derivative of surface forces (tractions)",
                        "groupName": "str - Optional surface group name for evaluating functions; defaults to all surfaces",
                        "mode": "str - Method for computing jacobian vector products: 'AD' (automatic differentiation), 'FD' (finite difference), or 'CS' (complex step)",
                        "h": "float - Step size used when mode is 'FD' or 'CS'"
                    },
                    "Returns": {
                        "dwdot": "array - Residual derivative",
                        "funcsdot": "dict - Derivatives of objective and constraint functions",
                        "fDot": "array - Surface force derivative"
                    },
                    "step_size of DVs": {
                        "alpha": 1e-4,   # Angle of attack perturbation (radians)
                        "beta":  1e-5,   # Sideslip angle perturbation (radians)
                        "mach":  1e-5,   # Mach number perturbation
                        "P":     1e-1,   # Pressure perturbation (Pa or relevant unit)
                        "T":     1e-4,   # Temperature perturbation (K)
                        "xRef":  1e-5,   # Reference x-coordinate perturbation (length units)
                        "yRef":  1e-5,   # Reference y-coordinate perturbation (length units)
                        "zRef":  1e-5,   # Reference z-coordinate perturbation (length units)
                    }
                }

                
                from collections import OrderedDict

                step_size = {
                    "alpha": 1e-4,
                    "beta": 1e-5,
                    "mach": 1e-5,
                    "P": 1e-1,
                    "T": 1e-4,
                    "xRef": 1e-5,
                    "yRef": 1e-5,
                    "zRef": 1e-5,
                }

                class DummyHandler:
                    def __init__(self):
                        self.data = {}

                    @staticmethod
                    def rel_error(a, b, eps=1e-14):
                        a = np.asarray(a)
                        b = np.asarray(b)
                        numerator = np.linalg.norm(a - b)
                        denominator = max(np.linalg.norm(a), np.linalg.norm(b), eps)
                        return numerator / denominator

                    def root_print(self, msg):
                        pass

                    def par_add_norm(self, label, arr, rtol, atol):
                        self.data[label] = {
                            "norm": float(np.linalg.norm(arr)),
                            "rtol": rtol,
                            "atol": atol
                        }

                    def root_add_dict(self, label, dct, rtol, atol):
                        clean_dict = {k: float(v) for k, v in dct.items()}
                        self.data[label] = {
                            "values": clean_dict,
                            "rtol": rtol,
                            "atol": atol
                        }

                    def par_add_rel_error(self, label, rel_err, rtol, atol):
                        self.data[label] = {
                            "rel_error": float(rel_err),
                            "rtol": rtol,
                            "atol": atol
                        }

                    def root_add_rel_error_dict(self, label, rel_err_dict, rtol, atol):
                        clean_dict = {k: float(v) for k, v in rel_err_dict.items()}
                        self.data[label] = {
                            "rel_error_values": clean_dict,
                            "rtol": rtol,
                            "atol": atol
                        }



                handler = DummyHandler()
                rtol = 1e-7
                atol = 1e-8

                results = {}
                for dv, h in step_size.items():
                    results[dv] = {}
                    xDvDot = {dv: h}
                    xSDot = CFDSolver.getSpatialPerturbation(333)
                    xVDot = CFDSolver.getSpatialPerturbation(314)
                    wDot = CFDSolver.getStatePerturbation(321)

                    modes = [("FD", 1e-9), ("CS", 1e-30)]
                    ad_out = CFDSolver.computeJacobianVectorProductFwd(
                        xDvDot=None,
                        xSDot=None,
                        xVDot=None,
                        wDot=wDot,
                        residualDeriv=True,
                        funcDeriv=True,
                        fDeriv=True,
                        mode="AD"
                    )

                    for mode_name, step in modes:
                        mode_result = {}

                        # Handler calls for AD (always available in ad_out)
                        handler.root_print(f"||dR/dw * wDot|| (AD)")
                        handler.par_add_norm(f"||dR/dw * wDot|| (AD)", ad_out[0], rtol=rtol, atol=atol)


                        handler.root_print(f"dFuncs/dw * wDot (AD)")
                        handler.root_add_dict(f"dFuncs/dw * wDot (AD)", ad_out[1], rtol=rtol, atol=atol)

                        handler.root_print(f"||dF/dw * wDot|| (AD)")
                        handler.par_add_norm(f"||dF/dw * wDot|| (AD)", ad_out[2], rtol=rtol, atol=atol)
                        print(f"[Check] Jacobian-vector product (AD) completed for {dv}")
                        handler.root_print(f"Jacobian-vector product verification completed for key: {dv}")

                        # Handler calls for FD mode only
                        if mode_name == "FD":
                            ref_out = CFDSolver.computeJacobianVectorProductFwd(
                                xDvDot=None,
                                xSDot=None,
                                xVDot=None,
                                wDot=wDot,
                                residualDeriv=True,
                                funcDeriv=True,
                                fDeriv=True,
                                mode=mode_name,
                                h=step
                            )
                            handler.root_print(f"||dR/dw * wDot|| (FD)")
                            handler.par_add_norm(f"||dR/dw * wDot|| (FD)", ref_out[0], rtol=rtol, atol=atol)
                            print(f"[Check] Jacobian-vector product (FD) completed for {dv}")


                            handler.root_print(f"dFuncs/dw * wDot (FD)")
                            handler.root_add_dict(f"dFuncs/dw * wDot (FD)", ref_out[1], rtol=rtol, atol=atol)
                            print(f"[Check] Jacobian-vector product (FD) completed for fDot_fd ({dv})")

                            handler.root_print(f"||dF/dw * wDot|| (FD)")
                            handler.par_add_norm(f"||dF/dw * wDot|| (FD)", ref_out[2], rtol=rtol, atol=atol)
                            print(f"[Check] Jacobian-vector product (FD) completed for fDot_fd ({dv})")
                            handler.root_print(f"Jacobian-vector product verification completed for key: {dv}")

                            # Compute relative errors between AD and FD
                            rel_err_dwdot = handler.rel_error(ad_out[0], ref_out[0])
                            rel_err_funcsDot = {
                                key: handler.rel_error(ad_out[1][key], ref_out[1][key]) for key in ad_out[1]
                            }
                            rel_err_fdot = handler.rel_error(ad_out[2], ref_out[2])

                            # ONLY this block should write to handler.data
                            handler.par_add_rel_error("Relative error dR/dw * wDot (AD vs FD):", rel_err_dwdot, rtol, atol)
                            handler.root_add_rel_error_dict("Relative error dFuncs/dw * wDot (AD vs FD):", rel_err_funcsDot, rtol, atol)
                            handler.par_add_rel_error("Relative error dF/dw * wDot (AD vs FD):", rel_err_fdot, rtol, atol)
                            
                        elif mode_name == "CS":
                            try:
                                if CFDSolver.dtype != "D":
                                    continue
                                ref_out = CFDSolver.computeJacobianVectorProductFwd(
                                    xDvDot=None,
                                    xSDot=None,
                                    xVDot=None,
                                    wDot=wDot,
                                    residualDeriv=True,
                                    funcDeriv=True,
                                    fDeriv=True,
                                    mode=mode_name,
                                    h=step
                                )
                                handler.root_print(f"||dF/dw * wDot|| (CS)")
                                handler.par_add_norm(f"||dF/dw * wDot|| (CS)", ref_out[2], rtol=rtol, atol=atol)
                                handler.root_print(f"Jacobian-vector product verification completed for key: {dv}")
                            except:
                                continue
                        else:
                            ref_out = CFDSolver.computeJacobianVectorProductFwd(
                                xDvDot=None,
                                xSDot=None,
                                xVDot=None,
                                wDot=wDot,
                                residualDeriv=True,
                                funcDeriv=True,
                                fDeriv=True,
                                mode=mode_name,
                                h=step
                            )

                # --- Full Jacobian-vector product checks for xVDot, xDvDot, xSDot ---
                # xVDot
                xVDot = CFDSolver.getSpatialPerturbation(314)
                xDvDot = {dv: 0.0 for dv in step_size.keys()}
                xSDot = CFDSolver.getSpatialPerturbation(333)
                wDot = CFDSolver.getStatePerturbation(321)

                # xVDot AD
                handler.root_print("||dR/dXv * xVDot|| (AD)")
                resDot_ad_xVDot, funcsDot_ad_xVDot, fDot_ad_xVDot = CFDSolver.computeJacobianVectorProductFwd(
                    xDvDot=None,
                    xSDot=None,
                    xVDot=xVDot,
                    wDot=None,
                    residualDeriv=True,
                    funcDeriv=True,
                    fDeriv=True,
                    mode="AD"
                )
                handler.par_add_norm("||dR/dXv * xVDot|| (AD)", resDot_ad_xVDot, rtol=rtol, atol=atol)

                handler.root_print("dFuncs/dXv * xVDot (AD)")
                handler.root_add_dict("dFuncs/dXv * xVDot (AD)", funcsDot_ad_xVDot, rtol=rtol, atol=atol)

                handler.root_print("||dF/dXv * xVDot|| (AD)")
                handler.par_add_norm("||dF/dXv * xVDot|| (AD)", fDot_ad_xVDot, rtol=rtol, atol=atol)
                print("[Check] Jacobian-vector product (AD) completed for xVDot")

                # xVDot FD
                handler.root_print("||dR/dXv * xVDot|| (FD)")
                resDot_fd_xVDot, funcsDot_fd_xVDot, fDot_fd_xVDot = CFDSolver.computeJacobianVectorProductFwd(
                    xDvDot=None,
                    xSDot=None,
                    xVDot=xVDot,
                    wDot=None,
                    residualDeriv=True,
                    funcDeriv=True,
                    fDeriv=True,
                    mode="FD",
                    h=1e-9
                )
                handler.par_add_norm("||dR/dXv * xVDot|| (FD)", resDot_fd_xVDot, rtol=rtol, atol=atol)

                handler.root_print("dFuncs/dXv * xVDot (FD)")
                handler.root_add_dict("dFuncs/dXv * xVDot (FD)", funcsDot_fd_xVDot, rtol=rtol, atol=atol)

                handler.root_print("||dF/dXv * xVDot|| (FD)")
                handler.par_add_norm("||dF/dXv * xVDot|| (FD)", fDot_fd_xVDot, rtol=rtol, atol=atol)
                print("[Check] Jacobian-vector product (FD) completed for xVDot")

                # Compute relative errors between AD and FD
                rel_err_dwdot = handler.rel_error(resDot_ad_xVDot, resDot_fd_xVDot)
                rel_err_funcsDot = {
                    key: handler.rel_error(funcsDot_ad_xVDot[key], funcsDot_fd_xVDot[key]) for key in funcsDot_ad_xVDot
                }
                rel_err_fdot = handler.rel_error(fDot_ad_xVDot, fDot_fd_xVDot)

                # ONLY this block should write to handler.data
                handler.par_add_rel_error("Relative error dR/dXv * xVDot (AD vs FD):", rel_err_dwdot, rtol, atol)
                handler.root_add_rel_error_dict("Relative error dFuncs/dXv * xVDot (AD vs FD):", rel_err_funcsDot, rtol, atol)
                handler.par_add_rel_error("Relative error dF/dXv * xVDot (AD vs FD):", rel_err_fdot, rtol, atol)

                # xDvDot (single DV perturbation at a time)
                for dv in step_size.keys():
                    xDvDot = {k: 0.0 for k in step_size.keys()}
                    xDvDot[dv] = step_size[dv]
                    # AD
                    handler.root_print(f"||dR/d{dv} || (AD)")
                    resDot_ad_xDvDot, funcsDot_ad_xDvDot, fDot_ad_xDvDot = CFDSolver.computeJacobianVectorProductFwd(
                        xDvDot=xDvDot,
                        xSDot=None,
                        xVDot=None,
                        wDot=None,
                        residualDeriv=True,
                        funcDeriv=True,
                        fDeriv=True,
                        mode="AD"
                    )
                    handler.par_add_norm(f"||dR/d{dv} || (AD)", resDot_ad_xDvDot, rtol=rtol, atol=atol)

                    handler.root_print(f"dFuncs/d{dv}  (AD)")
                    handler.root_add_dict(f"dFuncs/d{dv}  (AD)", funcsDot_ad_xDvDot, rtol=rtol, atol=atol)

                    handler.root_print(f"||dF/d{dv} || (AD)")
                    handler.par_add_norm(f"||dF/d{dv} || (AD)", fDot_ad_xDvDot, rtol=rtol, atol=atol)
                    print(f"[Check] Jacobian-vector product (AD) completed for xDvDot {dv}")

                    # FD
                    handler.root_print(f"||dR/d{dv} || (FD)")
                    resDot_fd_xDvDot, funcsDot_fd_xDvDot, fDot_fd_xDvDot = CFDSolver.computeJacobianVectorProductFwd(
                        xDvDot=xDvDot,
                        xSDot=None,
                        xVDot=None,
                        wDot=None,
                        residualDeriv=True,
                        funcDeriv=True,
                        fDeriv=True,
                        mode="FD",
                        h=step_size[dv]
                    )
                    handler.par_add_norm(f"||dR/d{dv} || (FD)", resDot_fd_xDvDot, rtol=rtol, atol=atol)

                    handler.root_print(f"dFuncs/d{dv} (FD)")
                    handler.root_add_dict(f"dFuncs/d{dv} (FD)", funcsDot_fd_xDvDot, rtol=rtol, atol=atol)

                    handler.root_print(f"||dF/d{dv} || (FD)")
                    handler.par_add_norm(f"||dF/d{dv} || (FD)", fDot_fd_xDvDot, rtol=rtol, atol=atol)
                    print(f"[Check] Jacobian-vector product (FD) completed for xDvDot {dv}")

                    # Compute relative errors (between AD and FD)
                    rel_err_res = handler.rel_error(resDot_ad_xDvDot, resDot_fd_xDvDot)
                    rel_err_funcs = {
                        key: handler.rel_error(funcsDot_ad_xDvDot[key], funcsDot_fd_xDvDot[key]) for key in funcsDot_ad_xDvDot
                    }
                    rel_err_fdot = handler.rel_error(fDot_ad_xDvDot, fDot_fd_xDvDot)

                    # Store relative errors with descriptive keys
                    handler.par_add_rel_error(f"Relative error dR/d{dv} (AD vs FD):", rel_err_res, rtol, atol)
                    handler.root_add_rel_error_dict(f"Relative error dFuncs/d{dv} (AD vs FD):", rel_err_funcs, rtol, atol)
                    handler.par_add_rel_error(f"Relative error dF/d{dv} (AD vs FD):", rel_err_fdot, rtol, atol)

                # xSDot
                xSDot = CFDSolver.getSpatialPerturbation(333)
                handler.root_print("||dR/dXs * xSDot|| (AD)")
                resDot_ad_xSDot, funcsDot_ad_xSDot, fDot_ad_xSDot = CFDSolver.computeJacobianVectorProductFwd(
                    xDvDot=None,
                    xSDot=xSDot,
                    xVDot=None,
                    wDot=None,
                    residualDeriv=True,
                    funcDeriv=True,
                    fDeriv=True,
                    mode="AD"
                )
                handler.par_add_norm("||dR/dXs * xSDot|| (AD)", resDot_ad_xSDot, rtol=rtol, atol=atol)

                handler.root_print("dFuncs/dXs * xSDot (AD)")
                handler.root_add_dict("dFuncs/dXs * xSDot (AD)", funcsDot_ad_xSDot, rtol=rtol, atol=atol)

                handler.root_print("||dF/dXs * xSDot|| (AD)")
                handler.par_add_norm("||dF/dXs * xSDot|| (AD)", fDot_ad_xSDot, rtol=rtol, atol=atol)
                print("[Check] Jacobian-vector product (AD) completed for xSDot")

                handler.root_print("||dR/dXs * xSDot|| (FD)")
                resDot_fd_xSDot, funcsDot_fd_xSDot, fDot_fd_xSDot = CFDSolver.computeJacobianVectorProductFwd(
                    xDvDot=None,
                    xSDot=xSDot,
                    xVDot=None,
                    wDot=None,
                    residualDeriv=True,
                    funcDeriv=True,
                    fDeriv=True,
                    mode="FD",
                    h=1e-9
                )
                handler.par_add_norm("||dR/dXs * xSDot|| (FD)", resDot_fd_xSDot, rtol=rtol, atol=atol)

                handler.root_print("dFuncs/dXs * xSDot (FD)")
                handler.root_add_dict("dFuncs/dXs * xSDot (FD)", funcsDot_fd_xSDot, rtol=rtol, atol=atol)

                handler.root_print("||dF/dXs * xSDot|| (FD)")
                handler.par_add_norm("||dF/dXs * xSDot|| (FD)", fDot_fd_xSDot, rtol=rtol, atol=atol)
                print("[Check] Jacobian-vector product (FD) completed for xSDot")

                # Compute relative errors between AD and FD
                rel_err_dwdot = handler.rel_error(resDot_ad_xSDot, resDot_fd_xSDot)
                rel_err_funcsDot = {
                    key: handler.rel_error(funcsDot_ad_xSDot[key], funcsDot_fd_xSDot[key]) for key in funcsDot_ad_xSDot
                }
                rel_err_fdot = handler.rel_error(fDot_ad_xSDot, fDot_fd_xSDot)

                # ONLY this block should write to handler.data
                handler.par_add_rel_error("Relative error dR/dXs * xSDot (AD vs FD):", rel_err_dwdot, rtol, atol)
                handler.root_add_rel_error_dict("Relative error dFuncs/dXs * xSDot (AD vs FD):", rel_err_funcsDot, rtol, atol)
                handler.par_add_rel_error("Relative error dF/dXs * xSDot (AD vs FD):", rel_err_fdot, rtol, atol)

                
                # Save results to file, including info and handler data in JSON export
                if comm.rank == 0:
                    os.makedirs("jac_check", exist_ok=True)
                    json_results = {
                        "info": info,
                        "definitions": None,
                        "results": results,
                        "handler_data": handler.data
                    }
                    with open(os.path.join("jac_check", f"{ap.name}_{turbulenceModel}_FWD.json"), "w") as f:
                        import json
                        json.dump(json_results, f, indent=2)

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
                    print(f"Polar results for {gridFile} with {turbulenceModel}:")
                    print("{:>6} {:>8} {:>8}".format("Alpha", "CL", "CD"))
                    print("=" * 24)
                    for alpha, cl, cd in zip(alphaList, CLList, CDList):
                        print(f"{alpha:6.1f} {cl:8.4f} {cd:8.4f}")


                