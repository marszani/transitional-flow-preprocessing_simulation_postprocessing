Here is your provided content rewritten clearly and precisely in Markdown format, keeping all equations and code intact:

⸻

Code Implementation

The Fortran routines implementing the Langtry–Menter \gamma–Re_{\theta} transition model have been developed and integrated into the ADflow solver framework. The transition model is compatible with all available solvers, including DADI, multigrid, and ANK. It is coupled to the flow variables and solved simultaneously via the coupled ANK formulation (CANK), and also supports the second-order accurate and coupled second-order accurate implicit formulations (SANK and CSANK) adflowSolvers.

In the Python interface pyADflow.py, key modifications enable specifying freestream turbulence intensity (TI_{\infty}) and freestream eddy viscosity ratio \left(\mu_t / \mu\right)_{\infty}, allowing users to control these parameters directly through ADflow options.

⸻

Code Documentation

Clear and comprehensive documentation is essential to ensure well-designed code remains usable, maintainable, and extensible over time. This section provides detailed documentation for implementing the Langtry–Menter transition model, excluding minor integrations in specific solver routines.

The goal is to present documentation straightforwardly and descriptively, enabling future developers to understand, modify, or build upon the codebase effectively. This effort should facilitate ongoing development and foster collaboration within the community.

⸻

GammaRethetaModel module

The GammaRethetaModel module implements the transition model based on the momentum-thickness Reynolds number, Re_{\theta t}, used in the Langtry-Menter correlation-based transition framework. It provides all routines to compute, diffuse, and scale transition variables such as intermittency \gamma and transition onset Reynolds number Re_{\theta t}.

Core Functionality
	•	solve_local_Re_thetatt_eq: Solves an implicit equation for local transition Reynolds number Re_{\theta t} via a secant method. This includes:
	•	Computing local velocity magnitude and gradients.
	•	Evaluating turbulence intensity and relevant flow parameters.
	•	Iterative nonlinear correlation solving with convergence checks and clipping for numerical stability.
	•	GammaRethetaSource: Computes source terms for transition transport equations:
	•	Loops over all interior cells.
	•	Calculates velocity gradients, vorticity, turbulence intensity, and local Re_{\theta t}.
	•	Applies empirical correlations and blending functions to evaluate transition onset and progression.
	•	GammaRethetaViscous: Computes viscous diffusion terms for transition variables in \xi, \eta, \zeta directions:
	•	Loops over interior grid points excluding ghost cells.
	•	Computes direction-specific volume metrics, e.g.:
\begin{aligned}
\text{voli} &= \frac{1}{\text{vol}(i,j,k)}, \\
\text{volmi} &= \frac{2}{\text{vol}(i,j,k) + \text{vol}(i,j,k-1)}, \quad \text{etc.}
\end{aligned}
	•	Uses central differences to approximate second derivatives:
\frac{\partial}{\partial \zeta} \left( \mu \frac{\partial \phi}{\partial \zeta} \right) \approx \frac{1}{\rho} \left[ c_{1m} \phi_{k-1} - c_{10} \phi_{k} + c_{1p} \phi_{k+1} \right]
	•	Two distinct contributions:
	•	\gamma-equation diffusion scaled by laminar and eddy viscosity with calibration constant rLMsigmaf.
	•	Re_{\theta t}-equation diffusion scaled similarly with rLMsigmathetat.
	•	Updates scratch array with directional contributions.
	•	GammaRethetaResScale: Scales residuals for solver consistency:
	•	Iterates interior points.
	•	Scales residuals by reference cell volume and applies masking:
dw(i,j,k) = - \text{volRef}(i,j,k) \cdot \text{scratch}(i,j,k) \cdot \text{iblank}(i,j,k)

⸻

Implementation Notes
	•	Secant method avoids analytical derivatives, converging fast (under 10 iterations).
	•	Routines tightly coupled to local flow conditions for accurate transitional modeling.
	•	Numerical stability via limiting, clipping, ghost-cell treatment.
	•	Module integrates seamlessly into turbulence framework to modulate source terms during transition.

⸻

SST_block_residuals and SST_block_residuals_d Subroutines with Transition Model Support (module SST)

These subroutines compute residuals for the SST turbulence model, optionally coupled with the Gamma-Rethetat transition model. The _d variant uses Automatic Differentiation (AD) for sensitivity and adjoint evaluations.

Overview
	•	SST_block_residuals: Computes block residuals.
	•	SST_block_residuals_d: Computes residuals and derivatives via AD hooks.
	•	Both evaluate transition terms if enabled, then standard SST turbulence terms.

Transition Model: Gamma-Re_{\theta t}

When transitionModel = gammaRetheta:
	•	Compute local strain rate and vorticity squared (strainNorm2, prodWmag2).
	•	Call GammaRethetaSource for transition source terms.
	•	Add advection terms via turbAdvection.
	•	Compute viscous diffusion via GammaRethetaViscous.
	•	Scale residuals via GammaRethetaResScale.

AD-aware versions are called in _d.

SST Turbulence Model Terms (always executed)
	•	kwCDterm: Cross-diffusion term in \omega-equation.
	•	f1SST: SST blending function F_1.
	•	prodSmag2, prodWmag2, prodKatoLaunder: Turbulence production.
	•	SSTSource: k and \omega source terms.
	•	turbAdvection: Advection terms.
	•	SSTViscous: Viscous diffusion terms.
	•	SSTResScale: Scales residuals by cell volume.

Automatic Differentiation (AD) Version
	•	_d mirrors primal routine with differentiated subroutines.
	•	Needed for gradient optimization and adjoint methods.
	•	Structure consistency ensures AD correctness.

Cleanup and Memory Management
	•	Temporary Jacobian/flux array qq allocated/deallocated based on cleanUp.

⸻

SSTSource Subroutine: SST Source Term Evaluation with Transition Model Support (module SST)

Computes source terms and partial Jacobians for k-\omega SST turbulence in each cell, supporting standard SST and Gamma-Re_{\theta} transition model.
	•	Blending constants \gamma_1, \gamma_2, \beta_1, \beta_2 blended using Menter function f_1:

\gamma = f_1 \gamma_1 + (1 - f_1) \gamma_2, \quad \beta = f_1 \beta_1 + (1 - f_1) \beta_2
	•	Turbulent kinetic energy production P_k:

P_k = \texttt{rev} \cdot S \cdot \frac{1}{\rho}

Limited by:

P_k \leftarrow \min(P_k, \texttt{pklim} \cdot D_k), \quad D_k = \beta^* k \omega
	•	Transition model activation:

Define wall-distance Reynolds number:

Re_w = \frac{\rho \omega y^2}{\mu}

Wake and transition suppression functions:

F_{\text{wake}} = e^{-(Re_w/10^5)^2}, \quad
F_{\theta t} = \min\left( \max\left(F_{\text{wake}} e^{-(y/\delta)^4}, 1 - \left(\frac{c \gamma - 1}{c - 1}\right)^2 \right), 1 \right)

Separation intermittency:

\gamma_{\text{sep}} = \min\left(s_1 \max\left(0, \frac{Re_S}{3.235 Re_{\theta c}} - 1 \right) F_{\text{reattach}}, 2\right) F_{\theta t}

Effective intermittency:

\gamma_{\text{eff}} = \max\left( \gamma, \gamma_{\text{sep}} \right)
	•	Modification to \omega-equation:

P_{\omega} = \frac{\gamma_{\text{eff}}}{t_{\gamma}} \left[ \max\left( 0, \frac{\gamma}{\gamma_{\text{eff}}} - c_{\text{lim}} \right) + \frac{1}{2} \right]
	•	Partial derivatives (Jacobian terms) computed consistently for implicit solves.

⸻

Jacobian Entries (approximate)

Only source terms are differentiated:

\frac{\partial S_k}{\partial k} = \beta^* \omega, \quad
\frac{\partial S_\omega}{\partial \omega} = 2 \beta \omega

Stored in local block Jacobian array qq(i,j,k,*,*).

⸻

AD-Compatible Looping

Loops are written to support forward and reverse automatic differentiation (Tapenade) using loop macros.

⸻

f1SST Subroutine: SST Blending Function f_1 Computation (module SST)

Computes f_1 over all owned cells and first-layer halos. Calculates cross-diffusion term stored in scratch(:,:,:,icd).

Computation Range Setup

Update range [iBeg:iEnd, jBeg:jEnd, kBeg:kEnd] adjusted based on boundary conditions, excluding halo layers.

Cross-Diffusion Term and f_1 Evaluation

Define:

t_1 =
\begin{cases}
\frac{\sqrt{k}}{0.09 \omega y}, & k > 0 \\
0, & \text{otherwise}
\end{cases}, \quad
t_2 = \frac{500 \nu}{\rho \omega y^2}

Compute:

t_1 \leftarrow \max(t_1, t_2)

Cross-diffusion based term:

t_2 = \frac{2k}{\max(\epsilon, CD) \cdot y^2}

where CD = \texttt{scratch(:,:,:,icd)}, \epsilon is a small constant or density-based depending on use2003SST.

Smooth minimum:

\texttt{arg1} = \min(t_1, t_2)

f_1 = \tanh(\texttt{arg1}^4)

Transition Model Modification

If transitionModel == GammaRetheta:

Re_y = \frac{\rho y \sqrt{k}}{\nu}, \quad f_3 = \exp\left( -\left( \frac{Re_y}{120} \right)^8 \right)

f_1 \leftarrow \max(f_1, f_3)

Halo Cell Update (Neumann BC)

Halo values set by copying interior values:

f_1|{\text{halo}} = f_1|{\text{adjacent interior}}

Applies to all domain boundaries.

Automatic Differentiation Support

Compatible with Tapenade reverse mode, supporting vectorized 3D loop traversal:

#ifdef TAPENADE_REVERSE
! uses 1D indexing over full 3D range
#endif

Output

Values stored in:

\texttt{scratch(i,j,k,if1SST)} = f_1

Used in SSTSource and coefficient blending.

⸻

Wall Boundary Condition Treatment for Gamma-Re_\theta Transition Model (module TurbBCRoutines)

Imposes Neumann (zero normal derivative) BCs on \gamma and Re_\theta at solid walls. Ensures no artificial flux through the wall.

Active only if transitionModel == GammaRetheta. For each block boundary face (BCFaceID(nn)), solver selects direction (iMin, iMax, etc.) and loops over face-aligned indices (icBeg:icEnd, jcBeg:jcEnd).

Ghost-cell values of \gamma and Re_\theta (iTransition1, iTransition2) are set equal to adjacent interior cell values using face-normal derivative Jacobian arrays:
	•	bmti1, bmti2 for iMin, iMax
	•	bmtj1, bmtj2 for jMin, jMax
	•	bmtk1, bmtk2 for kMin, kMax

Example on iMin face:

bmti1(i, j, iTransition1, iTransition1) = &
     bmti1(i+1, j, iTransition1, iTransition1)
bmti1(i, j, iTransition2, iTransition2) = &
     bmti1(i+1, j, iTransition2, iTransition2)

This enforces:

\frac{\partial \gamma}{\partial n} = 0, \quad \frac{\partial Re_\theta}{\partial n} = 0

Pattern repeats for all faces by copying from adjacent interior cells accordingly.

⸻

Freestream Boundary Conditions for Transition Model in initializeFlow

Active when:

transitionModel = GammaRetheta

Converts freestream turbulence intensity ratio (e.g., 0.02 for 2%) to percentage:

Tu\_inf = \texttt{turbIntensityInf} \times 100

Sets freestream values in wInf:
	•	Intermittency:

w\_inf[\texttt{iTransition1}] = 1.0
	•	Transition onset momentum thickness Reynolds number Re_{\theta t}:

w\inf[\texttt{iTransition2}] =
\begin{cases}
331.50 \cdot (Tu\infty - 0.5658)^{-0.671}, & Tu_\infty > 1.3 \\
1173.51 - 589.428 \cdot Tu_\infty + 0.2196 \cdot Tu_\infty^{-2}, & \text{otherwise}
\end{cases}

⸻

If you want me to split into sections or add any formatting details, just ask!