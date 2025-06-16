# Jacobian-Vector Product Verification

**File:** `compute_jacobian_vector_product_fwd.py`

## Purpose

Verify the correctness of Jacobian-vector products (JVP) computed via Automatic Differentiation (AD) by comparing them against Finite Difference (FD) and Complex Step (CS) approximations for different perturbation modes in the CFD solver.

---

## Perturbation Modes

- `wDot` : State vector perturbation
- `xVDot` : Vertex coordinate perturbation
- `xDvDot` : Design variable perturbation (one variable at a time)
- `xSDot` : Surface coordinate perturbation

---

## Methodology

For each perturbation mode:

1. **Compute JVP using AD:**  
   Call `computeJacobianVectorProductFwd` with `mode` set to the perturbation type and `mode_vector` set to the perturbation vector.

2. **Compute JVP using FD:**  
   Use `computeJacobianVectorProductFwd` with `mode="FD"` for finite difference approximation.

3. **(Optional) Compute JVP using CS:**  
   If solver supports complex arithmetic (`dtype == "D"`), compute JVP using `mode="CS"`. Exceptions during CS computation are ignored.

4. **Calculate Norms:**  
   Compute norms of residuals, functions, and forces for AD results and store them.

5. **Calculate Relative Errors:**  
   Compute relative errors between AD and FD results for residuals, functions, and forces.

6. **Store Results:**  
   Use a handler object to log norms and relative errors with appropriate tolerance settings.

---

## Key Functions

- `computeJacobianVectorProductFwd(x, mode_vector, mode)`:  
  Returns JVP outputs (`R` residuals, `F` functions, `Jf` forces) for given perturbation `mode` and vector.

- `handler.par_add_norm(label, value, rtol, atol)`:  
  Store a norm value with relative and absolute tolerance.

- `handler.par_add_rel_error(label, rel_error, rtol, atol)`:  
  Store a relative error value with tolerances.

- `handler.root_add_dict(label, dict_values, rtol, atol)`:  
  Store dictionary of function outputs with tolerances.

- `handler.root_add_rel_error_dict(label, rel_error_dict, rtol, atol)`:  
  Store relative errors for multiple function components.

- `handler.rel_error(val1, val2)`:  
  Computes relative error between two values or arrays.

---

## Data Output

- On process rank 0, results including detailed norms, function values, relative errors, and metadata are saved as JSON files under the `jac_check` directory.

---

## Notes

- FD and CS computations serve as references to validate the AD-based JVP.
- CS mode is only computed if complex arithmetic is enabled.
- The verification process handles multiple design variables independently when using `xDvDot`.