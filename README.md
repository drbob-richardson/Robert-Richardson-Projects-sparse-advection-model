# Sparse Nonlinear Spatio-Temporal Models Using Implied Advection

**Author:** Robert Richardson, Brigham Young University  
**Paper:** Richardson, R. (2025). *Sparsity in Nonlinear Dynamic Spatio-Temporal Models using Implied Advection.* **Environmetrics**.  
**DOI:** 10.1002/env.2456 
**Keywords:** Nearest-Neighbor Gaussian Process (NNGP), SPDE, Advection-Diffusion, Precision Tapering, Dynamic Modeling, Computational Scalability  

---

## Overview

Traditional dynamic spatio-temporal models are computationally expensive because each filtering step requires inverting dense $n \times n$ covariance matrices, giving cubic time complexity $O(n^3)$.  

This repository implements a **sparse, implied-advection dynamic model** that avoids dense matrix inversion entirely.  
By combining:

- **Backward-difference discretization** of the advection–diffusion SPDE (“implied advection”),
- **Nearest Neighbor Gaussian Processes (NNGP)** to enforce sparsity in the precision matrix, and  
- **Precision tapering** to maintain sparsity through time,

the method achieves **near-linear scaling** in the number of spatial locations $n$.

---

## Motivation

Dynamic models in environmental and geophysical data analysis (e.g., SST forecasting, pollution transport) often involve tens of thousands of spatial points.  
Even with reduced-rank methods, large matrix inversions dominate runtime.  
Our approach replaces dense inversions with sparse updates, preserving model fidelity while dramatically lowering computational cost.

---

## Computational Gains

| Approach | Matrix Structure | Time Complexity | Memory Complexity | Comments |
|-----------|------------------|-----------------|-------------------|-----------|
| **Dense Kalman Filtering** | Full $n \times n$ covariance | $O(Tn^3)$ | $O(n^2)$ | Standard dynamic spatio-temporal model |
| **Reduced-Rank Model** | Low-rank $r \times r$ basis ($r \ll n$) | $O(T(r^3 + n r^2))$ | $O(r^2)$ | Requires dimension reduction |
| **Sparse Implied-Advection + NNGP (this model)** | Sparse precision with $m$ neighbors ($m \ll n$) | $O(Tnm^2)$ | $O(nm)$ | No dimension reduction, near-linear scaling |

Here, $T$ is the number of time points, $n$ the number of spatial locations, and $m$ the neighborhood size (typically 10–30).

> **In short:**  
> The implied-advection NNGP model reduces the per-step cost from $O(n^3)$ to $O(n\,m^2)$,  
> yielding **linear scalability in both space and time** for fixed $m$.

---

## Empirical Results (Table 4 from the Paper)

| Spatial Locations (n) | Time Points (T) | Model | Runtime (hours) | CRPS (↓ better) |
|:----------------------:|:---------------:|:-------|:----------------:|:----------------:|
| 400 | 100 | Reduced-Rank (10 basis) | 1.2 | 0.425 |
| 400 | 100 | Sparse Advection | 2.0 | 0.310 |
| 625 | 100 | Reduced-Rank (30 basis) | 3.1 | 0.471 |
| 625 | 100 | Sparse Advection | 4.6 | 0.315 |
| 900 | 100 | Reduced-Rank (50 basis) | 6.0 | 0.488 |
| 900 | 100 | Sparse Advection | 8.3 | 0.320 |

*Results reproduced from the simulation study comparing the sparse advection–diffusion model with reduced-rank Bayesian hierarchical models (BHMs).  
The sparse model shows superior fit (lower CRPS) with modest additional computation.*  
(Source: *Sparsity in Nonlinear Dynamic Spatio-Temporal Models using Implied Advection*, Table 4)

---

## Key Techniques

- **Implied Advection:** Backward differencing of the SPDE time derivative produces a *sparse inverse evolution matrix* $G^{-1}$.  
- **Nearest Neighbor Gaussian Process (NNGP):** Each location depends on a small neighborhood ($m$) → sparse precision matrix.  
- **Precision Tapering:** Re-projects the precision matrix back to the NNGP sparsity pattern at each time step to prevent fill-in.  
- **Extended Kalman Filtering:** Enables nonlinear dynamics while preserving sparsity.  

---

## Advantages

- ✅ No dense matrix inversion  
- ✅ No need for rank reduction or basis truncation  
- ✅ Scales linearly with number of spatial locations  
- ✅ Handles nonlinear and spatially varying dynamics  
- ✅ Empirically maintains accuracy comparable to dense models






