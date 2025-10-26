# QG_planar_front

Quasi-Geostrophic Pseudo-Spectral Code for Simulating Atmospheric Fronts and Cyclogenesis  
Authors: Antonio Segalini & Juan Vázquez Portillo 
Department of Earth Sciences, Uppsala University

---

## Overview

**QG_planar_front** is a computational framework designed to simulate atmospheric fronts using the quasi-geostrophic (QG) approximation. This repository contains the numerical implementation accompanying the paper *"Quasi-geostrophic simulation of atmospheric fronts"* published in the Journal of Advances in Modeling Earth Systems (JAMES).

The model employs a pseudo-spectral approach to solve the quasi-geostrophic equations on a β-plane, capturing the dynamics of frontogenesis and cyclogenesis in mid-latitude weather systems. Despite the simplifications inherent to the QG framework, the code successfully reproduces key features observed in extratropical cyclones, including the transition between **Norwegian** and **Shapiro-Keyser** cyclone models.

This implementation is particularly valuable for researchers interested in atmospheric dynamics, as it offers **computational efficiency** compared to primitive equation (PE) models while retaining the essential physics governing baroclinic instability and frontal evolution. The model's efficiency makes it suitable for conducting large ensembles of simulations, sensitivity studies, and climate change investigations.

<p align="center">
  <img src="https://media3.giphy.com/media/v1.Y2lkPTc5MGI3NjExamRoMXN1d2d3aWRrZmUzMnJnN3pvMDA0dDI5M3E3cjhpM2VkbWc3MCZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/Ix2zxUCFpb3t9rH5tX/giphy.gif" width="500">
</p>

---

## Key Features

- **Quasi-Geostrophic Framework**: Efficiently filters high-frequency gravity waves while preserving synoptic-scale dynamics
- **Pseudo-Spectral Methods**: Fourier transforms in horizontal directions, Chebyshev polynomials in vertical direction
- **Frontogenesis Simulation**: Captures the formation and intensification of atmospheric fronts from baroclinic instability
- **Cyclone Life Cycles**: Reproduces both Norwegian (occlusion) and Shapiro-Keyser (frontal fracture) cyclone evolution patterns
- **MPI Parallelization**: Modern Fortran implementation with efficient parallel computing capabilities
- **Configurable Physics**: Adjustable damping, hyperviscosity, grid resolution, and initial perturbation parameters
- **Spectral Dealiasing**: 3/2-rule implementation prevents aliasing errors in nonlinear computations
- **Validated Against PE Models**: Results compared with established primitive equation simulations from literature

---

## Code Structure

### `QG_plane_front.f90`
The main solver program that:
- Sets up the quasi-geostrophic problem with β-plane geometry and baroclinic shear.
- Initializes the computational grid (Fourier in horizontal, Chebyshev in vertical) and atmospheric profiles.
- Solves the potential vorticity equation using pseudo-spectral methods with MPI parallelization.
- Implements time integration via Euler (first step) and Leapfrog schemes with semi-implicit treatment of dissipation.
- Handles PV inversion to compute streamfunction and velocity fields from PV anomalies.
- Outputs velocity and vorticity fields at specified intervals.

**Key Editable Parameters are found in `INPUT.txt`**
Inside this program:
- `Ujet`: **Jet velocity** (currently hardcoded as 60 m/s)
- `Tdamping_days`: **Damping timescale in days** (currently hardcoded as 2.7557)

### `INPUT.txt`

All simulation parameters are specified in `INPUT.txt`. Below is a detailed description of each parameter:

| Parameter | Type | Unit | Description |
|-----------|------|------|-------------|
| `fileName_root` | string | - | Output directory name for results |
| `Nx` | integer | - | Horizontal grid resolution in x-direction (zonal) |
| `Ny` | integer | - | Horizontal grid resolution in y-direction (meridional) |
| `Nz` | integer | - | Vertical grid resolution |
| `Lx` | real | m | Domain size in x-direction (zonal) |
| `Ly` | real | m | Domain size in y-direction (meridional) |
| `Lz` | real | m | Domain height (vertical extent) |
| `DT` | real | s | Time step for integration |
| `f0` | real | 1/s | Coriolis parameter at reference latitude |
| `beta` | real | 1/(m·s) | Meridional gradient of Coriolis parameter (β-plane) |
| `nu` | real | m²ˢ/s | Hyperviscosity coefficient |
| `p` | integer | - | Order of hyperviscosity (typically 4 or 8; ∇²ᵖ dissipation) |
| `Initial_perturbation` | real | - | Amplitude of initial Gaussian perturbation (dimensionless, 0-1) |
| `kx_initial` | integer | - | Zonal wavenumber of initial perturbation |
| `ky_initial` | integer | - | Meridional wavenumber of initial perturbation |
| `Nsave` | integer | - | Save output every Nsave timesteps |
| `Tfinal` | real | s | Total simulation time |
| `dealiasing_condition` | logical | - | Enable spectral dealiasing (.true. or .false.) |
| `nonlinear_condition` | logical | - | Enable nonlinear terms (.true.) or linear mode (.false.) |
| `friction_bottom` | real | - | Surface friction coefficient at bottom boundary |

**The three different cases treated in the scientific article can be found at `INPUT_Run1.txt`, `INPUT_Run2.txt` and `INPUT_Run3.txt`.**

### `PROFILES.py`
A Python utility for generating atmospheric background profiles and computing the model's characteristic timescales.

**Functions:**
- Generates vertical profiles of background zonal wind (U), density (ρ), and buoyancy frequency (N²).
- Computes vertical derivatives required by the QG model (dU/dz, d²U/dz², dρ/dz, dN²/dz).
- Exports profiles to `PROFILES.txt` for use in the Fortran solver.
- Calculates the **model timescale (Ts)** from physical parameters.

**Critical Note on Parameter Linkage:**

In the QG framework, **parameters are interdependent** and cannot be freely chosen independently. The model's characteristic timescale is governed by:

\[
T_s = \frac{\eta \cdot N \cdot H}{f \cdot U_{max}}
\]

where:
- η: Non-dimensional scaling factor
- N: Buoyancy frequency
- H: Domain height
- f: Coriolis parameter
- U_max: Maximum jet velocity

**Important:** Users must ensure consistency between:
1. `N` and `N2` in the profile
2. `Umax` in this script (currently hardcoded as 60 m/s)
3. `LAMBDA` in this script (must match Fortran code: currently 10.0/Lz)
4. `f0` and `Lz` in the Fortran INPUT.txt

Changing any of these independently will break the model balance and lead to unrealistic dynamics. The output timescale calculation needs to be updated at `QG_plane_front.f90` for consistency.

---
