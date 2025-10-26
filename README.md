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

## Scientific Context

### The Quasi-Geostrophic Approximation

The QG approximation assumes a dominant balance between Coriolis and pressure gradient forces, making it ideal for studying:
- **Baroclinic instability** in stratified sheared flows
- **Synoptic-scale dynamics** (100-1000 km horizontal scales)
- **Meridional energy transport** in mid-latitude systems
- **Cyclogenesis** and frontal development

### Norwegian vs. Shapiro-Keyser Cyclone Models

The code captures the transition between two classical conceptual models of extratropical cyclone evolution:

**Norwegian Model**: Cold front overtakes warm front → occlusion → decay  
**Shapiro-Keyser Model**: Frontal fracture → bent-back warm front → warm-core seclusion

Our simulations demonstrate that these distinct pathways emerge naturally from baroclinic instability under different initial and environmental conditions, suggesting they represent different manifestations of the same underlying dynamics rather than fundamentally different physical mechanisms.

---

## Code Structure
