# 1D Shallow Water Equations: Finite Difference Solvers

## Overview
This repository contains a MATLAB implementation of Finite Difference methods for solving the one-dimensional Shallow Water Equations (SWE). The SWE represent a system of non-linear hyperbolic conservation laws used to model fluid flow below a free surface. 

The state vector $q=(h,m)^{T}$ consists of the water depth $h(x,t)$ and the horizontal discharge $m(x,t)$. The governing system is formulated as:

$$\begin{bmatrix}h\\ m\end{bmatrix}_{t}+\begin{bmatrix}m\\ \frac{m^{2}}{h}+\frac{1}{2}gh^{2}\end{bmatrix}_{x}=S(x,t)$$

## Numerical Methods Implemented
The project explores the stability, accuracy, and behavior of two classic explicit finite difference schemes:
* **Lax-Friedrichs Scheme:** A first-order method known for its robust stability but high numerical dissipation (smearing) around discontinuities.
* **Lax-Wendroff Scheme:** A second-order method that provides higher accuracy for smooth flows but introduces numerical dispersion (spurious oscillations or Gibbs phenomenon) near shocks.

Both solvers adaptively compute the time-step $k$ at each iteration using the Courant-Friedrichs-Lewy (CFL) condition to ensure mathematical stability:

$$k=CFL\frac{\Delta x}{\max_{i}(|u_{i}|+\sqrt{gh_{i}})}$$

## Experiments & Project Structure
Run the numbered MATLAB scripts sequentially to reproduce the numerical experiments:

* **`01_analytical_convergence_test.m`**: Validates the Lax-Friedrichs implementation against a manufactured exact analytical solution using a specific source term $S(x,t)$. It outputs a log-log plot demonstrating the convergence rate of the $L^{\infty}$ error as the spatial step $\Delta x$ decreases.
* **`02_homogeneous_smooth_waves.m`**: Solves the homogeneous SWE ($S=0$) using smooth sinusoidal initial conditions and periodic boundaries. It compares the numerical results against a high-resolution reference mesh to evaluate performance on smooth wave propagation.
* **`03_shock_wave_riemann_problem.m`**: Investigates a Riemann-type problem featuring a discontinuous initial discharge ($m$) simulating a sudden dam break or step change. This script highlights the fundamental differences between the two schemes:
  * Demonstrates the diffusive nature of **Lax-Friedrichs**.
  * Demonstrates the dispersive oscillations of **Lax-Wendroff** around the resulting shock waves.
