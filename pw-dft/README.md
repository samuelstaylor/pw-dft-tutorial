# 1D Plane-Wave Density Functional Theory (DFT) — Toy Model

This repository contains a **minimal 1D Kohn–Sham DFT solver** using a plane-wave basis. It is designed for educational purposes to help you understand the structure of plane-wave DFT codes, like Quantum ESPRESSO, in a simplified and hackable way.

---

## Goals

- Understand how plane waves, FFTs, and self-consistent field (SCF) iterations work together.
- Reproduce the structure (not full realism) of plane-wave DFT.
- Keep the implementation minimal and explicit for learning purposes.

---

## Table of Contents

- [1. What We Include](#1-what-we-include)
- [2. What We Omit (on purpose)](#2-what-we-omit-on-purpose)
- [3. Physics Background](#3-physics-background)
- [4. Equations Used](#4-equations-used)
- [5. Code Structure](#5-code-structure)
- [6. Visualization](#6-visualization)
- [7. Convergence and Output](#7-convergence-and-output)

---

## 1. What We Include

1. 1D periodic supercell
2. Plane-wave basis with kinetic energy cutoff
3. Local external potential (toy pseudopotential)
4. Hartree potential (1D Poisson equation, periodic)
5. LDA exchange-correlation (toy functional)
6. Self-consistent field (SCF) loop

---

## 2. What We Omit (on purpose)

- Nonlocal pseudopotentials
- k-point sampling (Gamma only)
- Spin and magnetism
- Advanced mixing schemes (Broyden, DIIS)
- Realistic exchange-correlation parameterizations

This allows the code to remain **educational, minimal, and hackable**.

---

## 3. Physics Background

- **Kohn–Sham DFT** solves for a system of non-interacting electrons that reproduce the exact ground-state density \(n(x)\) of an interacting system:

\[
\hat{H}\_{KS} \psi_i(x) = \epsilon_i \psi_i(x)
\]

with the effective Hamiltonian:

\[
\hat{H}_{KS} = -\frac{1}{2}\frac{d^2}{dx^2} + V_{\text{ext}}(x) + V*H[n](x) + V*{xc}[n](x)
\]

- **Plane-wave expansion** of Kohn–Sham orbitals:

\[
\psi*i(x) = \sum_G C*{iG} e^{i G x}
\]

where \(G = \frac{2\pi}{L} k\) are the reciprocal lattice vectors of the 1D supercell.

- **Self-consistent field (SCF) loop** iteratively updates the density:

\[
n(x) = \sum_i f_i |\psi_i(x)|^2
\]

until convergence is reached.

---

## 4. Equations Used

1. **External potential (toy Gaussian pseudopotential)**:

\[
V\_{\text{ext}}(x) = V_0 \exp\left[-\frac{(x-x_0)^2}{2\sigma^2}\right]
\]

2. **Hartree potential (1D Poisson equation)**:

\[
\frac{d^2 V_H(x)}{dx^2} = -4\pi n(x)
\]

in reciprocal space:

\[
V_H(G) = \frac{4 \pi n(G)}{G^2}, \quad G \neq 0
\]

3. **LDA-like exchange-correlation potential (toy functional)**:

\[
V*{xc}(n) = - \frac{4}{3} A*{xc} n^{1/3}
\]

4. **Density from Kohn–Sham orbitals**:

\[
n(x) = \sum_i f_i |\psi_i(x)|^2
\]

5. **Kohn–Sham Hamiltonian in plane-wave basis**:

\[
H*{GG'} = \frac{G^2}{2} \delta*{GG'} + V\_{\text{eff}}(G-G')
\]

---

## 5. Code Structure

1. **Import Libraries**: `numpy`, `scipy.linalg`, `matplotlib`, `fft` routines.
2. **Simulation Cell**: define number of electrons `Ne`, supercell length `L`, real-space grid `Nx`, and reciprocal-space grid `G`.
3. **Plane-Wave Cutoff**: select plane waves satisfying \( \frac{G^2}{2} \le E\_{\text{cut}} \).
4. **External Potential**: define Gaussian pseudopotential.
5. **Hartree Potential**: solve 1D Poisson equation in reciprocal space.
6. **Exchange-Correlation Potential**: simple LDA functional.
7. **Kohn–Sham Hamiltonian**: build Hamiltonian in plane-wave basis.
8. **Density from Orbitals**: compute real-space density from plane-wave coefficients.
9. **SCF Loop**: iterate until density convergence.
10. **Output**: print converged density and energies, plot density and potential components.

---

## 6. Visualization

- **Plane-wave grid in reciprocal space**: visualize selected G vectors.
- **Real-space plane-wave functions**: see oscillatory behavior of basis functions.
- **Self-consistent density**: plot converged \(n(x)\) along with scaled \(V\_{\text{ext}}\).
- **Potential components**: \(V*{\text{ext}}, V_H, V*{xc}, V\_{\text{eff}}\) plotted together for comparison.

---

## 7. Convergence and Output

- SCF convergence tracked using density difference:

\[
\Delta n = \| n*{\text{new}} - n*{\text{old}} \|\_2
\]

- Convergence typically achieved in ~40 iterations for the toy system.
- Density table printed QE-style: `x (a.u.)` vs `n(x)` for visualization and debugging.

---

### References / Suggested Reading

1. R. M. Martin, _Electronic Structure: Basic Theory and Practical Methods_, Cambridge University Press (2004).
2. R. Car and M. Parrinello, _Unified Approach for Molecular Dynamics and DFT_, Phys. Rev. Lett. 55, 2471 (1985).
3. Quantum ESPRESSO: [https://www.quantum-espresso.org](https://www.quantum-espresso.org) — for real plane-wave DFT applications.

---
