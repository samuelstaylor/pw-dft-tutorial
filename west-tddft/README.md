# 1D Plane-Wave Linear-Response TDDFT (TDA + Full TDDFT)

This repository contains a **toy 1D plane-wave LR-TDDFT implementation** inspired by WEST, including both the **Tamm–Dancoff approximation (TDA)** and **full TDDFT (Casida equations)**. It demonstrates:

1. Ground-state Kohn–Sham DFT.
2. Construction of the TDDFT response matrices.
3. Calculation of vertical excitation energies.
4. Comparison of DFT gaps, TDA TDDFT, and full TDDFT.
5. Simple plotting of spectra and a table of energies.

This code is designed as a **didactic tutorial** for understanding the mechanics of TDDFT.

---

## Table of Contents

- [1. Physics Background](#1-physics-background)
- [2. Equations Used](#2-equations-used)
- [3. Code Structure](#3-code-structure)
- [4. TDDFT Implementation](#4-tddft-implementation)
- [5. Plotting and Analysis](#5-plotting-and-analysis)
- [6. Answers to Key Questions](#6-answers-to-key-questions)

---

## 1. Physics Background

- **Kohn–Sham DFT** provides a way to compute the ground-state electronic density of a system by solving the effective one-particle Schrödinger equation:

$$
\hat{H}\_{KS} \psi_i(\mathbf{r}) = \epsilon_i \psi_i(\mathbf{r})
$$

where

$$

\hat{H}_{KS} = -\frac{1}{2} \nabla^2 + V_{\text{ext}}(\mathbf{r}) + V*H[n](\mathbf{r}) + V*{xc}[n](\mathbf{r})
$$

- **Linear-Response TDDFT** computes excitation energies by perturbing the ground-state density and evaluating the response. It solves **Casida's equations**:

$$
\begin{pmatrix} \mathbf{A} & \mathbf{B} \\ -\mathbf{B} & -\mathbf{A} \end{pmatrix}
\begin{pmatrix} X \\ Y \end{pmatrix} = \omega
\begin{pmatrix} X \\ Y \end{pmatrix}
$$

where the matrices are:

$$

A*{ia,jb} = (\epsilon_a - \epsilon_i)\delta*{ij}\delta*{ab} + \langle ia | f*{\text{Hxc}} | jb \rangle
$$

$$
B*{ia,jb} = \langle ia | f*{\text{Hxc}} | jb \rangle
$$

- **Tamm–Dancoff approximation (TDA)** simplifies TDDFT by ignoring the \(\mathbf{B}\) matrix:

$$

\mathbf{A} X = \omega X
$$

---

## 2. Equations Used

1. **Hartree potential (1D, reciprocal space):**

$$
V_H(x) = \text{Re}\Big[ \text{IFFT}\Big( \frac{4 \pi n(G)}{G^2} \Big) \Big]
$$

2. **LDA Exchange potential (toy 1D model):**

$$

V*{xc}(n) = - \frac{4}{3} A*{xc} n^{1/3}
$$

3. **TDDFT kernel (functional derivative):**

$$
f*{xc} = \frac{\delta V*{xc}}{\delta n} = - \frac{4}{9} A\_{xc} n^{-2/3}
$$

4. **Transition density for matrix elements:**

$$
\rho_{ia}(x) = \psi _i^*(x) \psi _a(x)
$$

$$
\langle ia | f*{\text{Hxc}} | jb \rangle = \int dx \, \rho*{ia}(x) f*{\text{Hxc}}(x) \rho*{jb}^*(x)
$$

5. **TDA excitation energies:**

$$
\omega\_{\text{TDA}} = \text{eig}(\mathbf{A})
$$

6. **Full TDDFT excitation energies:**

$$
\omega\_{\text{full}} = \sqrt{\text{eig}((\mathbf{A}-\mathbf{B})(\mathbf{A}+\mathbf{B}))}
$$

---

## 3. Code Structure

1. **Simulation Parameters**: electron count, box length, grid size, plane-wave cutoff.
2. **Reciprocal Space**: generate G-vectors and truncate plane waves by kinetic energy.
3. **External Potential**: Gaussian well to bind electrons.
4. **Hartree + XC**: functions for Hartree potential, LDA exchange, and TDDFT kernel.
5. **KS Hamiltonian**: build plane-wave matrix for kinetic + potential energy.
6. **Density**: compute real-space density from KS orbitals.
7. **SCF Loop**: self-consistent solution for ground-state density and KS eigenvalues.
8. **TDDFT Matrices**: construct TDA (A) and full TDDFT (A, B) matrices using transition densities.
9. **Solve TDDFT**: diagonalize TDA or full Casida matrices for excitation energies.
10. **Plotting & Table**: visualize DFT gaps, TDA, and full TDDFT excitations; print table.

---

## 4. TDDFT Implementation

- **TDA**: ignore coupling between excitations and de-excitations. Diagonalize \(\mathbf{A}\).
- **Full TDDFT**: include both \(\mathbf{A}\) and \(\mathbf{B}\), solve \((A-B)(A+B)X = \omega^2 X\).
- **Physics interpretation**: TDDFT captures **electron-hole interactions** and **density response** to a time-dependent perturbation, going beyond simple KS orbital differences.

---

## 5. Plotting and Analysis

- **DFT HOMO-LUMO gaps**: black circles, represent independent particle gaps.
- **TDA TDDFT**: blue triangles, excitation energies using linear-response with TDA.
- **Full TDDFT**: red squares, excitation energies including full coupling.
- First 5 excitations plotted for clarity.

---

## 6. Answers to Key Questions

**Q1: Why are there unoccupied energies and states calculated with DFT?**

- DFT provides all eigenvalues and orbitals, including unoccupied (virtual) states. These virtual orbitals are needed for **excited-state calculations** in TDDFT because excitations correspond to transitions from occupied (HOMO) to unoccupied (LUMO + higher) orbitals.

**Q2: Why is DFT LUMO gap different from TDDFT?**

- DFT gaps are **independent-particle KS gaps** and do not include **electron-hole interactions**. TDDFT incorporates the **response kernel f_Hxc**, which accounts for Coulomb and exchange-correlation effects in excitations, giving more accurate excitation energies.

**Q3: Why are TDDFT TDA and no-TDA different?**

- TDA neglects de-excitation (B matrix) contributions. Full TDDFT includes **both excitation and de-excitation couplings**, which slightly shifts excitation energies. In some systems, the difference is small, in others it can be significant.

**Q4: Why are all their energy values so close in this case?**

- The toy 1D system has **small electron-hole interaction effects** and a simple local potential. Therefore, the KS gaps, TDA, and full TDDFT excitations are similar.

**Q5: Where does the time-dependent part of TDDFT come in?**

- TDDFT solves the **time-dependent Kohn–Sham linear-response problem**. The matrices A and B are derived from the **density response to a time-dependent perturbation**. The eigenvalues ω correspond to frequencies of density oscillations, i.e., excitation energies.

**Q6: What is the main thing we are solving in TDDFT?**

- We are solving for the **linear response of the electron density** to small time-dependent perturbations, which yields **vertical excitation energies** and includes **electron-hole interactions** beyond the independent-particle picture.

---

### References / Suggested Reading

1. E. Runge, E. K. U. Gross, _Density-Functional Theory for Time-Dependent Systems_, Phys. Rev. Lett. 52, 997 (1984).
2. M. Casida, _Time-Dependent Density Functional Response Theory for Molecules_, in _Recent Developments and Applications in Modern DFT_, Elsevier (1995).
3. WEST code: [https://west-code.org](https://west-code.org) – inspiration for plane-wave LR-TDDFT implementations.

---

$$
$$
