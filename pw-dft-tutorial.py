
# 1D Plane-Wave Density Functional Theory (DFT) — Toy Model
# --------------------------------------------------------
# This notebook implements a *minimal*, educational 1D Kohn–Sham DFT solver
# using a plane-wave basis, closely mirroring the logic used in
# plane-wave pseudopotential codes like Quantum ESPRESSO.
#
# GOAL:
# - Understand how plane waves, FFTs, and self-consistency work together
# - Reproduce the *structure* (not the full realism) of PW-DFT
# - Keep everything simple, explicit, and hackable
#
# WHAT WE INCLUDE:
# - 1D periodic supercell
# - Plane-wave basis with kinetic energy cutoff
# - Local external potential (toy pseudopotential)
# - Hartree potential (1D Poisson, periodic)
# - LDA exchange-correlation (toy functional)
# - Self-consistent field (SCF) loop
#
# WHAT WE OMIT (on purpose):
# - Nonlocal pseudopotentials
# - k-point sampling (Gamma only)
# - Spin, magnetism
# - Advanced mixing (Broyden, DIIS)
# - Realistic XC parameterizations
#
# This is *conceptually faithful* to QE, but numerically tiny.

import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.linalg import eigh


# -----------------------------
# 1. Parameters, Simulation Cell, Grids, and Initialization
# -----------------------------

TURN_ON_HARTREE = True
TURN_ON_XC = True

Ne = 162  # number of electrons

def neutral_ground_state_occupation(Ne):
    n_fully_occupied_band  = Ne // 2
    n_singly_occupied_band = Ne % 2
    n_bands = n_fully_occupied_band + n_singly_occupied_band # total number of occupied bands
    occ = [2.0] * n_fully_occupied_band + [1.0] * n_singly_occupied_band
    return occ
occ = neutral_ground_state_occupation(Ne)
print("Occupation numbers:", occ)


# Length of 1D supercell (periodic)
L = 20.0  # atomic units

# Real-space grid (used for density & potentials)
Nx = 512
x = np.linspace(0, L, Nx, endpoint=False)
dx = L / Nx

# Reciprocal-space grid (plane waves)
G = 2 * np.pi * fftfreq(Nx, d=dx)


# -----------------------------
# 2. Plane-Wave Cutoff
# -----------------------------

# Kinetic energy cutoff (in a.u.)
Ecut = 80.0 # max is about 260 for Nx=512
            # 80 means max of 162 electrons

# Select plane waves satisfying |G|^2 / 2 <= Ecut
pw_mask = 0.5 * G**2 <= Ecut
G_pw = G[pw_mask]
Np = len(G_pw)

print(f"Number of plane waves: {Np}")


# -----------------------------
# 3. External Potential
# -----------------------------
# This plays the role of a *local pseudopotential* in QE.
# We use a smooth Gaussian well centered in the cell.

V0 = -5.0
sigma = 0.5
x0 = L / 2

V_ext = V0 * np.exp(-((x - x0) ** 2) / (2 * sigma**2))

# Fourier transform of external potential
V_ext_G = fft(V_ext) / Nx


# -----------------------------
# 4. Hartree Potential (1D)
# -----------------------------
# Solve Poisson equation in reciprocal space:
#   d^2 V_H / dx^2 = -4π n(x)
# → V_H(G) = 4π n(G) / G^2   (G ≠ 0)

def hartree_potential(n):
    nG = fft(n)
    VH_G = np.zeros_like(nG, dtype=complex)
    if not TURN_ON_HARTREE:
        return np.zeros_like(n)
    for i, g in enumerate(G):
        if abs(g) > 1e-12:
            VH_G[i] = 4.0 * np.pi * nG[i] / (g**2)  # factor 4π for a.u.
    VH = np.real(ifft(VH_G))
    return VH


# -----------------------------
# 5. Exchange-Correlation (LDA)
# -----------------------------
# Extremely simple LDA-like model:
#   Exc[n] = -A ∫ n^(4/3) dx
#   Vxc = dExc/dn = -(4/3) A n^(1/3)

A_xc = 0.75

def v_xc(n):
    if TURN_ON_XC:
        return -(4/3) * A_xc * np.maximum(n, 1e-10) ** (1/3)
        # TEMPORARY FOR DEBUGGING
    return np.zeros_like(n)

# -----------------------------
# 6. Kohn–Sham Hamiltonian
# -----------------------------
# Construct H(G,G') = T(G) δ + V_eff(G-G')

def ks_hamiltonian(Veff):
    Veff_G = fft(Veff) / Nx

    H = np.zeros((Np, Np), dtype=complex)

    for i, Gi in enumerate(G_pw):
        # Kinetic energy (diagonal)
        H[i, i] = 0.5 * Gi**2

        for j, Gj in enumerate(G_pw):
            dG = Gi - Gj
            idx = np.where(np.isclose(G, dG))[0]
            if idx.size > 0:
                H[i, j] += Veff_G[idx[0]]

    return H


# -----------------------------
# 7. Density from Orbitals
# -----------------------------

def density_from_orbitals(C, occ):
    n = np.zeros(Nx)
    for i in range(len(occ)):
        psiG = np.zeros(Nx, dtype=complex)
        psiG[pw_mask] = C[:, i]
        psi = ifft(psiG)
        n += occ[i] * np.abs(psi) ** 2
    return n


# -----------------------------
# 8. Self-Consistent Field Loop
# -----------------------------

# Initial guess: uniform density
n = np.ones(Nx) * (Ne / L)

alpha = 0.3  # mixing parameter
n_iter = 50

# Convergence parameters
conv_thr = 1e-6

for it in range(n_iter):
    VH = hartree_potential(n)
    Vxc = v_xc(n)

    Veff = V_ext + VH + Vxc

    H = ks_hamiltonian(Veff)
    eigvals, eigvecs = eigh(H)

    n_new = density_from_orbitals(eigvecs, occ)

    # Density difference (L2 norm)
    delta = np.linalg.norm(n_new - n) * dx

    # Simple linear mixing
    n = alpha * n_new + (1 - alpha) * n

    print(f"Iteration {it:3d}: E0 = {eigvals[0]:.6f}   delta_n = {delta:.3e}")

    if delta < conv_thr:
        print(f"SCF converged in {it+1} iterations (delta_n < {conv_thr})")
        break


# -----------------------------
# 9. Final Results
# -----------------------------

import matplotlib.pyplot as plt

plt.figure(figsize=(8,4))
plt.plot(x, n, label="Electron density")
plt.plot(x, V_ext / np.max(np.abs(V_ext)) * np.max(n), '--', label="V_ext (scaled)")
plt.legend()
plt.xlabel("x")
plt.ylabel("Density")
plt.title("1D Plane-Wave DFT: Self-Consistent Density")
plt.tight_layout()
plt.show()


# -----------------------------
# 10. Plot Kohn–Sham Potentials
# -----------------------------

plt.figure(figsize=(8,5))

plt.plot(x, V_ext, label="V_ext")
plt.plot(x, VH, label="V_Hartree")
plt.plot(x, Vxc, label="V_xc")
# Scale Hartree and XC potentials for visibility
#plt.plot(x, VH / np.max(np.abs(VH)), label="V_Hartree (scaled)")
#plt.plot(x, Vxc / np.max(np.abs(Vxc)), label="V_xc (scaled)")
plt.plot(x, Veff, "--", linewidth=2, label="V_eff = V_ext + V_H + V_xc")

plt.xlabel("x")
plt.ylabel("Potential (a.u.)")
plt.title("Kohn–Sham Potential Components")
plt.legend()
plt.tight_layout()
plt.show()



# -----------------------------
# 11. Print Density Table
# -----------------------------
# Print x and n(x) side-by-side, similar to a QE-style output table


print("x (a.u.) n(x)")
print("---------------------------")
for i in range(0, Nx, 10):
    print(f"{x[i]:10.6f} {n[i]:12.8f}")