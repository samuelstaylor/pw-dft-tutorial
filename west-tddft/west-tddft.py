"""
=============================================================
1D Plane-Wave LR-TDDFT (Tamm–Dancoff Approximation + Full TDDFT)
Toy Model Inspired by WEST Implementation

This code:
1) Performs ground-state Kohn–Sham DFT
2) Builds the LR-TDDFT matrices (TDA and full form)
3) Solves for vertical excitation energies
4) Plots spectra and prints a table of energies

Implements:
A_{ia,jb} = (eps_a - eps_i) δ_{ij}δ_{ab} + < ia | f_Hxc | jb >
B_{ia,jb} = < ia | f_Hxc | jb >

Only Γ-point
Only local Hartree + LDA kernel
No hybrid exchange
No spin
=============================================================
"""

import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# Conversion factor: 1 Hartree (atomic unit of energy) = 27.2114 eV
# All internal calculations are done in Hartree atomic units.
Hartree2eV = 27.2114

# ============================================================
# 1. Simulation Parameters
# ============================================================

Ne = 6                     # number of electrons (closed shell)
# We assume spin-restricted DFT: each spatial orbital holds 2 electrons.

L = 20.0                   # supercell length (a.u.)
# Defines the 1D simulation box (periodic) in Bohr units.
# 1 Bohr ≈ 0.529 Å → L = 20 Bohr ≈ 10.6 Å.

Nx = 512                   # number of real-space grid points
dx = L / Nx                # spacing between grid points
x = np.linspace(0, L, Nx, endpoint=False)
# Real-space coordinate grid

Ecut = 60.0                # plane-wave kinetic energy cutoff (Hartree)
max_scf_iter = 50          # maximum self-consistent field iterations
mixing = 0.3               # linear density mixing for SCF convergence

# ============================================================
# 2. Reciprocal Space Grid
# ============================================================

G = 2*np.pi*fftfreq(Nx, d=L/Nx)   # reciprocal lattice vectors
pw_mask = 0.5 * G**2 <= Ecut      # mask plane waves below cutoff
G_pw = G[pw_mask]                  # plane waves included in KS basis
Np = len(G_pw)                     # number of plane waves used
# This defines the size of the KS Hamiltonian

# ============================================================
# 3. External Potential (Gaussian well)
# ============================================================

V0 = -5.0       # depth of Gaussian potential
sigma = 0.7     # width of Gaussian
x0 = L/2        # center of Gaussian

V_ext = V0 * np.exp(-(x-x0)**2/(2*sigma**2))
# Attractive Gaussian well to bind electrons and create discrete bound states

# ============================================================
# 4. Hartree + XC (Local LDA)
# ============================================================

A_xc = 0.75  # exchange constant for toy 1D LDA

def hartree_potential(n):
    """Compute Hartree potential from density using FFT

    Steps:
    1) Transform density to reciprocal space
    2) Multiply by Coulomb kernel V(G) = 4π / G^2 (G != 0)
    3) Transform back to real space
    """
    nG = fft(n)
    VH_G = np.zeros_like(nG, dtype=complex)
    for i, g in enumerate(G):
        if abs(g) > 1e-12:  # skip G=0 to avoid singularity
            VH_G[i] = 4*np.pi * nG[i] / (g**2)
    return np.real(ifft(VH_G))  # imaginary part is numerical noise

def v_xc(n):
    """Local density approximation (LDA) exchange potential
    v_xc = δE_x/δn ~ n^(1/3) in 1D toy model
    """
    return -(4/3)*A_xc*np.maximum(n,1e-10)**(1/3)

def f_xc_kernel(n):
    """Derivative of v_xc wrt density for TDDFT kernel
    f_xc = dv_xc/dn
    This enters the response kernel for linear-response TDDFT
    """
    return -(4/9)*A_xc*np.maximum(n,1e-10)**(-2/3)

# ============================================================
# 5. Build KS Hamiltonian
# ============================================================

def ks_hamiltonian(Veff):
    """Construct plane-wave KS Hamiltonian H = T + Veff
    
    In plane-wave basis:
    - T(G) = 1/2 G^2 (diagonal kinetic energy)
    - Veff couples different G vectors (off-diagonal)
    """
    Veff_G = fft(Veff)/Nx  # Fourier transform of effective potential
    H = np.zeros((Np,Np),dtype=complex)
    for i,Gi in enumerate(G_pw):
        H[i,i] = 0.5*Gi**2  # kinetic energy
        for j,Gj in enumerate(G_pw):
            dG = Gi-Gj
            idx = np.where(np.isclose(G,dG))[0]
            if idx.size>0:
                H[i,j] += Veff_G[idx[0]]  # potential coupling
    return H

# ============================================================
# 6. Density from Orbitals
# ============================================================

def density_from_orbitals(C, occ):
    """Compute real-space density from KS orbitals
    
    Steps:
    1) Embed plane-wave coefficients into full FFT grid
    2) IFFT to get real-space orbitals
    3) Sum |psi|^2 weighted by occupations
    """
    n = np.zeros(Nx)
    for i,f in enumerate(occ):
        psiG = np.zeros(Nx,dtype=complex)
        psiG[pw_mask] = C[:,i]
        psi = ifft(psiG)
        n += f*np.abs(psi)**2
    return n

# ============================================================
# 7. Ground-State SCF
# ============================================================

occ = [2.0]*(Ne//2)        # closed-shell occupation
n = np.ones(Nx)*(Ne/L)     # initial uniform guess

for it in range(max_scf_iter):
    # Compute effective potential
    VH = hartree_potential(n)    # Hartree potential
    Vxc = v_xc(n)                # LDA exchange potential
    Veff = V_ext + VH + Vxc      # total effective potential
    
    # Solve KS eigenproblem H ψ = ε ψ
    H = ks_hamiltonian(Veff)
    eigvals, eigvecs = eigh(H)
    
    # Update density and mix for SCF stability
    n_new = density_from_orbitals(eigvecs, occ)
    n = mixing*n_new + (1-mixing)*n

print("Ground-state converged.")
print("KS eigenvalues (eV):", eigvals[:10]*Hartree2eV)

# Determine HOMO-LUMO gap
KS_HOMO = eigvals[0]
KS_HOMO_idx = 0
for i, ev in enumerate(eigvals):
    if ev > 0:
        break
    if ev > KS_HOMO:
        KS_HOMO = ev
        KS_HOMO_idx = i

print("KS HOMO eigenvalue (eV):", KS_HOMO*Hartree2eV)
print("KS HOMO-LUMO+{i} gap",(eigvals[KS_HOMO_idx+1:10] - KS_HOMO)*Hartree2eV)

# ============================================================
# 8. Build LR-TDDFT Matrices
# ============================================================

Nocc = Ne//2
Nvirt = 6
Ntrans = Nocc*Nvirt  # total number of particle-hole transitions

# ---- TDA matrix ----
A_TDA = np.zeros((Ntrans,Ntrans))

# ---- Full TDDFT matrices ----
A_full = np.zeros((Ntrans,Ntrans))
B_full = np.zeros((Ntrans,Ntrans))

# Precompute Hartree kernel in real space
fH = np.zeros(Nx)
for i,g in enumerate(G):
    if abs(g)>1e-12:
        fH[i] = 4*np.pi/(g**2)
fH = np.real(ifft(fH))

# XC kernel
fXC = f_xc_kernel(n)
fHxc = fH + fXC  # total TDDFT response kernel

# Build list of (i -> a) transitions
transitions = []
for i in range(Nocc):
    for a in range(Nocc,Nocc+Nvirt):
        transitions.append((i,a))

# Precompute real-space orbitals
psi_real = []
for n in range(Nocc + Nvirt):
    psiG = np.zeros(Nx, dtype=complex)
    psiG[pw_mask] = eigvecs[:, n]
    psi_real.append(ifft(psiG))

# Fill TDA, A_full, B_full matrices
for p,(i,a) in enumerate(transitions):
    for q,(j,b) in enumerate(transitions):
        # Real-space orbitals
        psi_i = psi_real[i]
        psi_a = psi_real[a]
        psi_j = psi_real[j]
        psi_b = psi_real[b]
        
        # Transition densities: rho_ia = ψ_i* ψ_a
        rho_ia = np.conjugate(psi_i) * psi_a
        rho_jb = np.conjugate(psi_j) * psi_b
        
        # Kernel contraction: < ia | f_Hxc | jb >
        kernel_term = np.sum(rho_ia * fHxc * np.conjugate(rho_jb)) * dx
        
        # ---- TDA ----
        if p==q:
            # diagonal term: independent particle energy difference
            A_TDA[p,q] += eigvals[a] - eigvals[i]
            # A_full also includes same term for consistency
            A_full[p,q] += eigvals[a] - eigvals[i]
        A_TDA[p,q] += np.real(kernel_term)
        
        # ---- Full TDDFT ----
        A_full[p,q] += np.real(kernel_term)
        B_full[p,q] += np.real(kernel_term)  # B = <ia|f_Hxc|jb>

# ============================================================
# 9. Solve TDDFT Eigenvalue Problems
# ============================================================

# ---- TDA ----
omega_TDA, X_TDA = eigh(A_TDA)

# ---- Full TDDFT (Casida equations) ----
# Solve (A-B)(A+B) X = omega^2 X
C = (A_full - B_full) @ (A_full + B_full)
eigvals2, X_full = eigh(C)
omega_full = np.sqrt(np.maximum(eigvals2, 0.0))  # only physical positive excitations

# ============================================================
# 10. Plot Spectrum (Improved)
# ============================================================

# Take first 5 excitations
DFT_gaps_eV = np.array([eigvals[KS_HOMO_idx+1+i] - KS_HOMO for i in range(5)]) * Hartree2eV
TDA_eV = omega_TDA[:5] * Hartree2eV
Full_eV = omega_full[:5] * Hartree2eV

plt.figure(figsize=(8,5))

# DFT gaps: black circles
plt.stem(DFT_gaps_eV, np.ones_like(DFT_gaps_eV), linefmt='k-', markerfmt='ko', basefmt=" ", label="DFT Gap")
# TDA TDDFT: blue triangles
plt.stem(TDA_eV, np.ones_like(TDA_eV)*1.1, linefmt='b-', markerfmt='b^', basefmt=" ", label="TDA TDDFT")
# Full TDDFT: red squares
plt.stem(Full_eV, np.ones_like(Full_eV)*1.2, linefmt='r-', markerfmt='rs', basefmt=" ", label="Full TDDFT")

plt.xlabel("Excitation Energy (eV)")
plt.ylabel("Intensity (arb.)")
plt.title("Excitation Energies Comparison: DFT vs TDA vs Full TDDFT")
plt.legend()
plt.ylim(0,1.5)
plt.show()

# ============================================================
# 11. Print Table of Excitation Energies
# ============================================================

print("\nExcitation Energies (eV):")
print(f"{'Index':>5} {'DFT Gap':>12} {'TDA TDDFT':>12} {'Full TDDFT':>12}")
for idx in range(min(5, Ntrans)):
    dft_gap = (eigvals[KS_HOMO_idx+1+idx]-KS_HOMO)*Hartree2eV if KS_HOMO_idx+1+idx < len(eigvals) else 0.0
    print(f"{idx+1:5d} {dft_gap:12.4f} {omega_TDA[idx]*Hartree2eV:12.4f} {omega_full[idx]*Hartree2eV:12.4f}")
