"""Sparse 2D effective-mass finite-difference eigenproblem with Dirichlet walls."""
import argparse
import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import eigsh

HBAR = 1.054571817e-34
M_E = 9.1093837139e-31
EV = 1.602176634e-19

def solve(n=51, width_nm=160., mass_ratio=0.08, curvature_mev_nm2=0.001, states=4):
    if n < 5 or width_nm <= 0 or mass_ratio <= 0 or states < 1 or states >= (n-2)**2:
        raise ValueError("Invalid grid, mass, or state count")
    dx_nm = width_nm / (n-1)
    x = np.linspace(-width_nm/2, width_nm/2, n)[1:-1]
    a = diags([-np.ones(n-3), 2*np.ones(n-2), -np.ones(n-3)], [-1,0,1], format="csr")
    lap = kron(a, eye(n-2, format="csr")) + kron(eye(n-2, format="csr"), a)
    xx, yy = np.meshgrid(x, x, indexing="ij")
    potential = curvature_mev_nm2 * (xx**2 + yy**2)
    kinetic_mev = HBAR**2/(2*mass_ratio*M_E*(dx_nm*1e-9)**2)/(EV*1e-3)
    h = kinetic_mev*lap + diags(potential.ravel(), format="csr")
    energies, vectors = eigsh(h, k=states, which="SA", tol=1e-10)
    order = np.argsort(energies)
    vectors = vectors[:,order] / (dx_nm*1e-9)  # integral |psi|^2 dx dy = 1
    return energies[order], vectors, dx_nm

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, default=51)
    p.add_argument("--width-nm", type=float, default=160.)
    p.add_argument("--mass-ratio", type=float, default=0.08)
    args = p.parse_args()
    energies, _, dx = solve(n=args.n, width_nm=args.width_nm, mass_ratio=args.mass_ratio)
    print(f"dx_nm={dx:.6g}; energies_meV={energies}; gap01_meV={energies[1]-energies[0]:.6g}")
if __name__ == "__main__":
    main()
