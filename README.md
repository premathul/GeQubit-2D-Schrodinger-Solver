# GeQubit-2D-Schrodinger-Solver

GeQubit-2D-Schrodinger-Solver is a numerical research project for solving two-dimensional effective-mass Schrödinger problems relevant to gate-defined quantum dots in Ge/SiGe heterostructures. The goal is to provide a transparent intermediate layer between analytic confinement models and full commercial device solvers. Rather than beginning with a highly specialized multiband Hamiltonian, the project first establishes a reliable numerical foundation for sparse finite-difference eigenproblems, potential landscapes, orbital energies, wavefunctions, expectation values, and mesh convergence.

The starting equation is
[
left[
-rac{hbar^2}{2m^*}
left(
rac{partial^2}{partial x^2}+
rac{partial^2}{partial y^2}
ight)+V(x,y)
ight]psi_n(x,y)=E_npsi_n(x,y).
]
The current implementation assumes a constant scalar effective mass and a rectangular uniform grid with Dirichlet boundaries. These assumptions make the solver easy to validate and provide a controlled baseline for later extensions to anisotropic masses, spatially varying material parameters, coupled bands, magnetic vector potentials, and realistic electrostatic potentials imported from device simulations.

The project uses sparse matrices through SciPy so that grids substantially larger than simple dense-matrix examples can be treated efficiently. The Hamiltonian is assembled from Kronecker-sum Laplacians, and the lowest eigenstates are extracted using sparse eigensolvers. Utility functions provide common synthetic potentials such as two-dimensional harmonic confinement and double-well structures, which are useful for benchmarking orbital spacing and tunnel-coupled states before realistic device potentials are introduced.

Installation is performed with
```bash
git clone https://github.com/premathul/GeQubit-2D-Schrodinger-Solver.git
cd GeQubit-2D-Schrodinger-Solver
python -m pip install -e .
```
and tests can be run using
```bash
python -m pip install -e .[dev]
pytest -q
```.

The long-term purpose of this repository is to support a reproducible chain from gate-defined electrostatic potential to orbital states and eventually to quantities needed by spin-qubit models, including dipole matrix elements, orbital splittings, tunnel coupling, electric-field expectation values, and effective spin parameters. The present implementation should therefore be viewed as the numerical core of a larger Ge qubit device-simulation stack.

## Runnable scientific baseline

The current implementation discretizes a scalar effective-mass Hamiltonian on the interior nodes of a square grid with zero wavefunction at the boundary. Its potential is an isotropic parabolic confinement expressed in meV, with coordinates in nm. The kinetic coefficient converts SI constants to meV, while each returned eigenvector is normalized so the two-dimensional numerical integral of |ψ|² equals one. Symmetry can make excited levels degenerate, so their individual eigenvectors may rotate within a degenerate subspace even when the energies converge.

Install `numpy` and `scipy`, then run `python src/main.py --n 51 --width-nm 160 --mass-ratio 0.08`. Repeat at larger grid sizes and wider boxes to distinguish discretization error from boundary confinement. This is a single-band orbital teaching model: a quantitative Ge hole-spin calculation requires valence-band structure, electrostatics, strain, heterointerfaces, spin-orbit coupling, and validated boundary conditions.

## Validation and scope

The calculations in `src/main.py` are transparent baseline models intended for reproducibility and extension. Inputs and assumptions should be reported alongside outputs; numerical agreement with a plotted trace alone does not validate a material-specific prediction. New physical terms should be accompanied by dimensional checks and independent limiting-case comparisons.

## Contact

**Athul Prem** — [GitHub profile](https://github.com/premathul). For scientific discussion or collaboration, open an issue in this repository or reach out through my GitHub profile.
