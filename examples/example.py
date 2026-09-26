import numpy as np
from gequbit_2d.core import solve_lowest_states
from gequbit_2d.potentials import harmonic_potential_2d_ev

x=np.linspace(-60e-9,60e-9,41)
y=np.linspace(-60e-9,60e-9,41)
m=0.08
omega=2*np.pi*100e9
V=harmonic_potential_2d_ev(x,y,m,omega,omega)
E,_=solve_lowest_states(x,y,V,m,4)
print("Lowest energies [meV]:", E*1e3)
