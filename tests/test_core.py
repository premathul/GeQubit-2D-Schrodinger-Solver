import numpy as np
from gequbit_2d.core import solve_lowest_states
from gequbit_2d.potentials import harmonic_potential_2d_ev

def test_ground_state_normalization_and_ordering():
    x=np.linspace(-40e-9,40e-9,21)
    y=np.linspace(-40e-9,40e-9,21)
    m=0.08
    V=harmonic_potential_2d_ev(x,y,m,2*np.pi*80e9,2*np.pi*80e9)
    E,psi=solve_lowest_states(x,y,V,m,3)
    assert np.all(np.diff(E)>=0)
    dx=x[1]-x[0]; dy=y[1]-y[0]
    assert np.isclose(np.sum(abs(psi[:,:,0])**2)*dx*dy,1.0,rtol=1e-6)
