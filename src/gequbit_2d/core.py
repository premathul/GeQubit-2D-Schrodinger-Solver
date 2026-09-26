import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

HBAR_J_S = 1.054571817e-34
M0_KG = 9.1093837139e-31
J_PER_EV = 1.602176634e-19

def hamiltonian_2d(x_m, y_m, potential_ev, effective_mass_m0):
    x=np.asarray(x_m,float); y=np.asarray(y_m,float); v=np.asarray(potential_ev,float)
    if v.shape != (y.size, x.size):
        raise ValueError("potential shape must be (len(y), len(x))")
    if x.size<3 or y.size<3 or effective_mass_m0<=0:
        raise ValueError("invalid grid or mass")
    dx=np.diff(x); dy=np.diff(y)
    if not np.allclose(dx,dx[0]) or not np.allclose(dy,dy[0]):
        raise ValueError("uniform grids required")
    tx=HBAR_J_S**2/(2*effective_mass_m0*M0_KG*dx[0]**2)/J_PER_EV
    ty=HBAR_J_S**2/(2*effective_mass_m0*M0_KG*dy[0]**2)/J_PER_EV
    nx=x.size; ny=y.size
    lx=sparse.diags([-tx*np.ones(nx-1),2*tx*np.ones(nx),-tx*np.ones(nx-1)],[-1,0,1])
    ly=sparse.diags([-ty*np.ones(ny-1),2*ty*np.ones(ny),-ty*np.ones(ny-1)],[-1,0,1])
    kinetic=sparse.kron(sparse.eye(ny),lx)+sparse.kron(ly,sparse.eye(nx))
    return kinetic+sparse.diags(v.ravel())

def solve_lowest_states(x_m, y_m, potential_ev, effective_mass_m0, n_states=4):
    H=hamiltonian_2d(x_m,y_m,potential_ev,effective_mass_m0)
    vals,vecs=eigsh(H,k=n_states,which="SA")
    order=np.argsort(vals)
    vals=vals[order]; vecs=vecs[:,order]
    dx=float(np.diff(np.asarray(x_m))[0]); dy=float(np.diff(np.asarray(y_m))[0])
    psi=vecs.reshape(len(y_m),len(x_m),n_states)
    norms=np.sqrt(np.sum(np.abs(psi)**2,axis=(0,1))*dx*dy)
    psi/=norms
    return vals,psi
