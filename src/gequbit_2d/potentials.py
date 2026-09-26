import numpy as np
from .core import M0_KG, J_PER_EV

def harmonic_potential_2d_ev(x_m, y_m, mass_m0, omega_x_rad_s, omega_y_rad_s):
    X,Y=np.meshgrid(np.asarray(x_m,float),np.asarray(y_m,float))
    m=mass_m0*M0_KG
    return 0.5*m*(omega_x_rad_s**2*X**2+omega_y_rad_s**2*Y**2)/J_PER_EV

def double_well_potential_ev(x_m, y_m, separation_m, curvature_ev_m4):
    X,Y=np.meshgrid(np.asarray(x_m,float),np.asarray(y_m,float))
    a=separation_m/2
    return curvature_ev_m4*(X**2-a**2)**2 + curvature_ev_m4*(0.25*Y**4)
