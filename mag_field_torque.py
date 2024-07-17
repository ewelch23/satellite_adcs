# mag_field_torque.py
# EW_AERO424_HW3.py
from numpy import sin, cos
import numpy as np
from constants import i, mag_mom, R_E, n, D_vec


# mhat_orb = R_BO(theta,psi,phi)@mhat_E
# mhatx = mhat_orb[0]
# mhaty = mhat_orb[1]
# mhatz = mhat_orb[2]

# rx = lambda t: cos(n*t)
# ry = lambda t: sin(n*t)*cos(i)
# rz = lambda t: sin(i)*sin(n*t)

dipolex = lambda t: mag_mom/R_E**3*cos(n*t)*sin(i)
dipoley = lambda t: mag_mom/R_E**3*2*sin(n*t)*sin(i)
dipolez = -mag_mom/R_E**3*cos(i)
dipole = lambda t: np.array([dipolex(t),dipoley(t),dipolez])
tau_mag = lambda t: np.cross(D_vec,dipole(t))
