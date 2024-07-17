# solar_torque.py
import numpy as np
from constants import plate_spec, plate_diffus, plate_cm, plate_areas
from constants import i, n, R_BO, R_1, R_3

P = 4.644e-6  # N/m^2
rt_asc = np.deg2rad(-70)
orbit_angle = lambda t: n*t

def tau_sp(theta,psi,phi,t):
    s_ecliptic = np.array([np.cos(2*np.pi/365*t),
                                np.sin(2*np.pi/365*t),
                                0])
    s_body = R_BO(theta,psi,phi)@R_3(n*t)@R_1(i)@R_3(rt_asc)@R_1(np.deg2rad(-23.5))@s_ecliptic

    n_hat = np.array([-np.cos(orbit_angle(t)),
                    np.sin(orbit_angle(t)),
                    0])

    F_sr = P*plate_areas*(n_hat*s_body)*((1-plate_spec)*s_body + \
                                (plate_spec+2/3*plate_diffus)*n_hat)

    tausp = np.cross(plate_cm, F_sr)
    return tausp
