# grav_grad.py
import numpy as np
from constants import I, n, R_BO

def M_g(theta,psi,phi):
    '''Returns the gravity gradient vector w.r.t. the Euler angles.'''
    R_11 = R_BO(theta,psi,phi)[0,0]
    R_21 = R_BO(theta,psi,phi)[1,0]
    R_31 = R_BO(theta,psi,phi)[2,0]
    Mg = 3*n**2* np.array([R_21*R_31*(I[2]-I[1]),
                        R_11*R_31*(I[0]-I[2]),
                        R_11*R_21*(I[1]-I[0])])
    return Mg