# constants.py
import numpy as np

initial = [-.0008, .0006, .00085, .09, -.07, .15]  # rad/sec, rad
R_E = 6378*1000  # m, radius of earth
I = [4500, 6000, 7000]  # kg-m^2
a = 7000*1000  # m
G = 6.674e-11 # N m^2/kg^2
M = 5.9722e24  # kg
D = 3  # Ampere-m^2
T = 2*np.pi*np.sqrt(a**3/G/M)  # sec
T_hours = T/3600
i = np.deg2rad(83)  # inclination
n = 2*np.pi/T  # mean motion of orbiting satellite
lat = np.deg2rad(0)  # latitude
long = np.deg2rad(0)  # longitude
omega_0 = .001  # rad/sec, for the Earth's rotation
long_ascnode = lambda t: long - omega_0*t

# standard rotation matrices for any angle
R_3 = lambda angle: np.array([[np.cos(angle),np.sin(angle),0],
              [-np.sin(angle),np.cos(angle),0],
              [0,0,1]])

R_1 = lambda angle: np.array([[1,0,0],
                 [0,np.cos(angle),np.sin(angle)],
                 [0,-np.sin(angle),np.cos(angle)]])

R_2 = lambda angle: np.array([[np.cos(angle), 0, -np.sin(angle)],
                 [0,1,0],
                 [np.sin(angle), 0, np.cos(angle)]])

# R_BO = lambda t: R_3(n*t)@R_1(i)@R_3(long_ascnode(t))

def R_BO(theta,psi,phi):  # rotation from orbit to body frame, 2-1-3
    r_bo = R_3(theta) @ R_1(psi) @ R_2(phi)
    return r_bo

# Earth's magnetic field constants
D_vec = np.array([D,D,D])  # Am^2
B_0 = 3.11e-5  # Tesla
mhat_E = np.array([[-0.05585],[.16986],[-0.98389]])*7.77e22
mag_mom = 7.95e15

# Solar panel constants
plate_areas = np.array([5,7,7])  # m^2
beta = np.deg2rad(30)
plate_normals = np.array([[0,0,1],
                          [0,np.sin(beta),np.cos(beta)],
                          [0, np.sin(beta), np.cos(beta)]])
plate_cm = np.array([[0.1,0.1,0],
                     [2,0,0],
                     [-2,0,0]])
plate_absorb = np.array([0.1,0.7,0.7])
plate_spec = np.array([0.8,0.2,0.2])
plate_diffus = np.array([0.1,0.1,0.1])