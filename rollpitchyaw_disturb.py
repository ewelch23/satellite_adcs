# rollpitchyaw_disturb.py

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from numpy import cos, sin, tan
import numpy as np
from constants import n, T, initial, I

###################################################
# Earth's magnetic field
from mag_field_torque import tau_mag
taux_mag = lambda t: tau_mag(t)[0]
tauy_mag = lambda t: tau_mag(t)[1]
tauz_mag = lambda t: tau_mag(t)[2]
###################################################
# Grav gradient
from grav_grad import M_g
def Mg_x(theta,psi,phi):
    Mgx = M_g(theta,psi,phi)[0]
    return Mgx
def Mg_y(theta,psi,phi):
    Mgy = M_g(theta,psi,phi)[1]
    return Mgy
def Mg_z(theta,psi,phi):
    Mgz = M_g(theta,psi,phi)[2]
    return Mgz
###################################################
# Solar radiation
from solar_torque import tau_sp
def tau_spx(theta,psi,phi,t):
    tauspx = tau_sp(theta,psi,phi,t)[0,2]
    return tauspx
def tau_spy(theta,psi,phi,t):
    tauspy = tau_sp(theta,psi,phi,t)[1,2]
    return tauspy
def tau_spz(theta,psi,phi,t):
    tauspz = tau_sp(theta,psi,phi,t)[2,2]
    return tauspz
###################################################

def sys(t, y):
    omega1,omega2,omega3,psi,phi,theta = y
    I1,I2,I3 = I
    omega1_dot = (I2-I3)*omega2*omega3/I1+(taux_mag(t)+Mg_x(theta,psi,phi)+tau_spx(theta,psi,phi,t))/I1
    omega2_dot = (I3-I1)*omega1*omega3/I2+(tauy_mag(t)+Mg_y(theta,psi,phi)+tau_spy(theta,psi,phi,t))/I2
    omega3_dot = (I1-I2)*omega1*omega2/I3+(tauz_mag(t)+Mg_z(theta,psi,phi)+tau_spz(theta,psi,phi,t))/I3
    phi_dot = ((omega1*sin(theta)+omega2*cos(theta))-n*sin(psi)*cos(phi))/cos(psi)
    psi_dot = ((omega1*cos(theta)-omega2*sin(theta))+n*sin(phi))
    theta_dot = omega3 + (omega1*sin(theta)+omega2*cos(theta))*tan(psi)-n*cos(phi)/cos(psi)
    return [omega1_dot,omega2_dot,omega3_dot,psi_dot,phi_dot,theta_dot]

orbit_period = np.linspace(0,3*T,200)
sol = solve_ivp(sys, (0,3*T), initial, dense_output=True)
omega_1,omega_2,omega_3,yaw,roll,pitch = sol.sol(orbit_period)

plt.figure(0)
plt.plot(orbit_period/T, yaw.T, color="darkblue")
plt.xlabel("T, Orbits")
plt.ylabel(r"$\psi$, radians")
plt.title("Yaw Time History")
plt.grid()
plt.savefig("yaw_disturb.png")
plt.close()

plt.figure(1)
plt.plot(orbit_period/T, roll.T, color="purple")
plt.xlabel("T, Orbits")
plt.ylabel(r"$\phi$, radians")
plt.title("Roll Time History")
plt.grid()
plt.savefig("roll_disturb.png")
plt.close()

plt.figure(2)
plt.plot(orbit_period/T, pitch.T, color="darkgreen")
plt.xlabel("T, Orbits")
plt.ylabel(r"$\theta$, radians")
plt.title("Pitch Time History")
plt.grid()
plt.savefig("pitch_disturb.png")
plt.close()