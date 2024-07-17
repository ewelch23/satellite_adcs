# rollpitchyaw.py

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from numpy import cos, sin, tan
import numpy as np

# initial = [omega1,omega2,omega3,psi,phi,theta]
initial = [-.0008, .0006, .00085, .09, -.07, .15]  # rad/sec, rad
I = [4500, 6000, 7000]  # kg-m^2
r = 7000*1000  # m
G = 6.674e-11 # N m^2/kg^2
M = 5.9722e24  # kg
T = 2*np.pi*np.sqrt(r**3/G/M)  # sec

def sys(t, y, i=I, T=T):
    omega0 = 2*np.pi/T
    omega1,omega2,omega3,psi,phi,theta = y
    I1,I2,I3 = i
    omega1_dot = (I2-I3)*omega2*omega3/I1
    omega2_dot = (I3-I1)*omega1*omega3/I2
    omega3_dot = (I1-I2)*omega1*omega2/I3
    phi_dot = ((omega1*sin(theta)+omega2*cos(theta))-omega0*sin(psi)*cos(phi))/cos(psi)
    psi_dot = ((omega1*cos(theta)-omega2*sin(theta))+omega0*sin(phi))
    theta_dot = omega3 + (omega1*sin(theta)+omega2*cos(theta))*tan(psi)-omega0*cos(phi)/cos(psi)
    return [omega1_dot,omega2_dot,omega3_dot,psi_dot,phi_dot,theta_dot]

t = np.linspace(0,3*T,200)
T_hours = T/3600
sol = solve_ivp(sys, [0,3*T], initial, dense_output=True)
omega1,omega2,omega3,yaw,roll,pitch = sol.sol(t)
plt.figure(0)
plt.plot(t, yaw.T, color="darkblue")
plt.xlabel("T, Seconds")
plt.ylabel(r"$\psi$, radians")
plt.title("Yaw Time History")
plt.grid()
plt.savefig("yaw_nodisturb.png")

plt.figure(1)
plt.plot(t, roll.T, color="purple")
plt.xlabel("T, Seconds")
plt.ylabel(r"$\phi$, radians")
plt.title("Roll Time History")
plt.grid()
plt.savefig("roll_nodisturb.png")

plt.figure(2)
plt.plot(t, pitch.T, color="darkgreen")
plt.xlabel("T, Seconds")
plt.ylabel(r"$\theta$, radians")
plt.title("Pitch Time History")
plt.grid()
plt.savefig("pitch_nodisturb.png")