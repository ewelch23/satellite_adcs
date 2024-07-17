import matplotlib.pyplot as plt
import numpy as np
from constants import T,n
from rollpitchyaw_disturb import yaw, roll, pitch
from mag_field_torque import tau_mag
from solar_torque import tau_sp
from grav_grad import M_g

psi, phi, theta = yaw, roll, pitch
orbit_period = np.linspace(0, 3*T, 200)
time_values = orbit_period/T

# Calculate torques for each type
torque_mag = np.array([tau_mag(t) for t in orbit_period])
Mg = np.array([M_g(theta[i], psi[i], phi[i]) for i in range(len(orbit_period))])

# Calculate solar torques with corrected time values
tau_sp_data = np.array([tau_sp(theta[i], psi[i], phi[i], orbit_period[i] * T) for i in range(len(orbit_period))])
tau_spx = np.zeros(len(orbit_period))
tau_spy = np.zeros(len(orbit_period))
tau_spz = np.zeros(len(orbit_period))

for i, submatrix in enumerate(tau_sp_data):
    non_zero_values = [value for row in submatrix for value in row if value != 0]
    if len(non_zero_values) >= 3:
        tau_spx[i], tau_spy[i], tau_spz[i] = non_zero_values[:3]
    else:
        tau_spx[i], tau_spy[i], tau_spz[i] = [0] * 3

# Plotting magnetic torques
plt.figure(figsize=(10, 6))
plt.plot(time_values, torque_mag[:,0], label='Magnetic X',color="darkgreen")
plt.plot(time_values, torque_mag[:,1], label='Magnetic Y',color="purple")
plt.plot(time_values, torque_mag[:,2], label='Magnetic Z',color="darkblue")
plt.xlabel('Time (Orbits)')
plt.ylabel('Torque')
plt.title('Magnetic Torques')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("em_torques.png")
plt.close()

# Plotting solar torques, approximated when theta=phi=0, psi=30deg
A_x = -1.3e-6
B_x = 0.22e-6
theta_approx = np.deg2rad(-110)
tau_spx_approx = lambda t: A_x + B_x*np.cos(n*t+theta_approx)
A_y = 1.3e-6
B_y = B_x
tau_spy_approx = lambda t: A_y + B_y*np.cos(n*t+np.deg2rad(70))
A_z = 0.16e-6
B_z = 0.15e-6
tau_spz_approx = lambda t: A_z + B_z*np.cos(n*t+np.deg2rad(-40))

plt.figure(figsize=(10, 6))
plt.plot(time_values, tau_spx_approx(orbit_period), label=r'$\tau_{x}$, approx',color="darkgreen")
# plt.plot(time_values, tau_spx, label=r'$\tau_x$')
plt.plot(time_values, tau_spy_approx(orbit_period), label=r'$\tau_{y}$, approx',color="purple")
# plt.plot(time_values, tau_spy, label=r'$\tau_y$')
plt.plot(time_values, tau_spz_approx(orbit_period), label=r'$\tau_z$, approx',color="darkblue")
# plt.plot(time_values, tau_spz, label=r'$\tau_z$')
plt.xlabel('Time (Orbits)')
plt.ylabel('Torque')
plt.title('Solar Torques')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("solar_torques.png")
plt.close()

# Plotting gravity gradient torques
plt.figure(figsize=(10, 6))
plt.plot(time_values, Mg[:, 0], label='Gravity Gradient X',color="darkgreen")
plt.plot(time_values, Mg[:, 1], label='Gravity Gradient Y',color="purple")
plt.plot(time_values, Mg[:, 2], label='Gravity Gradient Z',color="darkblue")
plt.xlabel('Time (orbits)')
plt.ylabel('Torque')
plt.title('Gravity Gradient Torques')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("gg_torques.png")
plt.close()

