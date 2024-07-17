#plot_eulers.py
# plots the Euler angles of the system with respect to time
# for each instance of disturbance torques active
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from numpy import sin,cos,tan
from constants import T,n,omega_0,I
from mag_field_torque import tau_mag
from solar_torque import tau_sp
from grav_grad import M_g

case_1_initial = [0,0,omega_0,0,0,0]  # omega1,omega2,omega3,psi,phi,theta
case_2_initial = [0,0,omega_0,0.1,0.15,0.2]  # omega1,omega2,omega3,psi,phi,theta

def sys(t,y,taux_mag,tauy_mag,tauz_mag,tau_spx,tau_spy,tau_spz,Mg_x,Mg_y,Mg_z):
    '''The dynamics system as stated in the slides, set up for solve_ivp.'''
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
torques = ["electro_magnetic","gravity_gradient","solar_pressure","none"]

def plot_euler(sys,torques,initial,case):
    '''Plots the euler angles for each instance of distubance torque active,
    using any system, any initial values, and both of the two cases given in the assignment.'''
    for torque in torques:
        fig, axs = plt.subplots(3, 1, figsize=(10, 18))
        fig.suptitle(f'Euler Angles vs Time for Torque: {torque}', fontsize=16)
        euler_labels = [r'$\psi$', r'$\phi$', r'$\theta$']
        if torque=="electro_magnetic":
            taux_mag = lambda t: tau_mag(t)[0]
            tauy_mag = lambda t: tau_mag(t)[1]
            tauz_mag = lambda t: tau_mag(t)[2]
            Mg_x = lambda theta,psi,phi: 0
            Mg_y = lambda theta,psi,phi: 0
            Mg_z = lambda theta,psi,phi: 0
            tau_spx = lambda theta,psi,phi,t: 0
            tau_spy = lambda theta,psi,phi,t: 0
            tau_spz = lambda theta,psi,phi,t: 0
        elif torque=="gravity_gradient":
            taux_mag = lambda t: 0
            tauy_mag = lambda t: 0
            tauz_mag = lambda t: 0
            Mg_x = lambda theta,psi,phi: M_g(theta,psi,phi)[0]
            Mg_y = lambda theta,psi,phi: M_g(theta,psi,phi)[1]
            Mg_z = lambda theta,psi,phi: M_g(theta,psi,phi)[2]
            tau_spx = lambda theta,psi,phi,t: 0
            tau_spy = lambda theta,psi,phi,t: 0
            tau_spz = lambda theta,psi,phi,t: 0
        elif torque=="solar_pressure":
            taux_mag = lambda t: 0
            tauy_mag = lambda t: 0
            tauz_mag = lambda t: 0
            Mg_x = lambda theta,psi,phi: 0
            Mg_y = lambda theta,psi,phi: 0
            Mg_z = lambda theta,psi,phi: 0
            tau_spx = lambda theta,psi,phi,t: tau_sp(theta,psi,phi,t)[0,2]
            tau_spy = lambda theta,psi,phi,t: tau_sp(theta,psi,phi,t)[1,2]
            tau_spz = lambda theta,psi,phi,t: tau_sp(theta,psi,phi,t)[1,2]
        else:  # no disturbance torques
            taux_mag = lambda t: 0
            tauy_mag = lambda t: 0
            tauz_mag = lambda t: 0
            Mg_x = lambda theta,psi,phi: 0
            Mg_y = lambda theta,psi,phi: 0
            Mg_z = lambda theta,psi,phi: 0
            tau_spx = lambda theta,psi,phi,t: 0
            tau_spy = lambda theta,psi,phi,t: 0
            tau_spz = lambda theta,psi,phi,t: 0
    
        sol = solve_ivp(sys,(0,3*T),initial,args=(taux_mag,tauy_mag,tauz_mag,tau_spx,tau_spy,tau_spz,Mg_x,Mg_y,Mg_z),dense_output=True)
        omega_1,omega_2,omega_3,psi,phi,theta = sol.sol(orbit_period)
        euler_values = [psi, phi, theta]
        color_values = ["darkblue","purple","darkgreen"]

        for i, ax in enumerate(axs):
            ax.plot(orbit_period / T, euler_values[i], label=euler_labels[i],color=color_values[i])
            ax.set_xlabel('Time [orbits]')
            ax.set_ylabel('Angle [rad]')
            ax.legend()
            ax.grid(True)

        plt.tight_layout()
        plt.savefig(str(torque)+"_"+(case)+".png")
        plt.close()

cases = ["case1","case2"]

for case in cases:
    if case == 'case1':
        initials = case_1_initial
    elif case == 'case2':
        initials = case_2_initial
    else:
        initials = [0,0,0,0,0,0]
    plot_euler(sys,torques,initials,case)
