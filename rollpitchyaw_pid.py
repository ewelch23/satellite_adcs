# rollpitchyaw_pid.py
import time
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from numpy import cos, sin, tan
import numpy as np
from constants import n, T, initial, I, omega_0

###################################################
# Solar radiation
A_x = -1.3e-6
B_x = 0.22e-6
tau_spx = lambda t: A_x + B_x*np.cos(n*t)
A_y = 1.3e-6
B_y = B_x
tau_spy = lambda t: A_y + B_y*np.cos(n*t)
A_z = 0.16e-6
B_z = 0.15e-6
tau_spz = lambda t: A_z + B_z*np.cos(n*t)

###################################################
orbit_period = np.linspace(0,2*T,200)
I1,I2,I3 = I
zeta = 0.95
k_theta = 0
tau_theta = 0  # these aren't really changing anything
k_psi = 1
tau_psi = 1
k_phi = 1
tau_phi = 1

def sys(t, y):
    omega1,omega2,omega3,psi,phi,theta,control_psi,control_phi,control_theta = y
    control_psi = 0
    control_phi = 0
    control_theta = 0
    control_psi_update = k_psi*(tau_psi*psi+control_psi)  # 1
    control_phi_update =  k_phi*(tau_phi*phi+control_phi) # 2
    control_theta_update =  k_theta*(tau_theta*theta+control_theta) # 3
    
    omega1_update = tau_spx(t)/I1 - (control_psi_update+(I3-I2)*omega2*omega3+(control_theta_update*omega2-control_phi_update*omega3))/I1
    omega2_update = tau_spy(t)/I2 - (control_phi_update+(I1-I3)*omega1*omega3+(control_psi_update*omega3-control_theta_update*omega1))/I2
    omega3_update = tau_spz(t)/I3 - (control_theta_update+(I2-I1)*omega1*omega2+(control_phi_update*omega1-control_psi_update*omega2))/I3

    phi_update = ((omega1_update*sin(theta)+omega2_update*cos(theta))-n*sin(psi)*cos(phi))/cos(psi)
    psi_update = ((omega1_update*cos(theta)-omega2_update*sin(theta))+n*sin(phi))
    theta_update = omega3_update + (omega1_update*sin(theta)+omega2_update*cos(theta))*tan(psi)-n*cos(phi)/cos(psi)

    return omega1_update,omega2_update,omega3_update,psi_update,phi_update,theta_update,control_psi_update,control_phi_update,control_theta_update

cases = ["slew30","detumble"]
start = time.time()
for case in cases:
    if case == "slew30":
        initial = [0,0,omega_0,0,0,0,np.deg2rad(30),0,0] # omega1,omega2,omega3,psi,phi,theta,psi_control,phi_control,theta_control
        
    elif case == "detumble":
        initial = [0,0,omega_0,0.1,0.15,0.2,0,0,0]
        
    else: 
        initial = [0,0,0,0,0,0,0,0,0]
        
    sol = solve_ivp(sys, (0,2*T), initial, dense_output=True)
    omega_1,omega_2,omega_3,yaw,roll,pitch,control_psi,control_phi,control_theta = sol.sol(orbit_period)
    plt.figure(3)
    plt.plot(orbit_period/60, np.rad2deg(yaw.T), color="darkblue")
    plt.xlabel("T, (min)")
    plt.ylabel(r"$\psi$, deg")
    plt.title("Yaw Time History")
    plt.grid()
    plt.savefig("yaw_control_"+case+".png")
    plt.close()

    plt.figure(4)
    plt.plot(orbit_period/60, np.rad2deg(roll.T), color="purple")
    plt.xlabel("T, (min)")
    plt.ylabel(r"$\phi$, deg")
    plt.title("Roll Time History")
    plt.grid()
    plt.savefig("roll_control_"+case+".png")
    plt.close()

    plt.figure(5)
    plt.plot(orbit_period/60, np.rad2deg(pitch.T), color="darkgreen")
    plt.xlabel("T, (min)")
    plt.ylabel(r"$\theta$, deg")
    plt.title("Pitch Time History")
    plt.grid()
    plt.savefig("pitch_control_"+case+".png")
    plt.close()

end = time.time()
total_time = (end-start)
print("Runtime: "+str(total_time)+" seconds.")