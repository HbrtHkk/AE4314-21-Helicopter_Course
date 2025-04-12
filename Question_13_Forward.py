from Question_12_13_Hover import induced_hover
from Question_12_Forward import vi_bar
from Question_13_BEM_Hover import CL, CLbar0 
import numpy as np 
import matplotlib.pyplot as plt

#------------------------- QUESTION 1.3 FORWARD FLIGHT  -------------------------

def Pind(k, T, vi):
    return k*T*vi   

# Profile power and drag 
def Pprofile(sigma, cdp, rho, R, omega, mu):
    Pphov = ((sigma*cdp)/8)*rho*(omega*R)**3*(np.pi*R**2)
    P = Pphov*(1+4.65*mu**2)
    return P

# Parasite drag
def Pparasite(rho, V, CdS):
    return 0.5*rho*V**3*CdS

def tailrotor_power(k, T, vit, rho0, Rt, sigmat, Cdp, omegat, mut):
    pindt= 1.1*k*T*vit
    ppdt = (sigmat*Cdp*rho0*np.pi*(Rt**2)*((omegat*Rt)**3)*(1+4.65*mut**2))/8
    return pindt + ppdt

def tailCL(T, sigma, V, rho, S):
    Ct = T/0.5*rho*V**2*S
    CLbar = sigma*Ct/6.6
    return CLbar

#------------------------- STANDARD PARAMETERS -------------------------
Vcruise = 123*0.5144 #m/s cruise
Vmax = 170*0.5144 #m/s max

T0 = 288.15  # Temperature (K)
rho0 = 1.225  # Density at sea level (kg/m^3)
g = 9.81 # m/s^2 gravitational acceleration
mass = 4536 # kg mass of the helicopter MTOW
W = mass * g # N helicopter weight
diam = 13.41 # m Rotor diameter
R = diam/2 # m Rotor radius 
c = 0.76 # m chord length
rotor_RPM = 324 
rotor_speed = rotor_RPM * 2 * np.pi/60  # rad/s angular rotor speed
N = 2 # Number of blades
cR = c/R

#------------------------- Defining Characteristics -------------------------
print(" ")
print("Defining Characteristics needed for graphical values")
print("The blade tip speed (Mach) is: ", np.round(rotor_speed*R/340, 2), "[-]")
print(r'The \bar{C_{L}} is :', np.round(CLbar0, 2), "[-]")
print(" ")
#----------------------------------------------------------------------------
# Velocity range 
V = np.linspace(0, Vmax, len(vi_bar))

P_tot = []
P_ind = []
P_profile = []
P_parasite = []
P_t = []

for i in np.arange(0, len(V)):
    # MAIN ROTOR PARAMETERS
    lt = 7.901 # [m] Distance between rotor and tail rotor
    
    viff = vi_bar[i][0]*induced_hover(W, rho0, R)
    k = 1.15
    CdS = 2.415 # [m^2] Graphically determined 
    CDp = 0.01 # [-] Graphically determined 
    sigma = N*c/(np.pi*R)
    mu = V[i]/(rotor_speed*R)

    #### TAIL ROTOR PARAMETERS
    D_tail = 2.59 # m tail rotor diameter
    Rt = D_tail/2 # m tail rotor radius
    ct = 0.30 # m tail rotor chord
    kt = 1.4
    Omega_tail = 1654 # RPM tail rotor
    rotor_speed_t = Omega_tail * 2 * np.pi / 60
    sigmat = N*ct/(np.pi*Rt)
    mut = V[i]/(rotor_speed_t*Rt) # Advance ratio 

    # Rotor power
    P_rotor = Pind(k, W, viff) + Pprofile(sigma, CDp, rho0, R, rotor_speed, mu)+ Pparasite(rho0, V[i], CdS)

    T = P_rotor/(rotor_speed*lt)
    vit = induced_hover(T, rho0, Rt)

    # Tail rotor power
    P_tail = tailrotor_power(kt, T, vit, rho0, Rt, sigmat, CDp, rotor_speed_t, mut)

    # Total power
    P_tot.append((P_rotor+P_tail)/10**6) 
    
    # Power components
    P_ind.append(Pind(k, W, viff)/10**6)
    P_profile.append(Pprofile(sigma, CDp, rho0, R, rotor_speed, mu)/10**6)
    P_parasite.append(Pparasite(rho0, V[i], CdS)/10**6)
    P_t.append(P_tail/10**6)

    # Percentage of tail power to overall power
    print("% of tail power to overall power:", np.round(P_tail*100/(P_rotor+P_tail), 2), "%")

# Determining the point of that gives the best range
dP_dV = np.gradient(P_tot, V)

P_over_V = P_tot / V  

# Finding the point where P/V ≈ dP/dV
best_idx = np.argmin(np.abs(P_over_V - dP_dV))

# Best range velocity and power
V_best_range = V[best_idx]
P_best_range = P_tot[best_idx]

# Plot the tangent line from (0,0)
slope = P_best_range / V_best_range
line = slope * V  # y = mx


plt.plot(V, P_tot, color = 'black', label = 'Total Power')
plt.plot(V, P_ind, color = '#D81B60', linestyle = '--', label = 'Induced Power')
plt.plot(V, P_profile, color = '#FFB000', linestyle = '--', label = 'Profile Power')
plt.plot(V, P_parasite, color = '#785EF0', linestyle = '--', label = 'Parasite Power')
plt.plot(V, P_t, color = 'purple', linestyle = '--', label = 'Tail Rotor Power')
plt.legend()
plt.xlabel('Forward Velocity [m/s]')
plt.ylabel('Power [MW]')
plt.grid()
plt.show()

fig, axs = plt.subplots(1, 2, figsize = (16, 8))
axs[0].plot(V / 0.51444, P_tot, color = 'black', label = 'Total Power')
axs[0].plot((V / 0.51444), P_ind, color = '#D81B60', linestyle = '--', label = 'Induced Power')
axs[0].plot((V / 0.51444), P_profile, color = '#FFB000', linestyle = '--', label = 'Profile Power') 
axs[0].plot((V / 0.51444), P_parasite, color = '#785EF0', linestyle = '--', label = 'Parasite Power')
axs[0].plot((V / 0.51444), P_t, color = 'purple', linestyle = '--', label = 'Tail Rotor Power')
axs[0].set_xlabel('Forward Velocity [knots]')
axs[0].set_ylabel('Power [MW]')
#axs[0].set_title('Power vs Forward Velocity (BEM)')
axs[0].legend()
axs[0].grid()

axs[1].plot(V / 0.51444, line, linestyle="dashed", color="orange")
axs[1].scatter(V_best_range / 0.51444, P_best_range, color='orange', label="Point B", zorder=3)
axs[1].scatter(V[P_tot.index(min(P_tot))] / 0.51444, min(P_tot), color = '#D81B60', label = 'Point A', zorder=3)
axs[1].plot(V / 0.51444, P_tot, color = 'black', label = 'Total Power')
axs[1].set_ylim(0, max(P_tot))
axs[1].set_xlabel('Forward Velocity [knots]')
axs[1].set_ylabel('Power [MW]')
#axs[1].set_title('Power vs Forward Velocity (BEM)')
axs[1].legend()
axs[1].grid()
plt.show()

print("best range velocity is: ", np.round(V_best_range/0.51444, 2), "kts")
print("best range power is: ", np.round(P_best_range, 2), "MW")
print("best endurance velocity is: ", np.round(V[np.argmin(P_tot)]/0.51444, 2), "kts")
print("best endurance power is: ", np.round(np.min(P_tot), 2), "MW")
