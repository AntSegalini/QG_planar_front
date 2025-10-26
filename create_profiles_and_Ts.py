import numpy as np

Lz=1e4
z=np.linspace(0,Lz,100)

U=10*(z/Lz)
dUdz=10/Lz*np.ones_like(z)
d2Udz2=np.zeros_like(z)

N2=1e-4
g=9.80665
rho=np.ones_like(z)
drhodz=np.zeros_like(z)
# rho=np.exp(-N2/g*z);drhodz=-N2/g*rho  # Exponential profile
# T=288.15-0.0065*z
# rho=1.225*(T/288.15)**4.25;drhodz=-0.0065/288.15*4.25*1.225*(T/288.15)**3.25

N2=N2*np.ones_like(z)
dN2dz=np.zeros_like(z)

# Write profiles to file
np.savetxt('PROFILES.txt',np.column_stack((z,U,dUdz,d2Udz2,rho,drhodz,N2,dN2dz)),header='z [m], U [m/s], dUdz [1/s], d2Udz2 [1/m*s], rho [kg/m3], drhodz [kg/m4], N2 [1/s2], dN2dz [1/m*s2]')


##########################################
# PARAMETERS FOR TS CALCULATION
##########################################

eta = 20 
N = 0.01 # Desired buoyancy frequency
H = 9000 # Input height in meters
f = 0.00012 # Coriolis parameter
LAMBDA = 10 # shear of the background (horizontally averaged) zonal velocity (same value as in fortran code)
###########################
Umax = 60 # Maximum velocity in m/s (same value as in fortran code)
###########################
print(f"Ts = {eta * N * H / f/Umax/3600/24} days") # TO CHANGE AT FORTRAN CODE
 