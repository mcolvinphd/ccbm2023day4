from math import pi, cos, sqrt
dt=0.2
mass=10
period=10
omega=2*pi/period
k=mass*omega**2
tmax=period
x=+1
v=0
t=0
rms_error=0
steps=0
## Add lines here
print("dt=%f RMS error=%f"%(dt, sqrt(rms_error/steps)))

