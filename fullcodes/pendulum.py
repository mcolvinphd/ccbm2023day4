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
while t<tmax:
    f=-k*x
    a=f/mass
    x+=v*dt
    v+=a*dt
    x_exact=cos(2*pi*t/period)
    rms_error+=(x-x_exact)**2
    steps+=1
    t+=dt
print("dt=%f RMS error=%f"%(dt, sqrt(rms_error/steps)))
