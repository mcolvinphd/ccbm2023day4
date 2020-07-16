from math import sqrt
from drawtraj import drawtraj
def force(x,y,m,mstar):
    r2=x**2+y**2
    r32=r2*sqrt(r2)
    fx=-x*m*mstar/r32
    fy=-y*m*mstar/r32
    return fx,fy
def integrate(x,y,vx,vy,fx,fy,m,dt):
    ## Insert lines from lecture notes here

# Main part of the program
mstar=100
m=1
nsteps=100000
dt=0.01
r=50
x,y=0,r
vx,vy=1.2,0
trajx,trajy=[],[]
## Insert lines from lecture notes here
drawtraj(trajx,trajy,5*r)
