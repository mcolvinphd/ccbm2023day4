# Simple Soft Sphere Molecular Dynamics Program for Chem 160/260
##import mdfunctions here
##Insert line from lecture notes
from mdfunctions import *
import math
import random
import copy

# Function to produce a random unit vector
def random_unit():
    z=random.random()*2.0-1
    theta=random.random()*2.*math.pi
    x=math.sqrt(1-z*z)*math.cos(theta)
    y=math.sqrt(1-z*z)*math.sin(theta)
    return [x,y,z]

# Set simulation parameters
natoms=8**3
mass=39.948     #AMU (g/mole)
sigma=.341      #nm
epsilon=0.9661  #kj/mole equals nm^2 AMU/ps^2
T=240.           #K
tstep=0.001     #ps
tstep2=tstep*tstep #ps^2
nsteps=50      #Simulation steps
e_intvl=2      #Interval between calculating E

# Convert units and precalculate constants
rm=2**(1/6)*sigma
sep=1.*rm
sigma6=sigma**6
sigma12=sigma**12
k=0.00831446           #KJ/mol*K equals nm^2 AMU/ps^2*K
speed2=3.*k*T/mass  #nm^2/ps^2
speed=math.sqrt(speed2)
side=int(natoms**(1./3.)-.01)+1

# Create lists to hold values
pos=[]
pos_new=[]
vel=[]
vel_new=[]
force=[]
force_new=[]

# Initialize all of the simulation lists
iatom=0;
for x in range(side):
    for y in range(side):
        for z in range(side):
            pos.append([x*sep,y*sep,z*sep])
            pos_new.append([0.,0.,0.])
            rvect=random_unit()
            vel.append([rvect[0]*speed,rvect[1]*speed,rvect[2]*speed])
            vel_new.append([0.,0.,0.])
            force.append([0.,0.,0.])
            force_new.append([0.,0.,0.])
            iatom+=1
            if (iatom==natoms):
                break
        else:              # Some fancy python to break from all 3 loops
            continue       # when we've added all the atoms we need
        break
    else:
        continue
    break

# Run dynamics
for step in range(nsteps):
    #Update positions using velocity verlet integration and zero out force matrix
    for iatom in range(natoms):
        for dim in range(3):
            pos_new[iatom][dim]=pos[iatom][dim]+vel[iatom][dim]*tstep+force[iatom][dim]/mass*tstep2
            force_new[iatom][dim]=0.0

    #Calculate force on each atom
    for iatom in range(natoms-1):    
        for jatom in range(iatom+1,natoms):      #Loop over unique pairs of atoms
            r=dist(pos_new[iatom],pos_new[jatom])
            r7=r**(-7)
            r13=r**(-13)
            forcemag=24.*epsilon*(-2.*sigma12*r13+sigma6*r7)
            # Resolve force into vector components
            unit=vector(pos_new[iatom],pos_new[jatom])
            for idim in range(3):
                force_new[iatom][idim]-=unit[idim]*forcemag
                force_new[jatom][idim]+=unit[idim]*forcemag

    #Update velcities using velocity verlet integration
    for iatom in range(natoms):
        for dim in range(3):
            vel_new[iatom][dim]=vel[iatom][dim]+0.5/mass*(force[iatom][dim]+force_new[iatom][dim])*tstep

    #Calculate energy 
    if step%e_intvl==0:
        PE=0.
        KE=0.
        for iatom in range(natoms-1):    
            KE+=mass*sqmag(vel_new[iatom])/2.0
            for jatom in range(iatom+1,natoms):      #Loop over unique pairs of atoms
                r=dist(pos_new[iatom],pos_new[jatom])
                r6=r**(-6)
                r12=r**(-12)
                PE+=4.*epsilon*(sigma12*r12-sigma6*r6)
        KE+=mass*sqmag(vel_new[natoms-1])/2.0  # Add final contribution that's missed in loops
        print("     At step %5d PE is %7.2f  KE is %7.2f TE is %7.2f"%(step,PE,KE,PE+KE))

    #Now shift forward a timestep copying lists
    pos=copy.deepcopy(pos_new)
    vel=copy.deepcopy(vel_new)
    force=copy.deepcopy(force_new)

        
