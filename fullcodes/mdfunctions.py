import math
# Function to produce distance between two list points
def dist(a,b):
    return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)

# Function to produce square magnitude of vector
def sqmag(a):
    return (a[0]**2+a[1]**2+a[2]**2)

# Function to return unit vector pointing between two points
def vector(a,b):
    r=dist(a,b)
    x=(a[0]-b[0])/r
    y=(a[1]-b[1])/r
    z=(a[2]-b[2])/r
    return [x,y,z]
