import plotly.graph_objects as go
from math import pi,sqrt
import numpy as np

# Defining the primitive cell and tight-binding hopping terms for the system (angstroms and eV)
C = 17.743/3
A = 2.830
TNN = 1
TNNN = 0.15
TZZ = 0.042

# Defining the crystallographic vectors for the unit cell
_AV = np.array([1,0,0])*A
_BV = np.array([np.cos(2/3*pi),np.cos(1/6*pi),0])*A
_CV = np.array([0,0,1])*C

def Pd_energy(kx,ky,kz):
    """ PdCoO2 energy dispersion
    
    Arguments:
        kx {float} -- wave vector
        ky {float} -- wave vector
        kz {float} -- wave vector
    
    Returns:
        float -- energy (eV)
    """
    k = np.moveaxis(np.array([kx,ky,kz]),0,-1)

    energy_nn = -2*TNN*(np.cos(np.dot(k,_AV))+np.cos(np.dot(k,_BV))+np.cos(np.dot(k,_AV+_BV)))-2*TZZ*np.cos(np.dot(k,_CV))
    energy_nnn = -2*TNNN*(np.cos(np.dot(k,_AV))**2+np.cos(np.dot(k,_BV))**2+np.cos(np.dot(k,_AV+_BV))**2)
    energy_t = energy_nn + energy_nnn
    return energy_t

#################################################
C = C*3
def hexagon_vertices(shift=0, scaling = 4*pi/(3*A)):
    # Reciprocal lattice vectors
    v1 = np.array([1,0,0])*scaling
    v2 = np.array([1/2,sqrt(3)/2,0])*scaling
    vertices_vec = np.array([v1,v2,v2-v1,-v1,-v2,v1-v2,v1])+shift
    return vertices_vec.T

line_settings = dict(
    showlegend=False,
    mode = 'lines',
    line_color = 'black',
    line_width = 4
    )

def hexagon_plot(shift=0,scaling=4*pi/(3*A)):
    X, Y, Z = hexagon_vertices(shift=shift,scaling=scaling)
    
    trace = go.Scatter3d(
        x=X, y=Y, z=Z,
        **line_settings
    )
    
    return trace

def hexagon_brillouin():
    vertices = hexagon_vertices(shift = np.array([0,0,-pi/C])).T
    data = []
    for i in vertices[:-1]:
        line = np.array([i,i+np.array([0,0,2*pi/C])])
        line = line.T
        
        trace = go.Scatter3d(
            x=line[0], y=line[1], z=line[2],
            **line_settings
        )
        
        data.append(trace)
    data.append(hexagon_plot(shift=np.array([0,0,pi/C])))
    data.append(hexagon_plot(shift=np.array([0,0,-pi/C])))
    return data

########################################

# Produces the lattice to be fed into plotly
def lattice_generator(a1,a2,N=6):
    grid = np.arange(-N//2,N//2,1)
    xGrid, yGrid = np.meshgrid(grid,grid)
    return np.transpose(np.reshape(np.kron(xGrid.flatten(),a1),(-1,2))+np.reshape(np.kron(yGrid.flatten(),a2),(-1,2)))

        