import numpy as n
from numpy import (sqrt, linspace, asarray,
                    empty, hstack , vstack,
                    pi, meshgrid)
from sys import argv
import matplotlib.pyplot as p
from mpl_toolkits import mplot3d

'''
Applies eccentricity to mesh
Requires the eccentricity be specified (as number, not percentage!)

ecc = (tmax - tmin) / (tmax + tmin)

Define: tmax = tg + D
        tmin = tg - D
Then:   D = ecc*tg
        d(rho) = D*( OD(z) - rho )/(OD(z) - ID)
        Xnew = X + d
Therefore, (x = +OD/2, y = 0) is the thinnest wall, and
(x = -OD/2, y=0) is the thickest wall.
'''

try:
    ecc = float(argv[1])
except IndexError:
    ecc = 0

# [0] rho, [1]Z, [2]Theta, [3] NodeNum
NC = n.load('./ConstructionFiles/nc_all.npy')
nc_cors = n.load('./ConstructionFiles/nc_cors.npy')
nc_med = n.load('./ConstructionFiles/nc_med.npy')
nc_fine = n.load('./ConstructionFiles/nc_fine.npy')

ID = n.min(nc_med[:,0])
if ID != n.min(nc_cors[:,0]):
    raise ValueError('Bad node coords.  Not all IDs the same.')
ODtop = n.max(nc_cors[:,0]) 
# Outstanding question:  Do we want the ecc shift to base on a the perfect
# thickness, or the imperfection-reduced thickness?
tg = n.max(nc_fine[:,0]) - ID 
Lg = .4/2
R = .125
coord_end_chamf = n.sqrt( R**2 - (ODtop-(ID+tg+R))**2 ) + Lg  # y-coord where chamfer reaches ODtop
D = ecc*tg

del nc_cors
del nc_fine
del nc_med

def CalcOD(Y):
    X = empty(Y.shape)
    rgn1 = (0<=Y) & (Y<=Lg)
    rgn2 = (Lg<Y) & (Y<=coord_end_chamf)
    rgn3 = (coord_end_chamf<Y)
    X[rgn1] = ID + tg
    X[rgn2] = -sqrt( R**2 - (Y[rgn2]-Lg)**2 ) + (ID+tg+R)
    X[rgn3] = ODtop
    return X

X = NC[:,0]*n.cos(NC[:,2])  
Y = NC[:,0]*n.sin(NC[:,2])

OD_all = CalcOD(NC[:,1])

d = D*(OD_all - NC[:,0])/(OD_all - ID)
X+=d

data = n.vstack((NC[:,-1],X,Y,NC[:,1])).T
#n.savetxt('nodes_for_abaqus.dat',X=data,fmt='%.0f, %.12f, %.12f, %.12f')
n.save('./ConstructionFiles/abaqus_nodes.npy',data)
