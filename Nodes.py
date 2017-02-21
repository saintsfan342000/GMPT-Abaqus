import numpy as n
from numpy import (sqrt, linspace, asarray,
                    empty, hstack , vstack,
                    pi, meshgrid)
from sys import argv
import matplotlib.pyplot as p
from mpl_toolkits import mplot3d

'''
Node coordinates for 3D TT model
Specify Lg, Ltop, ODtop, ID, tg, numleth either on command line or in script.
Works by laying out the nodes in "orthogonal cylindrical" coordinates, 
then at the end using x=rcos(q),y-rsin(q) to map to a cylinder in cartesian space

In addition to a unique node number, each node gets an i,j,k index based on its position
in r, q, z space (referred to herein as x,z,y, respectively).

x or th, [:,0] = thru thickness, or radial coordinate
y, [:,1] = Axial coord
z or q, [:,2] = Angle coordinate theta
'''
try:
    worthless, Lg, Ltop, ODtop, ID, tg, R, num_el_fine_r, dt = argv
    Lg = float(Lg)  # Half-length of gage section
    Ltop = float(Ltop)  # Length of thick section above radius/chamf
    ODtop = float(ODtop)/2    # Outer RADIUS of thick section
    ID = float(ID)/2  # Inner RADIUS
    tg = float(tg)  # mean thickness of test sxn.  Must mean = 0.5*(tmin + tmax)
    R = float(R)    # Radius of chamf
    num_el_fine_r = int(num_el_fine_r)  # Num elements thru the test section thickness
    dt = float(dt) # Specify dt/tg.  
except:
    Lg = 4 / 2
    Ltop = 0.5  # Length of thick section above radius/chamf
    ODtop = 1.9685/2    # Radius of thick section
    ID = 1.6285/2 # Inner Radius
    tg = .05
    R = .125    # Radius of chamf
    num_el_fine_r = 3 # Num elements thru the test section thicknes
    dt = 0

angle = 1*pi    # 2pi for full tube
coord_end_chamf = sqrt( R**2 - (ODtop-(ID+tg+R))**2 ) + Lg  # y-coord where chamfer reaches ODtop
ztop = coord_end_chamf + Ltop
start_ref1 = 2*Lg/3  # y-coord of last node in fine-mesh elements, where the refine starts
start_ref2 = coord_end_chamf  # y-coord of last node in med-mesh elements, where the refine starts

def CalcOD(Y):
    X = empty(Y.shape)
    rgn1 = (0<=Y) & (Y<=Lg)
    rgn2 = (Lg<Y) & (Y<=coord_end_chamf)
    rgn3 = (coord_end_chamf<Y)
    X[rgn1] = ID + tg
    X[rgn2] = -sqrt( R**2 - (Y[rgn2]-Lg)**2 ) + (ID+tg+R)
    X[rgn3] = ODtop
    return X

# Fine mesh definitions
elht_fine = tg / num_el_fine_r # Elemt height equal to element width
num_el_fine_z = int(start_ref1 // elht_fine + 1)   # Number of fine-mesh elements in z-dir
num_node_fine_z = num_el_fine_z + 1    # Num nodes in z-dir up to START of refine
num_node_fine_r = num_el_fine_r + 1
dq_fine = elht_fine/ID  # Estimate angular step of an element
num_el_fine_q = int(angle//dq_fine)
if num_el_fine_q%9 != 0: 
    num_el_fine_q-=(num_el_fine_q%9)
# Adjust num if angle != 2pi
if angle != 2*pi:
   num_el_fine_q+=1
dq_fine = angle/(num_el_fine_q-1)   # Corrected angular step
# medium mesh definitions
num_el_med_r = int(num_el_fine_r/3)
num_node_med_r = num_el_med_r + 1
elht_med = 3*elht_fine
start_med_z = start_ref1 + 2*elht_med
num_el_med_z = int((start_ref2 - start_med_z)//elht_med)
num_node_med_z = num_el_med_z + 1
# Coarse mesh definitions
num_el_cors_r = num_el_med_r # Only refining height wise and circumfwise, not thru-thickness
num_node_cors_r = num_el_cors_r + 1
elht_cors = 3*elht_med
start_cors_z = start_ref2 + 1*elht_cors   # Again, not a thickness
num_el_cors_z = int((ztop - start_cors_z) // elht_cors)   # Number of elements up height in coarse rgn
num_node_cors_z = num_el_cors_z + 1   # Num nodes in z-dir from start of coarse to top

# Layout the fine-mesh nodes
zspace = linspace(0, start_ref1, num_node_fine_z)   # Linspace of z-coordinates of fine nodes, up to start_ref1
zindspace = n.arange(len(zspace))   # Corresponding indices
if angle == 2*pi:
    qspace = n.linspace(0, angle-dq_fine, num_el_fine_q)  # Linspace of theta-coordinates
else:
    qspace = n.linspace(0, angle, num_el_fine_q)
qindspace = n.arange(len(qspace))
rindspace = n.arange(num_node_fine_r)
OD_fine = CalcOD(zspace)
nc_fine = empty((0,3))
ni_fine = empty((0,3),dtype=int)
for k,z in enumerate(zspace):
    # Node coords
    ro = ID+(.5*(dt*tg)*(1+n.cos(-pi*z/start_ref1))) # Here's the imperfection! 
    r = linspace(ro, OD_fine[k], num_node_fine_r)
    nc_fine = vstack(( nc_fine, vstack(map(n.ravel,n.meshgrid(r,z,qspace))).T ))
    # Node indices
    nod = asarray([zindspace[k]])
    ni_fine = vstack(( ni_fine, vstack(map(n.ravel,n.meshgrid(rindspace,k,qindspace))).T ))

nodenums = (n.arange(len(nc_fine))+1)[:,None]
nc_fine = hstack((nc_fine,nodenums))
ni_fine = hstack((ni_fine,nodenums))

# Now create the 3d-indexed array to store the node numbers
ni3d_fine = n.empty((len(rindspace),len(zindspace),len(qindspace)))
for v in ni_fine:
    ni3d_fine[v[0],v[1],v[2]] = v[3]

# Medium nodes    
zspace = linspace(start_med_z, start_ref2, num_node_med_z)
zindspace = n.arange(len(zspace))
qspace = qspace[::3] # Keep every third!
qindspace = n.arange(len(qspace))
rindspace = n.arange(num_node_med_r)
OD_med = CalcOD(zspace)
nc_med = empty((0,3))
ni_med = empty((0,3),dtype=int)
for k,z in enumerate(zspace):
    # Node coords
    r = linspace(ID, OD_med[k], num_node_med_r)
    nc_med = vstack(( nc_med, vstack(map(n.ravel,n.meshgrid(r,z,qspace))).T ))
    # Node indices
    nod = asarray([zindspace[k]])
    ni_med = vstack(( ni_med, vstack(map(n.ravel,n.meshgrid(rindspace,k,qindspace))).T ))

nodenums = (n.arange(len(nc_med)) + n.max(nodenums))[:,None]+1
nc_med = hstack((nc_med,nodenums))
ni_med = hstack((ni_med,nodenums))

ni3d_med = n.empty((len(rindspace),len(zindspace),len(qindspace)))
for v in ni_med:
    ni3d_med[v[0],v[1],v[2]] = v[3]

# coarse nodes    
zspace = linspace(start_cors_z, ztop, num_node_cors_z)   # Linspace of z-coordinates of cors nodes, up to start_ref1
zindspace = n.arange(len(zspace))
qspace = qspace[::3] # Keep every third!
qindspace = n.arange(len(qspace))
rindspace = n.arange(num_node_cors_r)
OD_cors = CalcOD(zspace)
nc_cors = empty((0,3))
ni_cors = empty((0,3),dtype=int)
for k,z in enumerate(zspace):
    # Node coords
    r = linspace(ID, OD_cors[k], num_node_cors_r)
    nc_cors = vstack(( nc_cors, vstack(map(n.ravel,n.meshgrid(r,z,qspace))).T ))
    # Node indices
    nod = asarray([zindspace[k]])
    ni_cors = vstack(( ni_cors, vstack(map(n.ravel,n.meshgrid(rindspace,k,qindspace))).T ))

nodenums = (n.arange(len(nc_cors)) + n.max(nodenums))[:,None]+1
nc_cors = hstack((nc_cors,nodenums))
ni_cors = hstack((ni_cors,nodenums))        

ni3d_cors = n.empty((len(rindspace),len(zindspace),len(qindspace)))
for v in ni_cors:
    ni3d_cors[v[0],v[1],v[2]] = v[3]
    
# First refine (start_ref1) circumferential (q) nodes
z = start_ref1 + elht_med/2
# locs are those fine nodes that are at the top z-coord of the fine section
locs = nc_fine[:,1] == n.max(nc_fine[:,1])
r = n.sort( n.unique( nc_fine[locs,0] )) 
q = n.sort( n.unique( nc_fine[locs,2] ))
nc_ref1_q = vstack(map(n.ravel,n.meshgrid(r,z,q))).T
# Indices
z = n.array([0])
r = n.arange(len(r))
q = n.arange(len(q))
ni_ref1_q = vstack(map(n.ravel,n.meshgrid(r,z,q))).T
# Now take out node with q-index = 0,3,6,9,...
nc_ref1_q = nc_ref1_q[ ni_ref1_q[:,2]%3 != 0, :]
ni_ref1_q = ni_ref1_q[ ni_ref1_q[:,2]%3 != 0, :]
# nodenums
nodenums = (n.arange(len(nc_ref1_q)) + n.max(nodenums))[:,None]+1
nc_ref1_q = hstack((nc_ref1_q,nodenums))
ni_ref1_q = hstack((ni_ref1_q,nodenums))        
#3d index array
ni3d_ref1_q = n.empty((len(r),len(z),len(q)))
for v in ni_ref1_q:
    ni3d_ref1_q[v[0],v[1],v[2]] = v[3]

# First refine (start_ref1)  mid-level nodes (btwn the circumf. and the thickness refine)
# thickness (r) like fine, circumf. (q) like med
z = start_ref1 + elht_med
r = n.sort( n.unique( nc_fine[locs,0] ))  # locs defined above to account for imperfection
q = n.sort( n.unique( nc_med[:,2] ))
nc_ref1_mid = vstack(map(n.ravel,n.meshgrid(r,z,q))).T
# Indices
z = n.array([0])
r = n.arange(len(r))
q = n.arange(len(q))
ni_ref1_mid = vstack(map(n.ravel,n.meshgrid(r,z,q))).T
#nodenums
nodenums = (n.arange(len(nc_ref1_mid)) + n.max(nodenums))[:,None]+1
nc_ref1_mid = hstack((nc_ref1_mid,nodenums))
ni_ref1_mid = hstack((ni_ref1_mid,nodenums))
#3d index array
ni3d_ref1_mid = n.empty((len(r),len(z),len(q)))
for v in ni_ref1_mid:
    ni3d_ref1_mid[v[0],v[1],v[2]] = v[3]

# First refine (start_ref1) thickness (th, r-dir'n) nodes
z = start_ref1 + 3*elht_med/2
r = n.sort( n.unique( nc_fine[locs,0] )) # locs defined above to account for imperfection
q = n.sort( n.unique( nc_med[ni_med[:,1]==0,2] ))
nc_ref1_r = vstack(map(n.ravel,n.meshgrid(r,z,q))).T
# Indices
z = n.array([0])
r = n.arange(len(r))
q = n.arange(len(q))
ni_ref1_r = vstack(map(n.ravel,n.meshgrid(r,z,q))).T
# Now take out node with r-index = 0,3,6,9,...
nc_ref1_r = nc_ref1_r[ ni_ref1_r[:,0]%3 != 0, :]
ni_ref1_r = ni_ref1_r[ ni_ref1_r[:,0]%3 != 0, :]
#nodenums
nodenums = (n.arange(len(nc_ref1_r)) + n.max(nodenums))[:,None]+1
nc_ref1_r = hstack((nc_ref1_r,nodenums))
ni_ref1_r = hstack((ni_ref1_r,nodenums))     
#3d index array
ni3d_ref1_r = n.empty((len(r),len(z),len(q)))
for v in ni_ref1_r:
    ni3d_ref1_r[v[0],v[1],v[2]] = v[3]


# Second refine (start_ref2) circumferential (z) nodes
z = start_ref2 + elht_cors/2
# med crosses the radius from test sxn to ends, so unique won't work unless we
# also specify that we want the nodes at the top of the med section
r = n.sort( n.unique( nc_med[nc_med[:,1]==start_ref2,0] ))
q = n.sort( n.unique( nc_med[nc_med[:,1]==start_ref2,2] ))
nc_ref2_q = vstack(map(n.ravel,n.meshgrid(r,z,q))).T
# Indices
z = n.array([0])
r = n.arange(len(r))
q = n.arange(len(q))
ni_ref2_q = vstack(map(n.ravel,n.meshgrid(r,z,q))).T
# Now take out node with q-index = 0,3,6,9,...
nc_ref2_q = nc_ref2_q[ ni_ref2_q[:,2]%3 != 0, :]
ni_ref2_q = ni_ref2_q[ ni_ref2_q[:,2]%3 != 0, :]
#nodenums
nodenums = (n.arange(len(nc_ref2_q)) + n.max(nodenums))[:,None]+1
nc_ref2_q = hstack((nc_ref2_q,nodenums))
ni_ref2_q = hstack((ni_ref2_q,nodenums))   
#3d index array
ni3d_ref2_q = n.empty((len(r),len(z),len(q)))
for v in ni_ref2_q:
    ni3d_ref2_q[v[0],v[1],v[2]] = v[3]

# All node coords
NC = vstack((nc_fine,nc_med,nc_cors,nc_ref1_q,nc_ref1_mid,nc_ref1_r,nc_ref2_q))
NI = vstack((ni_fine,ni_med,ni_cors,ni_ref1_q,ni_ref1_mid,ni_ref1_r,ni_ref2_q))

'''
for k in 1:
    ## View all nodes at q = 0, looking down the circumferential dir
    %matplotlib
    p.figure()
    rgn = NI[:,2]==0
    p.plot(NC[rgn,0],NC[rgn,1],'.')
    for k,i in enumerate(NC[rgn,0]):
        p.text(i,NC[rgn][k,1],'({:.0f},{:.0f})'.format(NI[rgn][k,0],NI[rgn][k,1]),ha='center',va='center',alpha=0.0)

    # Plot all thru-thickness nodes but just a small, fast-to-plot rgn of the other two directions
    # and write the node's indices 
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)
    rgn = (ni_fine[:,1] <= 1)&(ni_fine[:,2]<=3)
    ax.plot(nc_fine[rgn,0],nc_fine[rgn,2],nc_fine[rgn,1],'b.',alpha=0.6)
    C = ni_fine[rgn].copy()
    for k,i in enumerate(nc_fine[rgn]):
        ax.text(i[0],i[2],i[1],'{},{},{}'.format(C[k,0],C[k,2],C[k,1]))
        
    # Demonstrate ref1
    p.close('all')
    %matplotlib
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)
    rgn = ((ni_fine[:,1] == n.max(ni_fine[:,1])) &
            (ni_fine[:,2] <= 9) &
            (ni_fine[:,0] <= 3))
    qmax = n.max(nc_fine[rgn,2])
    rmax = n.max(nc_fine[rgn,0])
    ax.plot(nc_fine[rgn,0],nc_fine[rgn,2],nc_fine[rgn,1],'b.',alpha=0.5)
    for k,i in enumerate(nc_fine[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(ni_fine[rgn][k,0],
                                            ni_fine[rgn][k,1],
                                            ni_fine[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='b',fontsize=10)
    rgn = ((nc_ref1_q[:,2] <= qmax)&
            (nc_ref1_q[:,0] <= rmax))
    ax.plot(nc_ref1_q[rgn,0],nc_ref1_q[rgn,2],nc_ref1_q[rgn,1],'r.',alpha=0.5)
    for k,i in enumerate(nc_ref1_q[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(ni_ref1_q[rgn][k,0],
                                            ni_ref1_q[rgn][k,1],
                                            ni_ref1_q[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='r',fontsize=10)
    rgn = ((nc_ref1_mid[:,2] <= qmax)&
            (nc_ref1_mid[:,0] <= rmax))
    ax.plot(nc_ref1_mid[rgn,0],nc_ref1_mid[rgn,2],nc_ref1_mid[rgn,1],'g.',alpha=0.5)
    for k,i in enumerate(nc_ref1_mid[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(ni_ref1_mid[rgn][k,0],
                                            ni_ref1_mid[rgn][k,1],
                                            ni_ref1_mid[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='g',fontsize=10)
    rgn = ((nc_ref1_r[:,2] <= qmax)&
            (nc_ref1_r[:,0] <= rmax))
    ax.plot(nc_ref1_r[rgn,0],nc_ref1_r[rgn,2],nc_ref1_r[rgn,1],'k.',alpha=0.5)
    for k,i in enumerate(nc_ref1_r[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(ni_ref1_r[rgn][k,0],
                                            ni_ref1_r[rgn][k,1],
                                            ni_ref1_r[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='k',fontsize=10)
    rgn = ((nc_med[:,2] <= qmax)&
            (nc_med[:,0] <= rmax)&
            (ni_med[:,1] == n.min(ni_med[:,1])))
    ax.plot(nc_med[rgn,0],nc_med[rgn,2],nc_med[rgn,1],'.',alpha=1,color='orange')
    for k,i in enumerate(nc_med[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(ni_med[rgn][k,0],
                                            ni_med[rgn][k,1],
                                            ni_med[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=1,color='orange',fontsize=10)
    p.tight_layout()        

    # Demonstrate ref2
    p.close('all')
    %matplotlib
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Y',labelpad=30)
    rgn = ((ni_med[:,1] == n.max(ni_med[:,1])) &
            (ni_med[:,2] <= 9) &
            (ni_med[:,0] <= 3))
    qmax = n.max(nc_med[rgn,2])
    rmax = n.max(nc_med[rgn,0])
    ax.plot(nc_med[rgn,0],nc_med[rgn,2],nc_med[rgn,1],'b.',alpha=0.5)
    for k,i in enumerate(nc_med[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(ni_med[rgn][k,0],
                                            ni_med[rgn][k,1],
                                            ni_med[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='b',fontsize=10)
    rgn = ((nc_ref2_q[:,2] <= qmax)&
            (nc_ref2_q[:,0] <= rmax))
    ax.plot(nc_ref2_q[rgn,0],nc_ref2_q[rgn,2],nc_ref2_q[rgn,1],'r.',alpha=0.5)
    for k,i in enumerate(nc_ref2_q[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(ni_ref2_q[rgn][k,0],
                                            ni_ref2_q[rgn][k,1],
                                            ni_ref2_q[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=0.5,color='r',fontsize=10)
    rgn = ((nc_cors[:,2] <= qmax)&
            (nc_cors[:,0] <= rmax)&
            (ni_cors[:,1] == n.min(ni_cors[:,1])))
    ax.plot(nc_cors[rgn,0],nc_cors[rgn,2],nc_cors[rgn,1],'.',alpha=1,color='g')
    for k,i in enumerate(nc_cors[rgn]):
        text = '{:.0f},{:.0f},{:.0f}'.format(ni_cors[rgn][k,0],
                                            ni_cors[rgn][k,1],
                                            ni_cors[rgn][k,2])
        ax.text(i[0],i[2],i[1],text,alpha=1,color='g',fontsize=10)        
    p.tight_layout()        
        

    # Demonstrate node index order
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)
    rgn = (ni_fine[:,1] <= 1)&(ni_fine[:,2]<=3)
    ax.plot(nc_fine[rgn,0],nc_fine[rgn,1],nc_fine[rgn,2],'b.',alpha=0.6)
    C = ni_fine[rgn].copy()
    for k,i in enumerate(nc_fine[rgn]):
        ax.text(i[0],i[1],i[2],'{},{},{}'.format(C[k,0],C[k,1],C[k,2]))
    p.tight_layout()
   
for k in [1]:
    # Plot from ID looking out (x == 0)
    rgn = NI[:,0] == 0
    p.plot(NC[rgn,2],NC[rgn,1],'.')
    rgn = (nc_ref1_q[:,0] == ID)
    p.plot(nc_ref1_q[rgn][:,2],nc_ref1_q[rgn][:,1],'r.')
    rgn = (nc_ref1_mid[:,0] == ID)
    p.plot(nc_ref1_mid[rgn][:,2],nc_ref1_mid[rgn][:,1],'g.')
    rgn = (ni_ref1_r[:,0] == 1)
    p.plot(nc_ref1_r[rgn][:,2],nc_ref1_r[rgn][:,1],'yo')
    rgn = (nc_ref2_q[:,0] == ID)
    p.plot(nc_ref2_q[rgn][:,2],nc_ref2_q[rgn][:,1],'r.')
    
    # Plot looking in hoop direction
    p.figure()
    rgn = NI[:,2] == 0
    p.plot(NC[rgn,0],NC[rgn,1],'.')
    rgn = (ni_ref1_q[:,2] == 1)
    p.plot(nc_ref1_q[rgn][:,0],nc_ref1_q[rgn][:,1],'r.')
    rgn = (ni_ref1_mid[:,2] == 0)
    p.plot(nc_ref1_mid[rgn][:,0],nc_ref1_mid[rgn][:,1],'g.')
    rgn = (ni_ref1_r[:,2] == 0)
    p.plot(nc_ref1_r[rgn][:,0],nc_ref1_r[rgn][:,1],'yo')
    rgn = (ni_ref2_q[:,2] == 1)
    p.plot(nc_ref2_q[rgn][:,0],nc_ref2_q[rgn][:,1],'r.')

    # 3D Polar
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    rgn = NC[:,2] <= pi/2
    c = NC[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'.',alpha=0.15)
    rgn = nc_ref1_q[:,2] <= pi/2
    c = nc_ref1_q[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'r.',alpha=0.5)
    rgn = nc_ref1_mid[:,2] <= pi/2
    c = nc_ref1_mid[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'g.',alpha=0.5)
    rgn = nc_ref1_r[:,2] <= pi/2
    c = nc_ref1_r[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'k.',alpha=0.5)
    rgn = nc_ref2_q[:,2] <= pi/2
    c = nc_ref2_q[rgn]
    x, y, z = c[:,0]*n.cos(c[:,2]), c[:,0]*n.sin(c[:,2]), c[:,1]
    ax.plot(x,y,z,'r.',alpha=0.5)

   
    NC_polar_X = NC[:,0]*n.cos(NC[:,2])
    NC_polar_Y = NC[:,0]*n.sin(NC[:,2])
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.plot(NC_polar_X,NC_polar_Y,NC[:,1],'.',alpha=0.3)
    ax.set_xlabel('x')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    p.show()

    p.figure()
    p.plot(NC_polar_X,NC_polar_Y,'.',alpha=0.3)
    p.axis('equal')
'''

nodelists = ['nc_cors', 'nc_fine', 'nc_med',
            'nc_ref1_mid', 'nc_ref1_r', 'nc_ref1_q',
            'nc_ref2_q', 'ni_cors', 'ni_fine',
            'ni_med', 'ni_ref1_mid', 'ni_ref1_r',
            'ni_ref1_q', 'ni_ref2_q', 'NC','NI',
            'ni3d_fine', 'ni3d_med', 'ni3d_cors', 'ni3d_ref1_mid', 
            'ni3d_ref1_r', 'ni3d_ref1_q', 'ni3d_ref2_q']
 
 
for k, F in enumerate(nodelists):
    if F.rfind('nc') == 0:
        n.save('./ConstructionFiles/{}'.format(F), eval(F))
    elif F.rfind('ni') == 0:
        n.save('./ConstructionFiles/{}'.format(F), eval(F))
    elif F.rfind('NC') == 0:
        #n.savetxt('./ConstructionFiles/{}.dat'.format('node_coords_all'),X=eval(F),fmt = '%.12f,%.12f,%.12f,%.0f')
        n.save('./ConstructionFiles/{}'.format('nc_all'),eval(F)) 
    elif F.rfind('NI') == 0:
        #n.savetxt('./ConstructionFiles/{}.dat'.format('node_indices_all'),X=eval(F),fmt = '%.0f,%.0f,%.0f,%.0f')
        n.save('./ConstructionFiles/{}'.format('ni_all'),eval(F).astype(int)) 
    elif F.rfind('ni3d_') == 0:
        n.save('./ConstructionFiles/{}'.format(F),eval(F).astype(int)) 
    elif F.rfind('abaqus') == 5:
        '''
        # node list for abaqus input file
        X = NC[:,0]*n.cos(NC[:,2])  
        Y = NC[:,0]*n.sin(NC[:,2])
        data = n.vstack((NC[:,-1],X,Y,NC[:,1])).T
        n.savetxt('nodes.dat',X=data,fmt='%.0f, %.12f, %.12f, %.12f')        
        '''
        pass
    else:
        raise('Something bad in nodelist')
