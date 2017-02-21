import numpy as n
from numpy import (sqrt, linspace, asarray,
                    empty, hstack , vstack, pi)
from sys import argv

'''
r, [:,0] = thru thickness/radial
z, [:,1] = Axial coord
q, [:,2] = Angle coordinate theta (radians)
'''

ni3d_fine = n.load('./ConstructionFiles/ni3d_fine.npy')
ni3d_med = n.load('./ConstructionFiles/ni3d_med.npy')
ni3d_cors = n.load('./ConstructionFiles/ni3d_cors.npy')
ni3d_ref1_mid = n.load('./ConstructionFiles/ni3d_ref1_mid.npy')
ni3d_ref1_r = n.load('./ConstructionFiles/ni3d_ref1_r.npy')
ni3d_ref1_q = n.load('./ConstructionFiles/ni3d_ref1_q.npy')
ni3d_ref2_q = n.load('./ConstructionFiles/ni3d_ref2_q.npy')
nc_cors = n.load('./ConstructionFiles/nc_cors.npy')

# Check if this is full ring (2pi) or half model,
# as this affects connectivity
dq = n.diff(nc_cors[:,2]).max()
if n.isclose(nc_cors[:,2].max(),(2*pi-dq),rtol=.001):
    fullring = True
else:
    fullring = False
del nc_cors

# Mesh definitions
num_el_fine_r = ni3d_fine.shape[0] - 1
num_el_fine_z = ni3d_fine.shape[1] - 1
num_el_fine_q = ni3d_fine.shape[2]
if not fullring:  num_el_fine_q-=1
num_el_fine_tot = num_el_fine_q*num_el_fine_z*num_el_fine_r
num_el_med_r = ni3d_med.shape[0] - 1
num_el_med_z = ni3d_med.shape[1] - 1
num_el_med_q = ni3d_med.shape[2]
if not fullring:  num_el_med_q-=1
num_el_med_tot = num_el_med_q*num_el_med_z*num_el_med_r
num_el_cors_r = ni3d_cors.shape[0] - 1
num_el_cors_z = ni3d_cors.shape[1] - 1
num_el_cors_q = ni3d_cors.shape[2]
if not fullring:  num_el_cors_q-=1
num_el_cors_tot = num_el_cors_q*num_el_cors_z*num_el_cors_r

########################################################    
################## Connnect on Fine ##################
########################################################
# Notice, the inner-most loop is "for i in num_el_fine_r",
# meaning I will connect thru-thickness first, then in q/circ, then in z/axial
elcon_fine = n.empty((num_el_fine_tot,8))
row = 0
for j in range(num_el_fine_z):
    for k in range(num_el_fine_q):
        for i in range(num_el_fine_r):
            index = n.array([[i,j+1,k], 
                        [i+1,j+1,k],
                        [i+1,j+1,k+1],
                        [i,j+1,k+1],
                        [i,j,k],
                        [i+1,j,k],
                        [i+1,j,k+1],
                        [i,j,k+1]]).astype(int)
            # Connect last q-nodes back to q=0
            if (k == num_el_fine_q -1) and fullring:
                index[index[:,2] == k+1,2] = 0
            elcon_temp = [0,0,0,0,0,0,0,0]
            for u,v in enumerate(index):
                elcon_temp[u] = ni3d_fine[v[0],v[1],v[2]]
            elcon_fine[row] = elcon_temp
            row+=1
elnums = (n.arange(len(elcon_fine))+1)[:,None]
elcon_fine = hstack((elcon_fine,elnums)) 
print('Connected on Fine') 
########################################################    
################## Connnect on Med ##################
########################################################
elcon_med = n.empty((num_el_med_tot,8))
row = 0
for j in range(num_el_med_z):
    for k in range(num_el_med_q):
        for i in range(num_el_med_r):
            index = n.array([[i,j+1,k], 
                        [i+1,j+1,k],
                        [i+1,j+1,k+1],
                        [i,j+1,k+1],
                        [i,j,k],
                        [i+1,j,k],
                        [i+1,j,k+1],
                        [i,j,k+1]]).astype(int)
            elcon_temp = [0,0,0,0,0,0,0,0]
            if (k == num_el_med_q -1) and fullring:
                index[index[:,2] == k+1,2] = 0            
            for u,v in enumerate(index):
                elcon_temp[u] = ni3d_med[v[0],v[1],v[2]]
            elcon_med[row] = elcon_temp
            row+=1
elnums = (n.arange(len(elcon_med)) + n.max(elnums))[:,None]+1
elcon_med = hstack((elcon_med,elnums)) 
print('Connected on Med')
########################################################    
################## Connnect on Coarse ##################
########################################################
elcon_cors = n.empty((num_el_cors_tot,8))
row = 0
for j in (n.arange(num_el_cors_z)):
    for k in range(num_el_cors_q):
        for i in range(num_el_cors_r):
            index = n.array([[i,j+1,k], 
                        [i+1,j+1,k],
                        [i+1,j+1,k+1],
                        [i,j+1,k+1],
                        [i,j,k],
                        [i+1,j,k],
                        [i+1,j,k+1],
                        [i,j,k+1]]).astype(int)
            if (k == num_el_cors_q -1) and fullring:
                index[index[:,2] == k+1,2] = 0                        
            elcon_temp = [0,0,0,0,0,0,0,0]
            for u,v in enumerate(index):
                elcon_temp[u] = ni3d_cors[v[0],v[1],v[2]]
            elcon_cors[row] = elcon_temp
            row+=1
elnums = (n.arange(len(elcon_cors)) + n.max(elnums))[:,None]+1
elcon_cors = hstack((elcon_cors,elnums)) 
print('Connected on Coarse')
########################################################    
########## Connnect Ref1.  Fine to ref1_q ##############
########## and ref1_q to ref1_mid ######################
########################################################
# I have for i in... inside the if k%3 check so that I connect in
# complete thru-thickness stacks of each of the four similar-shape
# elements, then proceed circumferentially
j_fine = ni3d_fine.shape[1] - 1
j_mid = 0
j_ref = 0
elcon_ref1_q = n.empty((0,8))   # Too difficult to figure out the shape I need before hand 
for k in range(num_el_fine_q):
    if k%3 == 0:
        for i in range(num_el_fine_r): 
            nodes = [ni3d_ref1_mid[i,j_mid,k//3],
                     ni3d_ref1_mid[i+1,j_mid,k//3],
                     ni3d_ref1_q[i+1,j_ref,k+1],
                     ni3d_ref1_q[i,j_ref,k+1],
                     ni3d_fine[i,j_fine,k],
                     ni3d_fine[i+1,j_fine,k],
                     ni3d_fine[i+1,j_fine,k+1],
                     ni3d_fine[i,j_fine,k+1]
                    ]
            elcon_ref1_q = vstack((elcon_ref1_q,nodes))
    elif k%3 == 1:
        for i in range(num_el_fine_r):
            nodes = [ni3d_ref1_q[i,j_ref,k],
                     ni3d_ref1_q[i+1,j_ref,k],
                     ni3d_ref1_q[i+1,j_ref,k+1],
                     ni3d_ref1_q[i,j_ref,k+1],
                     ni3d_fine[i,j_fine,k],
                     ni3d_fine[i+1,j_fine,k],
                     ni3d_fine[i+1,j_fine,k+1],
                     ni3d_fine[i,j_fine,k+1]
                    ]
            elcon_ref1_q = vstack((elcon_ref1_q,nodes))
    elif k%3 == 2:
        for i in range(num_el_fine_r):
            # Exception for the last circumferential set to connect back to q=0
            if (k == num_el_fine_q - 1) and fullring:
                nodes = [ni3d_ref1_q[i,j_ref,k],
                         ni3d_ref1_q[i+1,j_ref,k],
                         ni3d_ref1_mid[i+1,j_mid,0],
                         ni3d_ref1_mid[i,j_mid,0],
                         ni3d_fine[i,j_fine,k],
                         ni3d_fine[i+1,j_fine,k],
                         ni3d_fine[i+1,j_fine,0],
                         ni3d_fine[i,j_fine,0]
                        ]
        
            else:   
                nodes = [ni3d_ref1_q[i,j_ref,k],
                         ni3d_ref1_q[i+1,j_ref,k],
                         ni3d_ref1_mid[i+1,j_mid,(k+1)//3],
                         ni3d_ref1_mid[i,j_mid,(k+1)//3],
                         ni3d_fine[i,j_fine,k],
                         ni3d_fine[i+1,j_fine,k],
                         ni3d_fine[i+1,j_fine,k+1],
                         ni3d_fine[i,j_fine,k+1]
                        ]
            elcon_ref1_q = vstack((elcon_ref1_q,nodes))
    # In a separate loop, connect ref1_q up to ref1_mid
    if k%3 == 0:
        for i in range(num_el_fine_r):            
            # Exception for the last the final z elements which need to 
            # connect back to the 0th znodes
            if (k == num_el_fine_q - 3) and fullring:
                nodes = [ni3d_ref1_mid[i,j_mid,k//3],
                         ni3d_ref1_mid[i+1,j_mid,k//3],
                         ni3d_ref1_mid[i+1,j_mid,0],
                         ni3d_ref1_mid[i,j_mid,0],
                         ni3d_ref1_q[i,j_ref,k+1],
                         ni3d_ref1_q[i+1,j_ref,k+1],
                         ni3d_ref1_q[i+1,j_ref,k+2],
                         ni3d_ref1_q[i,j_ref,k+2]
                        ]
                elcon_ref1_q = vstack((elcon_ref1_q,nodes))
            else:
                nodes = [ni3d_ref1_mid[i,j_mid,k//3],
                         ni3d_ref1_mid[i+1,j_mid,k//3],
                         ni3d_ref1_mid[i+1,j_mid,(k+3)//3],
                         ni3d_ref1_mid[i,j_mid,(k+3)//3],
                         ni3d_ref1_q[i,j_ref,k+1],
                         ni3d_ref1_q[i+1,j_ref,k+1],
                         ni3d_ref1_q[i+1,j_ref,k+2],
                         ni3d_ref1_q[i,j_ref,k+2]
                        ]
                elcon_ref1_q = vstack((elcon_ref1_q,nodes))    

elnums = (n.arange(len(elcon_ref1_q)) + n.max(elnums))[:,None]+1
elcon_ref1_q = hstack((elcon_ref1_q,elnums)) 

########################################################    
########## Connnect Ref1_mid to ref1_r ################
########## and ref1_r to med ######################
########################################################    
# Again, connect thru-thickness then proceed circumferentially
# Four loops for the four differently-shaped elements here
# 
j_med = 0
j_mid = 0
j_ref = 0
elcon_ref1_r = n.empty((0,8))
for k in range(num_el_med_q):
    for i in range(num_el_fine_r):
        # Exception for the last the final z elements which need to 
        # connect back to the 0th znodes
        if (k == num_el_med_q - 1) and fullring:
            k = -1 # We can do this since all elements have same z boundaries
        if i%3 == 0:
            nodes = [ni3d_med[i//3,j_med,k],
                     ni3d_ref1_r[i+1,j_ref,k],
                     ni3d_ref1_r[i+1,j_ref,k+1],
                     ni3d_med[i//3,j_med,k+1],
                     ni3d_ref1_mid[i,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k+1],
                     ni3d_ref1_mid[i,j_mid,k+1]
                    ]
            elcon_ref1_r = vstack((elcon_ref1_r,nodes))
        elif i%3 == 1:
            nodes = [ni3d_ref1_r[i,j_ref,k],
                     ni3d_ref1_r[i+1,j_ref,k],
                     ni3d_ref1_r[i+1,j_ref,k+1],
                     ni3d_ref1_r[i,j_ref,k+1],
                     ni3d_ref1_mid[i,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k+1],
                     ni3d_ref1_mid[i,j_mid,k+1]
                    ]
            elcon_ref1_r = vstack((elcon_ref1_r,nodes))
        elif i%3 == 2:
            nodes = [ni3d_ref1_r[i,j_ref,k],
                     ni3d_med[(i+1)//3,j_med,k],
                     ni3d_med[(i+1)//3,j_med,k+1],
                     ni3d_ref1_r[i,j_ref,k+1],
                     ni3d_ref1_mid[i,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k],
                     ni3d_ref1_mid[i+1,j_mid,k+1],
                     ni3d_ref1_mid[i,j_mid,k+1]
                    ]
            elcon_ref1_r = vstack((elcon_ref1_r,nodes))                
        if i%3 == 0:
            nodes = [ni3d_med[i//3,j_med,k],
                     ni3d_med[(i+3)//3,j_med,k],
                     ni3d_med[(i+3)//3,j_med,k+1],
                     ni3d_med[i//3,j_med,k+1],
                     ni3d_ref1_r[i+1,j_ref,k],
                     ni3d_ref1_r[i+2,j_ref,k],
                     ni3d_ref1_r[i+2,j_ref,k+1],
                     ni3d_ref1_r[i+1,j_ref,k+1]
                    ]
            elcon_ref1_r = vstack((elcon_ref1_r,nodes))
            
elnums = (n.arange(len(elcon_ref1_r)) + n.max(elnums))[:,None]+1
elcon_ref1_r = hstack((elcon_ref1_r,elnums))  
print('Connected on Ref1')
########################################################    
########## Connect med to ref2, ref2 to cors
########################################################
j_med = ni3d_med.shape[1] - 1
j_cors = 0
j_ref = 0
elcon_ref2 = n.empty((0,8))
for k in range(num_el_med_q):
    if k%3 == 0:
        for i in range(num_el_med_r): 
            nodes = [ni3d_cors[i,j_cors,k//3],
                     ni3d_cors[i+1,j_cors,k//3],
                     ni3d_ref2_q[i+1,j_ref,k+1],
                     ni3d_ref2_q[i,j_ref,k+1],
                     ni3d_med[i,j_med,k],
                     ni3d_med[i+1,j_med,k],
                     ni3d_med[i+1,j_med,k+1],
                     ni3d_med[i,j_med,k+1]
                    ]
            elcon_ref2 = vstack((elcon_ref2,nodes))
    elif k%3 == 1:
        for i in range(num_el_med_r):
            nodes = [ni3d_ref2_q[i,j_ref,k],
                     ni3d_ref2_q[i+1,j_ref,k],
                     ni3d_ref2_q[i+1,j_ref,k+1],
                     ni3d_ref2_q[i,j_ref,k+1],
                     ni3d_med[i,j_med,k],
                     ni3d_med[i+1,j_med,k],
                     ni3d_med[i+1,j_med,k+1],
                     ni3d_med[i,j_med,k+1]
                    ]
            elcon_ref2 = vstack((elcon_ref2,nodes))
    elif k%3 == 2:
        for i in range(num_el_med_r):
            # Exception for the last circumferential set
            if (k == num_el_med_q - 1) and fullring:
                nodes = [ni3d_ref2_q[i,j_ref,k],
                         ni3d_ref2_q[i+1,j_ref,k],
                         ni3d_cors[i+1,j_cors,0],
                         ni3d_cors[i,j_cors,0],
                         ni3d_med[i,j_med,k],
                         ni3d_med[i+1,j_med,k],
                         ni3d_med[i+1,j_med,0],
                         ni3d_med[i,j_med,0]
                        ]
        
            else:   
                nodes = [ni3d_ref2_q[i,j_ref,k],
                         ni3d_ref2_q[i+1,j_ref,k],
                         ni3d_cors[i+1,j_cors,(k+1)//3],
                         ni3d_cors[i,j_cors,(k+1)//3],
                         ni3d_med[i,j_med,k],
                         ni3d_med[i+1,j_med,k],
                         ni3d_med[i+1,j_med,k+1],
                         ni3d_med[i,j_med,k+1]
                        ]
            elcon_ref2 = vstack((elcon_ref2,nodes))
    # In a separate loop, connect ref2_q up to ref2_mid
    if k%3 == 0:
        for i in range(num_el_med_r):            
            if (k == num_el_med_q - 3) and fullring:
                nodes = [ni3d_cors[i,j_cors,k//3],
                         ni3d_cors[i+1,j_cors,k//3],
                         ni3d_cors[i+1,j_cors,0],
                         ni3d_cors[i,j_cors,0],
                         ni3d_ref2_q[i,j_ref,k+1],
                         ni3d_ref2_q[i+1,j_ref,k+1],
                         ni3d_ref2_q[i+1,j_ref,k+2],
                         ni3d_ref2_q[i,j_ref,k+2]
                        ]
                elcon_ref2 = vstack((elcon_ref2,nodes))
            else:
                nodes = [ni3d_cors[i,j_cors,k//3],
                         ni3d_cors[i+1,j_cors,k//3],
                         ni3d_cors[i+1,j_cors,(k+3)//3],
                         ni3d_cors[i,j_cors,(k+3)//3],
                         ni3d_ref2_q[i,j_ref,k+1],
                         ni3d_ref2_q[i+1,j_ref,k+1],
                         ni3d_ref2_q[i+1,j_ref,k+2],
                         ni3d_ref2_q[i,j_ref,k+2]
                        ]
                elcon_ref2 = vstack((elcon_ref2,nodes))    
elnums = (n.arange(len(elcon_ref2)) + n.max(elnums))[:,None]+1
elcon_ref2 = hstack((elcon_ref2,elnums)) 
print('Connected on Ref2')

elcon = n.vstack((elcon_fine,elcon_med,elcon_cors,elcon_ref1_q,elcon_ref1_r,elcon_ref2,)).astype(int)
# save it, put the el numbers out front
# Reorder so that connectivity processes in a different direction to get good mesh
elcon = elcon[:,[-1,3,2,1,0,7,6,5,4]]
#n.savetxt('Elements.dat',X = n.hstack((elcon[:,[-1]],elcon[:,:-1])) ,fmt='%.0f',delimiter=',')
#n.savetxt('elements_for_abaqus.dat',X = elcon ,fmt='%.0f',delimiter=',')
n.save('./ConstructionFiles/abaqus_elements.npy',elcon)


'''
def plotall():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d
    fig = p.figure(figsize=(15,12))
    ax = fig.add_subplot(111,projection='3d')
    p.tight_layout()
    NC = n.load('./ConstructionFiles/nc_all.npy')
    X = NC[:,0]*n.cos(NC[:,2])
    Y = NC[:,0]*n.sin(NC[:,2])
    rgn = (NC[:,2]>=angle-1*dq_med) | (NC[:,2]<=1*dq_med)
    nodes = vstack((X[rgn],Y[rgn],NC[rgn,1],NC[rgn,-1])).T
    ax.plot(nodes[:,0],nodes[:,1],nodes[:,2],'b.',alpha=0.5)
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)    
    N = elcon_ref1_r[n.arange(-6,6),:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nodes[:,-1],i), nodes, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,1],NC[loc,2],'ko',alpha=0.3)
            line, = ax.plot(NC[loc,0],NC[loc,1],NC[loc,2],'ro',alpha=1)
            p.pause(0.5)
            line.remove()

def fineconcheck():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    fig = p.figure(figsize=(15,12))
    ax = fig.add_subplot(111,projection='3d')
    p.tight_layout()
    nc_fine = n.load('./ConstructionFiles/nc_fine.npy')            
    rgn = (nc_fine[:,1] <=4*elht_fine) & (nc_fine[:,2] <=8*dq_fine)
    ax.plot(nc_fine[rgn,0],nc_fine[rgn,2],nc_fine[rgn,1],'b.',alpha=0.4)
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)    
    N = elcon_fine[:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nc_fine[:,-1],i), nc_fine, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.3)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.25)
            line.remove()

def medcheckcon():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    fig = p.figure(figsize=(15,12))
    ax = fig.add_subplot(111,projection='3d')
    p.tight_layout()    
    nc_med = n.load('./ConstructionFiles/nc_med.npy')            
    rgn = (nc_med[:,1] >=2*elht_med) & (nc_med[:,2] >=angle - 3*dq_med)
    ax.plot(nc_med[rgn,0],nc_med[rgn,2],nc_med[rgn,1],'b.',alpha=0.5)
    N = elcon_med[-6:,:-1]
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nc_med[:,-1],i), nc_med, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.3)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.5)
            line.remove()

def corscheckcon():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    fig = p.figure(figsize=(15,12))
    ax = fig.add_subplot(111,projection='3d')
    p.tight_layout()
    ax.axis('off')
    nc_cors = n.load('./ConstructionFiles/nc_cors.npy')            
    rgn = (nc_cors[:,1] >=5*elht_cors) & (nc_cors[:,2] >=angle - 3*dq_cors)
    rgn = nc_cors[:,2] < angle
    ax.plot(nc_cors[rgn,0],nc_cors[rgn,2],nc_cors[rgn,1],'b.',alpha=0.3)
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)    
    N = elcon_cors[:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nc_cors[:,-1],i), nc_cors, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.05)
            line.remove()

def ref1zcehckcon():
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    p.close('all')
    nc_all = n.load('./ConstructionFiles/nc_all.npy')    
    nc_fine = n.load('./ConstructionFiles/nc_fine.npy')            
    nc_ref1_r = n.load('./ConstructionFiles/nc_ref1_r.npy')            
    nc_ref1_mid = n.load('./ConstructionFiles/nc_ref1_mid.npy')            
    nc_ref1_q = n.load('./ConstructionFiles/nc_ref1_q.npy')      
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    p.tight_layout()
    rgn = ((ni_fine[:,1] == n.max(ni_fine[:,1])) &
            (ni_fine[:,2] <= 9) &
            (ni_fine[:,0] <= 6))
    qmax = n.max(nc_fine[rgn,2])
    rmax = n.max(nc_fine[rgn,0])
    ax.plot(nc_fine[rgn,0],nc_fine[rgn,2],nc_fine[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_q[:,2] <= qmax)&
            (nc_ref1_q[:,0] <= rmax))
    ax.plot(nc_ref1_q[rgn,0],nc_ref1_q[rgn,2],nc_ref1_q[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_mid[:,2] <= qmax)&
            (nc_ref1_mid[:,0] <= rmax))
    ax.plot(nc_ref1_mid[rgn,0],nc_ref1_mid[rgn,2],nc_ref1_mid[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_r[:,2] <= qmax)&
            (nc_ref1_r[:,0] <= rmax))
    ax.plot(nc_ref1_r[rgn,0],nc_ref1_r[rgn,2],nc_ref1_r[rgn,1],'b.',alpha=0.5)
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)    
    N = elcon_ref1_q[4:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nc_all[:,-1],i), nc_all, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(0.5)
            line.remove()
    
def ref1zcehckcon_end(ptime=0.25):
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    p.close('all')
    nc_all = n.load('./ConstructionFiles/nc_all.npy')    
    nc_fine = n.load('./ConstructionFiles/nc_fine.npy')            
    nc_ref1_r = n.load('./ConstructionFiles/nc_ref1_r.npy')            
    nc_ref1_mid = n.load('./ConstructionFiles/nc_ref1_mid.npy')            
    nc_ref1_q = n.load('./ConstructionFiles/nc_ref1_q.npy')      
    dq = n.diff(nc_fine[:,2]).max()
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    p.tight_layout()    
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)
    rgn = ((ni_fine[:,1] == n.max(ni_fine[:,1])))   # Highest z-coord fine nodes
    rgn = rgn & (ni_fine[:,2]>=(ni_fine[:,2].max()-10))
    qmin = n.min(nc_fine[rgn,2])
    rmax = n.max(nc_fine[rgn,0])
    ax.plot(nc_fine[rgn,0],nc_fine[rgn,2],nc_fine[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_q[:,2] >= qmin)&
            (nc_ref1_q[:,0] <= rmax))
    ax.plot(nc_ref1_q[rgn,0],nc_ref1_q[rgn,2],nc_ref1_q[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_mid[:,2] >= qmin)&
            (nc_ref1_mid[:,0] <= rmax))
    ax.plot(nc_ref1_mid[rgn,0],nc_ref1_mid[rgn,2],nc_ref1_mid[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_r[:,2] >= qmin)&
            (nc_ref1_r[:,0] <= rmax))
    ax.plot(nc_ref1_r[rgn,0],nc_ref1_r[rgn,2],nc_ref1_r[rgn,1],'b.',alpha=0.5)
    N = elcon_ref1_q[-4*num_el_fine_r:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nc_all[:,-1],i), nc_all, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(ptime)
            line.remove()

def ref1thcheckcon(ptime=0.25):
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    p.close('all')
    nc_all = n.load('./ConstructionFiles/nc_all.npy')    
    nc_fine = n.load('./ConstructionFiles/nc_fine.npy')            
    nc_ref1_r = n.load('./ConstructionFiles/nc_ref1_r.npy')            
    nc_ref1_mid = n.load('./ConstructionFiles/nc_ref1_mid.npy')            
    nc_ref1_q = n.load('./ConstructionFiles/nc_ref1_q.npy')     
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    p.tight_layout()    
    rgn = ((ni_fine[:,1] == n.max(ni_fine[:,1])) &
            (ni_fine[:,2] <= 9) &
            (ni_fine[:,0] <= 6))
    qmax = n.max(nc_fine[rgn,2])
    rmax = n.max(nc_fine[rgn,0])
    #ax.plot(nc_fine[rgn,0],nc_fine[rgn,2],nc_fine[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_q[:,2] <= qmax)&
            (nc_ref1_q[:,0] <= rmax))
    ax.plot(nc_ref1_q[rgn,0],nc_ref1_q[rgn,2],nc_ref1_q[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_mid[:,2] <= qmax)&
            (nc_ref1_mid[:,0] <= rmax))
    ax.plot(nc_ref1_mid[rgn,0],nc_ref1_mid[rgn,2],nc_ref1_mid[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref1_r[:,2] <= qmax)&
            (nc_ref1_r[:,0] <= rmax))
    ax.plot(nc_ref1_r[rgn,0],nc_ref1_r[rgn,2],nc_ref1_r[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_med[:,2] <= qmax)&
            (nc_med[:,0] <= rmax)&
            (nc_med[:,1] == n.min(nc_med[:,1])))
    ax.plot(nc_med[rgn,0],nc_med[rgn,2],nc_med[rgn,1],'b.',alpha=0.5)    
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)    
    N = elcon_ref1_q[-8:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nc_all[:,-1],i), nc_all, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(ptime)
            line.remove()
    
def ref2checkcon(ptime=0.25):
    import matplotlib.pyplot as p
    from mpl_toolkits import mplot3d    
    p.close('all')
    nc_all = n.load('./ConstructionFiles/nc_all.npy')    
    ni_med = n.load('./ConstructionFiles/ni_med.npy')
    nc_cors = n.load('./ConstructionFiles/nc_cors.npy')                
    nc_ref2_q = n.load('./ConstructionFiles/nc_ref2_q.npy')                
    fig = p.figure()
    ax = fig.add_subplot(111,projection='3d')
    p.tight_layout()    
    rgn = ((ni_med[:,1] == n.max(ni_med[:,1])))
    qmax = n.max(nc_med[rgn,2])
    rmax = n.max(nc_med[rgn,0])
    ax.plot(nc_med[rgn,0],nc_med[rgn,2],nc_med[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_ref2_q[:,2] <= qmax)&
            (nc_ref2_q[:,0] <= rmax))
    ax.plot(nc_ref2_q[rgn,0],nc_ref2_q[rgn,2],nc_ref2_q[rgn,1],'b.',alpha=0.5)
    rgn = ((nc_cors[:,2] <= qmax)&
            (nc_cors[:,0] <= rmax)&
            (ni_cors[:,1] == 0))
    ax.plot(nc_cors[rgn,0],nc_cors[rgn,2],nc_cors[rgn,1],'b.',alpha=0.5)
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)    
    N = elcon_ref2[-4*num_el_med_r:,:-1]
    for j,i in enumerate(N):
        NC = n.compress(n.in1d(nc_all[:,-1],i), nc_all, axis=0)
        for k,z in enumerate(i):
            loc = n.nonzero(NC[:,-1] == z)[0]
            ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ko',alpha=0.15)
            line, = ax.plot(NC[loc,0],NC[loc,2],NC[loc,1],'ro',alpha=1)
            p.pause(ptime)
            line.remove()
'''
