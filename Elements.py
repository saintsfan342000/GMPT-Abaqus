import numpy as n
from numpy import (sqrt, linspace, asarray,
                    empty, hstack , vstack, pi)

'''
r, [:,0] = thru thickness/radial
z, [:,1] = Axial coord
q, [:,2] = Angle coordinate theta (radians)

IMPORTANT
In all cases, I connect as follows:
   See the ConstructionFiles/NodeIndexOrder.png
   The nodes (0,0,0), (0,1,0), (1,1,0), (1,0,0)
   and nodes (0,0,1), (0,1,1), (1,1,1), (1,0,1)
   form an element.  In all cases I connect in this order:
   (0,0,1)-->(1,0,1)-->(1,1,1)-->(0,1,1)-->(0,0,0)-->(1,0,0)-->(1,1,0)-->(0,1,0)
It turns out that this is not an acceptable order for abaqus.  So after all the
connectivity loops, I then reorder the elcon.

'''

fullring = False

# Load up the 3d arrays
# They all come in as ints 
nodelist = ['ni3d_fine', 'ni3d_med', 'ni3d_cors',
            'ni3d_lowmed', 'ni3d_ref1_mid', 'ni3d_ref1_r',
            'ni3d_ref1_q', 'ni3d_ref2_q', 'ni3d_reff_z']
for k,i in enumerate(nodelist):
    exec('{} = n.load("./ConstructionFiles/{}.npy")'.format(i,i))


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

q_lowmed = n.nonzero( ni3d_lowmed[0,0,:]!=ni3d_lowmed[0,0,0])[0][0]
num_el_lowmed_r = ni3d_lowmed.shape[0] - 1
num_el_lowmed_z = ni3d_lowmed.shape[1] - 1
# To account for the non-zero start
num_el_lowmed_q = ni3d_lowmed[:,:,q_lowmed:].shape[2] - 1
num_el_lowmed_tot = num_el_lowmed_r*num_el_lowmed_z*num_el_lowmed_q

########################################################    
################## Connnect on Fine ##################
########################################################
# Notice, the inner-most loop is "for i in num_el_fine_r",
# meaning I will connect thru-thickness first, then in q/circ, then in z/axial
elcon_fine = empty((num_el_fine_tot,9), dtype=int)
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
            elcon_fine[row,:-1] = elcon_temp
            row+=1
elnums = n.arange(elcon_fine.shape[0])+1
elcon_fine[:,-1] = elnums
print('Connected on Fine') 
########################################################    
################## Connnect on Med ##################
########################################################
elcon_med = empty((num_el_med_tot,9), dtype=int)
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
            elcon_med[row, :-1] = elcon_temp
            row+=1
elnums = n.arange(elcon_med.shape[0]) + n.max(elnums) + 1
elcon_med[:,-1] = elnums
print('Connected on Med')
########################################################    
################## Connnect on Coarse ##################
########################################################
elcon_cors = n.empty((num_el_cors_tot,9), dtype=int)
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
            elcon_cors[row, :-1] = elcon_temp
            row+=1
elnums = n.arange(elcon_cors.shape[0]) + n.max(elnums) + 1
elcon_cors[:,-1] = elnums 
print('Connected on Coarse')

########################################################    
################## Connnect on Lowmed ##################
################## And lowmed to ref1_mid ##################
########################################################
elcon_lowmed = n.empty((num_el_lowmed_tot,9), dtype=int)
row = 0
for j in (n.arange(num_el_lowmed_z)):
    for k in range(q_lowmed, q_lowmed+num_el_lowmed_q):
        for i in range(num_el_lowmed_r):
            index = n.array([[i,j+1,k], 
                            [i+1,j+1,k],
                            [i+1,j+1,k+1],
                            [i,j+1,k+1],
                            [i,j,k],
                            [i+1,j,k],
                            [i+1,j,k+1],
                            [i,j,k+1]]
                           )
            elcon_temp = [0,0,0,0,0,0,0,0]
            for u,v in enumerate(index):
                elcon_temp[u] = ni3d_lowmed[v[0],v[1],v[2]]
            elcon_lowmed[row, :-1] = elcon_temp
            row+=1
elnums = n.arange(elcon_lowmed.shape[0]) + elnums.max() + 1
elcon_lowmed[:,-1] = elnums

z_lowmed = num_el_lowmed_z
z_ref1_mid = 0
elcon_lowmed_mid = n.empty((num_el_lowmed_r*num_el_lowmed_q,9), dtype=int)
row=0
for k in range(q_lowmed, q_lowmed+num_el_lowmed_q):
    for i in range(num_el_lowmed_r):
        nodes = [ni3d_ref1_mid[i,z_ref1_mid,k],
                 ni3d_ref1_mid[i+1,z_ref1_mid,k],
                 ni3d_ref1_mid[i+1,z_ref1_mid,k+1],
                 ni3d_ref1_mid[i,z_ref1_mid,k+1],
                 ni3d_lowmed[i,z_lowmed,k],
                 ni3d_lowmed[i+1,z_lowmed,k],
                 ni3d_lowmed[i+1,z_lowmed,k+1],
                 ni3d_lowmed[i,z_lowmed,k+1]
                ]
        elcon_lowmed_mid[row, :-1] = nodes
        row+=1
elnums = n.arange(elcon_lowmed_mid.shape[0]) + elnums.max() + 1
elcon_lowmed_mid[:,-1] = elnums
print('Connected Lowmed and Lowmed to Mid')

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
elcon_ref1_q = n.empty((num_el_fine_r*num_el_fine_q*4//3,9), dtype=int)
row = 0
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
            elcon_ref1_q[row,:-1] = nodes
            row+=1
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
            elcon_ref1_q[row,:-1] = nodes
            row+=1
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
            elcon_ref1_q[row,:-1] = nodes
            row+=1
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
            elcon_ref1_q[row,:-1] = nodes    
            row+=1

elnums = n.arange(elcon_ref1_q.shape[0]) + n.max(elnums) + 1
elcon_ref1_q[:,-1] = elnums

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
elcon_ref1_r = n.empty((num_el_med_r*4*num_el_med_q,9), dtype=int)
row=0
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
            elcon_ref1_r[row, :-1] = nodes
            row+=1
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
            elcon_ref1_r[row, :-1] = nodes
            row+=1
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
            elcon_ref1_r[row, :-1] = nodes                
            row+=1
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
            elcon_ref1_r[row, :-1] = nodes
            row+=1
            
elnums = n.arange(elcon_ref1_r.shape[0]) + n.max(elnums) + 1
elcon_ref1_r[:,-1] = elnums  
print('Connected on Ref1')


#######################################################
########## Connect fine to lowmed and reff_z, and reff_z to low med
#######################################################
q_fine = ni3d_fine[0,0,:].max() - 1
q_reff = 0
q_lowmed = n.nonzero( ni3d_lowmed[0,0,:]!=ni3d_lowmed[0,0,0])[0][0]
elcon_reff = n.empty((4*num_el_fine_r*num_el_fine_z//3 , 9), dtype=int)
row = 0
for j in (n.arange(num_el_fine_z)):
    if j%3 == 0:
        for i in range(num_el_fine_r):
            # Bottom right-trapezoid
            nodes = [ni3d_fine[i,j+1,q_fine],
                     ni3d_fine[i+1,j+1,q_fine],
                     ni3d_reff_z[i+1,j+1,q_reff],
                     ni3d_reff_z[i,j+1,q_reff],
                     ni3d_fine[i,j,q_fine],
                     ni3d_fine[i+1,j,q_fine],
                     ni3d_lowmed[i+1,j//3,q_lowmed], 
                     ni3d_lowmed[i,j//3,q_lowmed]
                    ]
            elcon_reff[row, :-1] = nodes
            row+=1
    elif j%3 == 1:
        for i in range(num_el_fine_r):
            # Center square
            nodes = [ni3d_fine[i,j+1,q_fine],
                     ni3d_fine[i+1,j+1,q_fine],
                     ni3d_reff_z[i+1,j+1,q_reff],
                     ni3d_reff_z[i,j+1,q_reff],
                     ni3d_fine[i,j,q_fine],
                     ni3d_fine[i+1,j,q_fine],
                     ni3d_reff_z[i+1,j,q_reff],
                     ni3d_reff_z[i,j,q_reff]
                    ]
            elcon_reff[row, :-1] = nodes
            row+=1
    elif j%3 == 2:
        for i in range(num_el_fine_r):
            # Top right-trapezoid
            nodes = [ni3d_fine[i,j+1,q_fine],
                     ni3d_fine[i+1,j+1,q_fine],
                     ni3d_lowmed[i+1,(j+1)//3,q_lowmed],
                     ni3d_lowmed[i,(j+1)//3,q_lowmed],
                     ni3d_fine[i,j,q_fine],
                     ni3d_fine[i+1,j,q_fine],
                     ni3d_reff_z[i+1,j,q_reff],
                     ni3d_reff_z[i,j,q_reff]
                    ]
            elcon_reff[row, :-1] = nodes
            row+=1
    # In a separate loop, connect the large iscoceles trapezoid 
    # that connects reff_z to lowmed
    if (j%3 == 0):
        for i in range(num_el_fine_r):
            nodes = [ni3d_reff_z[i,j+2,q_reff],
                     ni3d_reff_z[i+1,j+2,q_reff],
                     ni3d_lowmed[i+1,j//3+1,q_lowmed],
                     ni3d_lowmed[i,j//3+1,q_lowmed],
                     ni3d_reff_z[i,j+1,q_reff],
                     ni3d_reff_z[i+1,j+1,q_reff],
                     ni3d_lowmed[i+1,j//3,q_lowmed],
                     ni3d_lowmed[i,j//3,q_lowmed]
                    ]
            elcon_reff[row, :-1] = nodes
            row+=1
elnums = n.arange(elcon_reff.shape[0]) + n.max(elnums) + 1
elcon_reff[:,-1] = elnums  


########################################################    
##### Connect corner fine to ref1_mid and lowmed  ######
########################################################
z_fine = num_el_fine_z
q_fine = num_el_fine_q
z_lowmed = num_el_lowmed_z
elcon_corner = n.empty((num_el_fine_r,9), dtype=int)
for z,i in enumerate(range(num_el_fine_r)):
    nodes = [ni3d_ref1_mid[i,0,q_lowmed-1],
             ni3d_ref1_mid[i+1,0,q_lowmed-1],
             ni3d_ref1_mid[i+1,0,q_lowmed],
             ni3d_ref1_mid[i,0,q_lowmed],
             ni3d_fine[i,z_fine,q_fine],
             ni3d_fine[i+1,z_fine,q_fine],
             ni3d_lowmed[i+1,z_lowmed,q_lowmed],
             ni3d_lowmed[i,z_lowmed,q_lowmed]
            ]
    elcon_corner[z,:-1] = nodes
elnums = n.arange(elcon_corner.shape[0]) + elnums.max() + 1
elcon_corner[:,-1] = elnums
print('Connected Fine to Lowmed')

########################################################    
########## Connect med to ref2, ref2 to cors
########################################################
j_med = ni3d_med.shape[1] - 1
j_cors = 0
j_ref = 0
elcon_ref2 = empty((num_el_cors_r*num_el_cors_q*4,9), dtype=int)
row = 0
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
            elcon_ref2[:,:-1] = nodes
            row+=1
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
            elcon_ref2[:,:-1] = nodes
            row+=1
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
                elcon_ref2[:,:-1] = nodes
            row+=1
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
            elcon_ref2[:,:-1] = nodes
            row+=1
elnums = n.arange(elcon_ref2.shape[0]) + n.max(elnums) + 1
elcon_ref2[:,-1] = elnums 
print('Connected on Ref2')

elcon = n.vstack((elcon_fine, elcon_med, elcon_cors,
                  elcon_lowmed, elcon_lowmed_mid, elcon_ref1_q,
                  elcon_ref1_r, elcon_reff, elcon_corner, 
                  elcon_ref2,))
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
    elht_fine = n.diff(n.unique(nc_fine[:,1]))[0]
    dq_fine = n.diff(n.unique(nc_fine[:,2]))[0]
    rgn = (nc_fine[:,1] <=4*elht_fine) & (nc_fine[:,2] <=8*dq_fine)
    ax.plot(nc_fine[rgn,0],nc_fine[rgn,2],nc_fine[rgn,1],'b.',alpha=0.4)
    ax.set_xlabel('X',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)    
    ax.view_init(5,0)
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
    ax.view_init(5,0)
    ax.set_xlabel('R',labelpad=30)
    ax.set_ylabel('$\\theta$',labelpad=30)
    ax.set_zlabel('Z',labelpad=30)    
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
    N = elcon_ref1_q[:,:-1]
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
