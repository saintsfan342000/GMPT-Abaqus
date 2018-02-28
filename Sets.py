import numpy as n
from numpy import (linspace, asarray, in1d, empty,
                    hstack , vstack, pi, compress)
from sys import argv

'''
Generates the node sets and element sets.
'''

fullring = False

cdt = float(argv[1])
Lg = 2

nodelist = ['nc_cors', 'nc_fine', 'nc_med',
            'nc_ref1_mid', 'nc_ref1_r', 'nc_ref1_q',
            'nc_ref2_q', 'ni_cors', 'ni_fine',
            'ni_med', 'ni_ref1_mid', 'ni_ref1_r',
            'ni_ref1_q', 'ni_ref2_q', 'nc_lowmed',
            'ni_lowmed', 'nc_reff_z', 'ni_reff_z',
            'nc_all', 'ni_all']

for k,name in enumerate(nodelist):
    exec('{} = n.load("./ConstructionFiles/{}.npy")'.format(name,name))

E = n.load('./ConstructionFiles/abaqus_elements.npy')

fid = open('./ConstructionFiles/abaqus_sets.txt','w')

################
## Node sets ###
################

# nset TopSurface
# A set of all nodes on the top surface of the model
fid.write('*nset, nset=NS_TOPSURFACE\n')
rng = (ni_cors[:,1] == ni_cors[:,1].max())
nodenums = compress(rng, ni_cors[:,3])
for i,no in enumerate(nodenums):
    if ((i+1)%16 == 0) or (i == len(nodenums)-1):
        fid.write('{:.0f}\n'.format(no))
    else:
        fid.write('{:.0f}, '.format(no))

# nset BottomSurface
# A set of all nodes on the bottom surface of the model
fid.write('*nset, nset=NS_BOTTOMSURFACE\n')
rng = (ni_fine[:,1] == 0 )
nodenums = compress(rng, ni_fine[:,3])
for i,no in enumerate(nodenums):
    if ((i+1)%16 == 0) or (i == len(nodenums)-1):
        fid.write('{:.0f}\n'.format(no))
    else:
        fid.write('{:.0f}, '.format(no))
rng = (ni_lowmed[:,1] == 0 )
nodenums = compress(rng, ni_lowmed[:,3])
for i,no in enumerate(nodenums):
    if ((i+1)%16 == 0) or (i == len(nodenums)-1):
        fid.write('{:.0f}\n'.format(no))
    else:
        fid.write('{:.0f}, '.format(no))

# nset Disp-Rot node
# A single node to calculate nominal disp and rot
# Right at the cusp of the radius
# index z = 0, index theta = 0, index_r = max
fid.write('*nset, nset=NS_DISPROT_HI\n')
rng = ( (ni_cors[:,1] == 0 ) & (ni_cors[:,2] == 0 ) &
        (ni_cors[:,0] == ni_cors[:,0].max()) )
nodenums = compress(rng, ni_cors[:,3])
if len(nodenums) != 1:
    raise ValueError('Seeking a single node for NS_DISPROT_LO, but len(nodenums)!=1')
else:
    fid.write('{}\n'.format(nodenums[0]))

# NS_DISP_LO:  A node at 0.5" above symplane on the front
fid.write('*nset, nset=NS_DISPROT_LO\n')
po = n.array([1.9685/2, 0.5, 0]) # r, z, q
loc = n.linalg.norm(nc_all[:,:-1] - po, axis=1).argmin()
nodenum = nc_all[loc, -1]
fid.write('{:.0f}\n'.format(nodenum))

# NS_DISP_LO_BACK:  A node at 0.5" above symplane on the back
fid.write('*nset, nset=NS_DISPROT_LO_BACK\n')
po = n.array([1.9685/2, 0.5, pi]) # r, z, q
loc = n.linalg.norm(nc_all[:,:-1] - po, axis=1).argmin()
nodenum = nc_all[loc, -1]
fid.write('{:.0f}\n'.format(nodenum))

#nset radial contraction
# A line of nodes running up along the OD
fid.write('*nset, nset=NS_RADIALCONTRACTION\n')
# index theta = 0, index_r = max. From fine, ref1_mid, med, cors
rng = (ni_fine[:,2] == 0) & (ni_fine[:,0] == ni_fine[:,0].max())
nodenums = compress(rng, ni_fine[:,3])
rng = (ni_ref1_mid[:,0] == ni_ref1_mid[:,0].max()) & (ni_ref1_mid[:,2] == 0)
nodenums = hstack(( nodenums, compress(rng, ni_ref1_mid[:,3]) ))
rng = (ni_med[:,2] == 0) & (ni_med[:,0] == ni_med[:,0].max())
nodenums = hstack(( nodenums, compress(rng, ni_med[:,3]) ))
rng = (ni_cors[:,2] == 0) & (ni_cors[:,0] == ni_cors[:,0].max())
nodenums = hstack(( nodenums, compress(rng, ni_cors[:,3]) ))
for i,no in enumerate(nodenums):
    if ((i+1)%16 == 0) or (i == len(nodenums)-1):
        fid.write('{:.0f}\n'.format(no))
    else:
        fid.write('{:.0f}, '.format(no))
# NS_AxialSymm
# A set of nodes on the two edges if we've done half-ring
if not fullring:
    fid.write('*nset, nset=NS_AXIALSYMM\n')
    rng = (nc_all[:,2]==0) | (nc_all[:,2] == nc_all[:,2].max())
    nodenums = nc_all[rng,3]
    for i,no in enumerate(nodenums):
        if ((i+1)%16 == 0) or (i == len(nodenums)-1):
            fid.write('{:.0f}\n'.format(no))
        else:
            fid.write('{:.0f}, '.format(no))
if cdt > 0:
    fid.write('*nset, nset=NS_CIRCIMPERF\n')
    # Get the ID and thickness
    ID = n.min(nc_med[:,0])
    tg = n.max(nc_fine[:,0]) - ID 
    rng = (ni_all[:,0]==0) & (nc_all[:,2] <= 2*tg/ID) & (nc_all[:,1] <= ID+tg/2)
    nodenums = nc_all[rng,3]
    for i,no in enumerate(nodenums):
        if ((i+1)%16 == 0) or (i == len(nodenums)-1):
            fid.write('{:.0f}\n'.format(no))
        else:
            fid.write('{:.0f}, '.format(no))

################
# Element sets #
################

# Elset_q
fid.write('*elset, elset=ES_Z\n')
# A line of elements on the OD running up the test section
# Place it on the thinnest wall-thickness area (y = 0, x = OD/2)
# We'll need to grab from fine, all ref1s, and med
# This could need modification if I change the coordinate of ref1
# Fine
rng = (ni_fine[:,0] == ni_fine[:,0].max()) & (ni_fine[:,2] == 0)
nodenums = compress(rng, ni_fine[:,3])
# Do it again, specifying one column of nodes over, so that I only get one colum of elements
nodenums2 = nodenums + 1
rng = (in1d(E[:,3],nodenums)) & (in1d(E[:,2],nodenums2))   # VERY VERY dependent on node-ordering!!
elnums = E[ rng, 0]
# ref1
rng = (ni_ref1_mid[:,0] == ni_ref1_mid[:,0].max()) & (ni_ref1_mid[:,2] == 0)
nodenums = compress(rng, ni_ref1_mid[:,3])
#nodenums = hstack(( nodenums, compress(rng, ni_ref1_mid[:,3]) ))
rng = (ni_ref1_q[:,0] == ni_ref1_q[:,0].max()) & (ni_ref1_q[:,2] == 1)
#nodenums = hstack(( nodenums, compress(rng, ni_ref1_mid[:,3]) ))
nodenums2 = compress(rng, ni_ref1_q[:,3])
rng = (in1d(E[:,3],nodenums)) & (in1d(E[:,2],nodenums2))
elnums = hstack(( elnums, E[ rng, 0] ))
# Med
rng = (ni_med[:,0] == ni_med[:,0].max()) & (ni_med[:,2] == 0) & (nc_med[:,1]<=Lg)
nodenums = compress(rng, ni_med[:,3])
nodenums2 = nodenums + 1
rng = (in1d(E[:,3],nodenums)) & (in1d(E[:,2],nodenums2))
elnums = hstack(( elnums, E[ rng, 0] ))
for i,el in enumerate(elnums):
    if ((i+1)%16 == 0) or (i == len(elnums)-1):
        fid.write('{:.0f}\n'.format(el))
    else:
        fid.write('{:.0f}, '.format(el))

# Elset_thickness
fid.write('*elset, elset=ES_THICKNESS\n')
# A line of elements on the sym-plane running thru-thickness
# Places along thinnest wall-thickness area (y = 0, x>0)
rng = (ni_fine[:,1] == 1) & (ni_fine[:,2] == 0) # == 1 b/c top nodes of element
nodenums = compress(rng, ni_fine[:,3])
nodenums2 = nodenums + 1
rng = (in1d(E[:,3],nodenums)) & (in1d(E[:,2],nodenums2))
elnums = E[rng, 0]
for i,el in enumerate(elnums):
    if ((i+1)%16 == 0) or (i == len(elnums)-1):
        fid.write('{:.0f}\n'.format(el))
    else:
        fid.write('{:.0f}, '.format(el))

# Elset_thickness_back
fid.write('*elset, elset=ES_THICKNESS_BACK\n')
# A line of elements on the sym-plane running thru-thickness
# Places along THICKNESS wall-thickness area (y = 0, x < 0)
# z-index = 1, and theta-coord closest to pi
rng = ((ni_lowmed[:,1] == 1) &
        (n.abs(nc_lowmed[:,2]-pi) == n.abs(nc_lowmed[:,2]-pi).min()) )
nodenums = compress(rng, ni_lowmed[:,3])
nodenums2 = nodenums - 1 # Minus to accommodate half model!
# Order of 2 and 3 is switched b/c nodenums2 is minus 1
rng = (in1d(E[:,2],nodenums)) & (in1d(E[:,3],nodenums2))
elnums = E[rng, 0]
for i,el in enumerate(elnums):
    if ((i+1)%16 == 0) or (i == len(elnums)-1):
        fid.write('{:.0f}\n'.format(el))
    else:
        fid.write('{:.0f}, '.format(el))

# Elset_thickness_side
fid.write('*elset, elset=ES_THICKNESS_SIDE\n')
# A line of elements on the sym-plane running thru-thickness
# Placed along  wall-thickness area (y = ymax, x = 0)
# z-index = 1, and theta-coord closest to pi/2
rng = ((ni_lowmed[:,1] == 1) &
        (n.abs(nc_lowmed[:,2]-pi/2) == n.abs(nc_lowmed[:,2]-pi/2).min()) )
nodenums = compress(rng, ni_lowmed[:,3])
nodenums2 = nodenums + 1
rng = (in1d(E[:,3],nodenums)) & (in1d(E[:,2],nodenums2))
elnums = E[rng, 0]
for i,el in enumerate(elnums):
    if ((i+1)%16 == 0) or (i == len(elnums)-1):
        fid.write('{:.0f}\n'.format(el))
    else:
        fid.write('{:.0f}, '.format(el))


# Elset_wholeid
fid.write('*elset, elset=ES_WHOLEID\n')
# Every element with a face on the ID
# Just take all nodes with thickness index == 0, and
# all elements that contain any such a node
rng = (ni_all[:,0] == 0)
nodenums = compress(rng, ni_all[:,3])
#rng = in1d(E[:,1:], nodenums)
rng = n.zeros_like(E.shape[0], dtype=bool)
for i in range(1,E.shape[1]):
    rng = rng | in1d(E[:,i], nodenums)
elnums = E[rng, 0]
# Now exclude the elements that connect ref1_r to bottom of med thru thickness, because they don't have a face on the ID
# Any element containing both a node with ni_ref1_r[:,0] == 1 and ni_ref1_r[:,0]==2
# An element must contain a node from each! of the following
omit_nodes1 = ni_ref1_r[ ni_ref1_r[:,0] == 1, 3 ] # Node numbers to exclude
omit_nodes2 = ni_ref1_r[ ni_ref1_r[:,0] == 2, 3 ] # Node numbers to exclude
# Now the elnums to keep
# Requires that element numbers increase sequentially, i.e. E[j,0] is el. no j+1
omit_els = n.array([n.any(in1d(omit_nodes1, E[el-1,1:])) &
                     n.any(in1d(omit_nodes2, E[el-1,1:]))
                    for el in elnums])
elnums = elnums[~omit_els]
for i,el in enumerate(elnums):
    if ((i+1)%16 == 0) or (i == len(elnums)-1):
        fid.write('{:.0f}\n'.format(el))
    else:
        fid.write('{:.0f} ,'.format(el))

fid.close()
