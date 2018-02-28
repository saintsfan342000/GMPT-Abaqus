''' Read ODB
job_name:   ODB name without odb extension
inst_name:  name of instance
step_name:  name of load step
ASSEMBLY.NS_RPTOP : The top reference point to which force and torque are applied
ASSEMBLY.NS_RPBOT : The bottom ref. pt. which is fixed
INSTANCE.NS_DISPROT_LO : The node just on the cusp of the the test section
INSTANCE.NS_DISPROT_HI : The node one above LO (in case I want to interpolate btwn them)
INSTANCE.ES_Z : The line of elements at theta=0
'''

import odbAccess as O
import numpy as np
import os
from sys import argv
pi = np.pi

job = argv[1]
half = True

t = 0.0500
R = 0.8391

if not os.path.isfile(job + '.odb'):
    raise ValueError('The specified job name "%s" does not exist.'%(job))

nset_rp_top = 'NS_RPTOP'
nset_rp_bot = 'NS_RPBOT'
nset_dr_lo = 'NS_DISPROT_LO'
nset_dr_back = 'NS_DISPROT_LO_BACK'
nset_radcont = 'NS_RADIALCONTRACTION'
elset_th_front = 'ES_THICKNESS'
elset_th_back = 'ES_THICKNESS_BACK'
elset_th_side = 'ES_THICKNESS_SIDE'


h_odb = O.openOdb(job + '.odb',readOnly=True)
h_inst = h_odb.rootAssembly.instances[ h_odb.rootAssembly.instances.keys()[0] ]
h_step = h_odb.steps[ h_odb.steps.keys()[0] ]
h_All_Frames = h_step.frames
num_incs = len(h_All_Frames)
frameNos = np.empty(num_incs, dtype=int)
for k,frame in enumerate(h_All_Frames):
    frameNos[k] = frame.frameId

h_nset_rp_top = h_inst.nodeSets[nset_rp_top]
h_nset_dr_lo = h_inst.nodeSets[nset_dr_lo]
h_nset_dr_back = h_inst.nodeSets[nset_dr_back]
h_nset_urprof = h_inst.nodeSets[nset_radcont]
h_elset_th = h_inst.elementSets[elset_th_front]
h_elset_th_b = h_inst.elementSets[elset_th_back]
h_elset_th_s = h_inst.elementSets[elset_th_side]
numel_th = len(h_elset_th.elements)
# Need to find the historyRegion key that contains the 
# NS_PCAV node to get CVOL and PCAV.In a riks analysis, 
# the LPF is also written in a history region that is
# the Assembly. We don't want that one.
for k,i in enumerate( h_step.historyRegions.keys() ):
    if i.rfind('Node') != -1:
        h_histrgn = h_step.historyRegions[i]
        break

F = np.empty( (num_incs) )
P = np.empty( (num_incs) )
V = np.empty( (num_incs) )
d_lo = np.empty( (num_incs) )
d_back = np.empty( (num_incs) )

sts = np.empty((num_incs,3))
stn = np.empty((num_incs,3))

# Grab undef coords of dr_lo and hi
Lg_lo = h_nset_dr_lo.nodes[0].coordinates[2]
Lg_back = h_nset_dr_back.nodes[0].coordinates[2]

# Some special handling for ur_prof since the jth value of U subset does not correspond to jth node in the set :(
numnode_prof = len(h_nset_urprof.nodes)
urnodes = np.array([i.label for i in h_nset_urprof.nodes])
ur_prof = np.empty((num_incs+2, numnode_prof))

for i in range(num_incs):
    F[i] = h_All_Frames[i].fieldOutputs['CF'].getSubset(region=h_nset_rp_top).values[0].data[2]
    # P[i] = h_All_Frames[i].fieldOutputs['P'].values[0].data
    d_lo[i] = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_dr_lo).values[0].data[2]
    d_back[i] = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_dr_back).values[0].data[2]
    tempsts, tempstn = 0, 0
    stresses = h_All_Frames[i].fieldOutputs['S'].getSubset(region=h_elset_th).values
    strains =  h_All_Frames[i].fieldOutputs['LE'].getSubset(region=h_elset_th).values
    for j in range(numel_th):
         tempsts += stresses[j].data[:3]
         tempstn += strains[j].data[:3]
    sts[i] = tempsts/numel_th
    stn[i] = tempstn/numel_th
    
    # For some reason, the jth node in the set didn't correspond to the .values[j].nodeLabel
    # So I have to do this ridiculous thing where I create a subset for each node in the set
    uvals = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_urprof)
    for j in range(numnode_prof):
        h_node = h_nset_urprof.nodes[j]
        if i == 0:
            ur_prof[0,j] = h_node.coordinates[2]
            ur_prof[1,j] = h_node.coordinates[0]
        nodeUvals = uvals.getSubset(region=h_node).values[0]
        ur_prof[i+2,j] = nodeUvals.data[0]
        if nodeUvals.nodeLabel != h_node.label:
            print "Warning...you're getting node mistmatch in the profile."
            
    
    #for j in range(numnode_prof):
    #    print h_nset_urprof.nodes[j].label, h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_urprof).values[j].nodeLabel


# Now get history data. 
z = 0
for k,i in enumerate(h_histrgn.historyOutputs['PCAV'].data):
    if k in frameNos:
        P[z] = i[1]
        z+=1
z = 0
for k,i in enumerate(h_histrgn.historyOutputs['CVOL'].data):
    if k in frameNos:
        V[z] = i[1]
        z+=1

h_odb.close()

# Convert disp to d/Lg
d_lo, d_back = d_lo/Lg_lo, d_back/Lg_back

if half == True:
    # Since half model, multiple F by 2
    F*=2

sig_x = F/(2*pi*R*t) + P*R/(2*t)
sig_q = P*R/t

def headerline(fname, hl):
    fid = open(fname, 'r')
    data = fid.read()
    fid.close()
    del fid
    fid = open(fname, 'w')
    fid.write(hl)
    if hl[-1] != '\n':
        fid.write('\n')
    fid.write(data)
    fid.close()

# Save
fname = '%s_results.dat'%(job)

np.savetxt(fname, X = np.vstack((F, P, sig_x, sig_q, d_lo, d_back, sts.T, stn.T, V)).T, fmt='%.6f', delimiter=', ')
hl = '#[0] Force (kip), [1]Pressure (ksi), [2]NomAxSts, [3]NomHoopSts,' 
hl += ' [4]d/Lg lo, [5]d/Lg Back, [6,7,8]S11,22,33, [9,10,11]LE11,22,33, [12]Vol'
headerline(fname, hl)

ur_prof = ur_prof[:,ur_prof[0].argsort()]
fname = '%s_UR.dat'%(job)
np.savetxt(fname, X=ur_prof, fmt='%.8f', delimiter=',')
hl = '#1st line: Undef Z-coord of nodes\n'
hl += '#2nd line: Undef r-coord of nodes\n'
hl += '#3+nth line: nth incremend Ur of node in column'
headerline(fname,hl)

