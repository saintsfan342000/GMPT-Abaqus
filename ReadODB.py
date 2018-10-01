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
from __future__ import division
import odbAccess as O
import numpy as np
import os
from sys import argv
pi = np.pi

job = argv[1]
if not os.path.isfile(job + '.odb'):
    raise ValueError('The specified job name "%s" does not exist.'%(job))
half = True

# Get the true R and t from ExptSummary
# For generality, recursively look for the ExpSumm.dat file
# If we don't find it, then use defaults
relpath = '.'
for i in range(10):
    if not os.path.exists('%s/ExptSummary.dat'%(relpath)):
        relpath = '../' + relpath
    else:
        key = np.genfromtxt('%s/ExptSummary.dat'%(relpath), delimiter=',')
        x = int(job.split('-')[-1])
        key = key[key[:,0] == x].ravel()
        R, t = key[6:8]
        break
else:
    t = 0.0500
    R = 0.8391
    print 'ExptSummart.dat not located. Using default R and t values'

nset_rp_top = 'NS_RPTOP'
nset_rp_bot = 'NS_RPBOT'
nset_dr_lo = 'NS_DISPROT_LO'
nset_dr_back = 'NS_DISPROT_LO_BACK'
nset_radcont = 'NS_RADIALCONTRACTION'
elset_th_front = 'ES_ANALZONE'
elset_th_back = 'ES_THICKNESS_BACK'
elset_th_side = 'ES_THICKNESS_SIDE'
elset_leprof = 'ES_LEPROF'

h_odb = O.openOdb(job + '.odb',readOnly=True)
h_inst = h_odb.rootAssembly.instances[ h_odb.rootAssembly.instances.keys()[0] ]
h_step = h_odb.steps[ h_odb.steps.keys()[0] ]
h_All_Frames = h_step.frames
num_incs = len(h_All_Frames)

frameNos = np.empty(num_incs, dtype=int)
#for k,frame in enumerate(h_All_Frames):
#    frameNos[k] = frame.frameId

h_nset_rp_top = h_inst.nodeSets[nset_rp_top]
h_nset_dr_lo = h_inst.nodeSets[nset_dr_lo]
h_nset_dr_back = h_inst.nodeSets[nset_dr_back]
h_nset_urprof = h_inst.nodeSets[nset_radcont]
h_elset_th = h_inst.elementSets[elset_th_front]
h_elset_th_b = h_inst.elementSets[elset_th_back]
h_elset_th_s = h_inst.elementSets[elset_th_side]
h_elset_leprof = h_inst.elementSets[elset_leprof]
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
peeq = np.empty((num_incs))
eqsts = np.empty((num_incs))

# Grab undef coords of dr_lo and hi
Lg_lo = h_nset_dr_lo.nodes[0].coordinates[2]
Lg_back = h_nset_dr_back.nodes[0].coordinates[2]

numnode_prof = len(h_nset_urprof.nodes)
ur_prof = np.empty((num_incs+2, numnode_prof))
numel_leprof = len(h_elset_leprof.elements)
le_prof = np.empty((numel_leprof,num_incs+1))

for i in range(num_incs):
    F[i] = h_All_Frames[i].fieldOutputs['CF'].getSubset(region=h_nset_rp_top).values[0].data[2]
    # P[i] = h_All_Frames[i].fieldOutputs['P'].values[0].data
    d_lo[i] = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_dr_lo).values[0].data[2]
    d_back[i] = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_dr_back).values[0].data[2]
    tempsts, tempstn, temppeeq, tempeqsts = 0, 0, 0 ,0
    stresses = h_All_Frames[i].fieldOutputs['S'].getSubset(region=h_elset_th).values
    strains =  h_All_Frames[i].fieldOutputs['LE'].getSubset(region=h_elset_th).values
    try:
        peeq_set =  h_All_Frames[i].fieldOutputs['SDV1'].getSubset(region=h_elset_th).values
        eqsts_set =  h_All_Frames[i].fieldOutputs['SDV2'].getSubset(region=h_elset_th).values
        for j in range(numel_th):
             tempsts += stresses[j].data[:3]
             tempstn += strains[j].data[:3]
             temppeeq += peeq_set[j].data
             tempeqsts = eqsts_set[j].data
    except KeyError:
        peeq_set =  h_All_Frames[i].fieldOutputs['PEEQ'].getSubset(region=h_elset_th).values
        for j in range(numel_th):
             tempsts += stresses[j].data[:3]
             tempstn += strains[j].data[:3]
             temppeeq += peeq_set[j].data
             S = stresses[j].data
             tempeqsts += (.5*( (S[0]-S[1])**2 + (S[1]-S[2])**2 + (S[2]-S[0])**2 +
                                6*(S[3]**2 + S[4]**2 + S[5]**2)  ) )**.5

    sts[i] = tempsts/numel_th
    stn[i] = tempstn/numel_th
    peeq[i] = temppeeq/numel_th
    eqsts[i] = tempeqsts/numel_th
    
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
    # Now treat the elements in elset_leprof the sameway

    levals_le = h_All_Frames[i].fieldOutputs['LE'].getSubset(region=h_elset_leprof)
    for j in range(numel_leprof):
        h_el = h_elset_leprof.elements[j]
        if i == 0:
            levals_coord =  h_All_Frames[i].fieldOutputs['COORD'].getSubset(region=h_elset_leprof)
            x,y,z = levals_coord.getSubset(region=h_el).values[0].data
            Rel = np.sqrt(x**2+y**2)
            qel = np.arctan2(y,x)
            le_prof[j,i] = Rel*qel/t
        elLE = levals_le.getSubset(region=h_el).values[0]
        LE = elLE.data
        le_prof[j,i+1] = np.sqrt((2.0/3.0)*(LE*LE).sum())
    if elLE.elementLabel != h_el.label:
        print "Warning.  You're getting an element mismatch in LEprofile"

    # Grab the frame number for the history data
    frameNos[i] = h_All_Frames[i].frameId


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
fname = '%s_NewResults.dat'%(job)
np.savetxt(fname, X = np.c_[F, P, V, sig_x, sig_q, d_lo, d_back, sts, stn, peeq, eqsts], fmt='%.6f', delimiter=', ')
hl = '#[0] Force (kip), [1]Pressure (ksi), [2]Vol, [3]NomAxSts, [4]NomHoopSts,' 
hl += ' [5]d/Lg lo, [6]d/Lg Back, [7,8,9]S11,22,33, [10,11,12]LE11,22,33, [13]PEEQ(SDV1), [14]EqSts(SDV2)'
headerline(fname, hl)

ur_prof = ur_prof[:,ur_prof[0].argsort()]
fname = '%s_UR.dat'%(job)
np.savetxt(fname, X=ur_prof, fmt='%.8f', delimiter=',')
hl = '#1st line: Undef Z-coord of nodes\n'
hl += '#2nd line: Undef r-coord of nodes\n'
hl += '#3+nth line: nth incremend Ur of node in column'
headerline(fname,hl)

le_prof = le_prof[ le_prof[:,0].argsort() ]
fname = '%s_LEprof.dat'%(job)
np.savetxt(fname, X=le_prof, fmt='%.6f', delimiter=', ')
hl = '# First column:  Undeformed Rq/to coord.\n'
hl+= '# (k+1)th column:  LEp'
headerline(fname,hl)

