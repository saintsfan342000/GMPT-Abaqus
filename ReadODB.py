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

try:
    job, R, t = argv[1:]
    job = job.split('.')[0] # In case I give .odb with the name
    R = float(R)
    t = float(t)
except IndexError:
    job = 'job'
    t = 0.05
    R = 1.6285/2 + t/2

if not os.path.isfile(job + '.odb'):
    raise ValueError('The specified job name "%s" does not exist.'%(job))

nset_rp_top = 'NS_RPTOP'
nset_rp_bot = 'NS_RPBOT'
nset_dr_lo = 'NS_DISPROT_LO'
nset_dr_hi = 'NS_DISPROT_HI'
nset_dr_new = 'NS_DISPROT_NEW'
elset_th = 'ES_THICKNESS'

h_odb = O.openOdb(job + '.odb',readOnly=True)
h_inst = h_odb.rootAssembly.instances[ h_odb.rootAssembly.instances.keys()[0] ]
h_step = h_odb.steps[ h_odb.steps.keys()[0] ]
h_All_Frames = h_step.frames
num_incs = len(h_All_Frames)

h_nset_rp_top = h_odb.rootAssembly.nodeSets[nset_rp_top]
h_nset_dr_lo = h_inst.nodeSets[nset_dr_lo]
h_nset_dr_hi = h_inst.nodeSets[nset_dr_hi]
h_nset_dr_new = h_inst.nodeSets[nset_dr_new]
h_elset_th = h_inst.elementSets[elset_th]

F = np.empty( (num_incs) )
P = np.empty( (num_incs) )
d_lo = np.empty( (num_incs) )
d_hi = np.empty( (num_incs) )

# Grab undef coords of dr_lo and hi
Lg_lo = h_nset_dr_lo.nodes[0].coordinates[2]
Lg_hi = h_nset_dr_hi.nodes[0].coordinates[2]


for i in range(num_incs):
    F[i] = h_All_Frames[i].fieldOutputs['CF'].getSubset(region=h_nset_rp_top).values[0].data[2]
    P[i] = h_All_Frames[i].fieldOutputs['P'].values[0].data
    d_lo[i] = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_dr_lo).values[0].data[2]
    d_hi[i] = h_All_Frames[i].fieldOutputs['U'].getSubset(region=h_nset_dr_hi).values[0].data[2]

h_odb.close()

# Convert disp to d/Lg
d_lo, d_hi = d_lo/Lg_lo, d_hi/Lg_hi

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
np.savetxt('FPD.dat',X = np.vstack((F, P, F/(2*pi*R*t), P*R/t, d_lo, d_hi)).T, fmt='%.15f', delimiter=', ')
headerline('FPD.dat','#[0] Force (kip), [1]Pressure (ksi), [2]Nom AxSts, [3]Nom Shear Sts, [4]d/Lg lo, [5]d/Lg hi\n')
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
