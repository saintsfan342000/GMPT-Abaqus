import numpy as n
from numpy import (linspace, asarray, in1d, empty,
                    hstack , vstack, pi, compress)
from scipy.spatial import Delaunay as trimesh
from sys import argv
from pandas import read_excel
import time

'''
Writes the input file
'''

fullring = False

if len(argv[1:]) == 12:
    expt = int(argv[1])
    inpname = argv[2]
    constit = argv[3]
    num_el_fine_th = int(argv[4])
    dt = float(argv[5])
    ecc = float(argv[6])
    ID = float(argv[7])
    R = float(argv[8])
    t = float(argv[9])
    alpha = float(argv[10])
    cdt = float(argv[11])
    impN = (argv[12])
else:
    raise ValueError('Wrong number of args. Require 11')

if n.isclose(alpha, 0.5):
    UAMP = False
else:
    UAMP = True

# Make sure we have a valid constitutive model
if not( constit in ['vm', 'VM', 'H8', 'h8', 'anis', 'ANIS']):
    raise ValueError("Bad constit given '{}'.\nMust be 'vm', 'VM', 'H8', 'anis', 'ANIS'.".format(constit))
    
press = 1400  # This is 2.8 ksi * 500
# The *500 is b/c abaqus behaves odd when the cloads are O(1) 
# The behavior is normal when the cloads are O(1000)
# I do 500 rather than 1000 just to have better resolution on the LPF
# The max LPF is very small accordingly 
K = (2*alpha-1)*(pi*R*R)  
force = K*press # force
# Disp. control:  Monitor axial disp
riks_DOF_num = 3
riks_DOF_val = .4

print('Cload 3 magnitude: {:.2f}'.format(force))
print('Dsload magnitude:  {:.2f}'.format(press))
print('Max displacement of Riks Node = {:.8f}'.format(riks_DOF_val))

# Load up the node and element lists
nodelist = n.load('./ConstructionFiles/abaqus_nodes.npy')
elemlist = n.load('./ConstructionFiles/abaqus_elements.npy')

fid =  open('{}.inp'.format(inpname),'w')

## Info for me
fid.write('** GMPT={}\n'.format(expt))
fid.write('** alpha = {:.3f}\n'.format(alpha))
fid.write('** num_el_fine_th = {}\n'.format(num_el_fine_th))
fid.write('** Axial Imperfection = {} #Not percentage\n'.format(dt))
fid.write('** Circumferential Imperfection = {} #Not percentage\n'.format(cdt))
fid.write('** Eccentricity = {:.2f} # Percentage\n'.format(100*ecc))
fid.write('** Rm = {:.4f}\n'.format(R))
fid.write('** tg = {:.4f}\n'.format(t))
fid.write('** ID = {}\n'.format(ID))
fid.write('** constit = "{}"\n'.format(constit))
fid.write('** Fluid Cavity and elements.  If alpha = 0.5, then flux control\n')
fid.write('** {} nodes, {} elements.\n'.format(nodelist.shape[0], elemlist.shape[0]))
z = time.localtime()
fid.write('** Circ. Imperf. Width = {}*tg\n'.format(impN))
fid.write('** Generated on {}/{}/{} at {}:{}\n'.format(z[1],z[2],z[0],z[3],z[4]))

## HEADING and PREPRINT
fid.write('*Heading\n' +
          'Z is Tube Axis\n'  +
          'Min Wall Thickness is (x>0,y=0)\n'
          )
fid.write('*Preprint, echo=NO, model=NO, history=NO, contact=NO\n')

###################
## PART and sets
###################
fid.write('****************************************\n')
fid.write('***************** PART *****************\n')
fid.write('****************************************\n')
fid.write('*part, name=PART\n')
# Nodes
fid.write('*node, nset=NS_ALLNODES\n')
for i in nodelist:
    fid.write('{:.0f}, {:.12f}, {:.12f}, {:.12f}\n'.format(i[0],i[1],i[2],i[3]))
# Elements
fid.write('*element, type=C3D8R, elset=ES_ALLELEMENTS\n')
for i in elemlist:
    fid.write('{}, {}, {}, {}, {}, {}, {}, {}, {}\n'.format(
               i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8])
              )
# Node and Element Sets
with open('./ConstructionFiles/abaqus_sets.txt','r') as setfid:
    sets = setfid.read()
    fid.write(sets)
    setfid.close()

# Reference points
nodenum, zcoord = nodelist[:,0].max()+1, nodelist[:,3].max()
fid.write('*node, nset=NS_RPTOP\n' +
          '{:.0f}, 0, 0, {:.12f}\n'.format(nodenum, zcoord)
          )
fid.write('*node, nset=RP_CAV\n' + 
          '{:.0f}, 0, 0, 0\n'.format(nodenum+1)
          )
# Now generate the surface elements that connect from top ID nodes to top RP
# From Abq. man 11.5.1: "The boundary of the fluid cavity is defined by an 
# element-based surface with normals pointing to the inside of the cavity"
# From Abq. man 32.7.1:  "For general surface elements the positive normal 
# direction is defined by the right-hand rule going around the nodes of the 
# element in the order that they are specified in the element definition
# SO looking at tube from top down, I need to process from (x,y)=(-ID,0) to (ID,0)

cnodes = n.load('./ConstructionFiles/ni_cors.npy')
# Inner diam and max z coord
rng = (cnodes[:,0] == 0) & (cnodes[:,1] == cnodes[:,1].max())
nodes = cnodes[rng,3] # node num
nodes = nodelist[nodes-1] # From master node list
# Reorder nodes.  Nodenum increases from q=0 to pi.  I need to go in reverse on connecting.
nodes = n.flipud(nodes)
fid.write('*element, type=SFM3D3, elset=ES_SURFELS\n')
e1 = elemlist[:,0].max()
for k,i in enumerate(nodes[:-1,0]):
    fid.write('{:d}, {:.0f}, {:.0f}, {:.0f}\n'.format(e1+k+1,
                nodenum, nodes[k,0], nodes[k+1,0])
         )
if fullring:
    # Connect back to zero.  Not tested.
    fid.write('{:d}, {:.0f}, {:.0f}, {:.0f}'.format(el+k+2),
                nodenum, nodes[k+1,0], nodes[0,0])

# Orientation, transformation, section
fid.write('*orientation, name=ANISOTROPY, system=cylindrical, definition=coordinates\n' +
          '0, 0, 0, 0, 0, 1\n'
          )
fid.write('*solid section, elset=ES_ALLELEMENTS, material=MATERIAL, orientation=ANISOTROPY\n')
fid.write('*Hourglass Stiffness\n' +
          '40.0, , , \n'
          )
fid.write('*transform, nset=NS_ALLNODES, type=C\n' +
          '0, 0, 0, 0, 0, 1\n'
          )
fid.write('*surface section, elset=ES_SURFELS\n')
fid.write('*end part\n')
# end part

###################
#### ASSEMBLY ####
###################
fid.write('****************************************\n')
fid.write('*************** Assembly ***************\n')
fid.write('****************************************\n')
fid.write('*assembly, name=ASSEMBLY\n')
fid.write('*instance, name=INSTANCE, part=PART\n')
fid.write('*end instance\n')

# transform for RPs
fid.write('*transform, nset=INSTANCE.NS_RPTOP, type=C\n' +
          '0, 0, 0, 0, 0, 1\n')
nodenum+=1
# Riks monitoring point, must be defined in the assembly
fid.write('** Riks displacement monitoring node must be defined in the assembly\n' + 
          '*nset, nset=RIKSMON, instance=INSTANCE\n' + 
          'NS_DISPROT_LO\n')
# Surfaces
fid.write('*surface, type=node, name=SURF_TOPSURFACE\n' +
          'INSTANCE.NS_TOPSURFACE\n'
          )
fid.write('*surface, type=node, name=SURF_BOTSURFACE\n' +
          'INSTANCE.NS_BOTTOMSURFACE\n'
          )
fid.write('*surface, type=element, name=SURF_IDSURFACE\n' +
          'INSTANCE.ES_WHOLEID, S6\n')
# Kinematic coupling
fid.write('*orientation, name=ORI_COUPLING, system=cylindrical, definition=coordinates\n' +
          '0, 0, 0, 0, 0, 1\n'
          )
fid.write('*coupling, constraint name=CONSTRAINT_TOPSURFACE, ref nod=INSTANCE.NS_RPTOP, surface=SURF_TOPSURFACE, orientation=ORI_COUPLING\n' +
          '*kinematic\n'
          )
fid.write('** Surface definition for the fluid cavity\n' + 
          '*surface, type=element, name=SURF_CAVITY\n' + 
          'INSTANCE.ES_SURFELS, SPOS\n' + 
          'INSTANCE.ES_WHOLEID, S6\n'
          )
fid.write('*end assembly\n')
# end assembly

###################
#### MATERIAL #####
###################
fid.write('****************************************\n')
fid.write('***************  MATERIAL **************\n')
fid.write('****************************************\n')

if constit in ['vm', 'VM']:
    with open('./ConstructionFiles/abaqus_material_VM.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()
elif constit in ['h8', 'H8']:
    with open('./ConstructionFiles/abaqus_material_H8.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()    
elif constit in ['ANIS', 'anis']:
    with open('./ConstructionFiles/abaqus_material_ANISp1.txt','r') as matfid:
        mat = matfid.read()
        fid.write(mat)
        matfid.close()

########################
### Fluid and cavity ###
########################
fid.write('****************************************\n')
fid.write('*************** FLUID CAV  *************\n')
fid.write('****************************************\n')
fid.write('*fluid behavior, name=FLUID\n')
fid.write('*fluid density\n' + 
          '1.0\n'
         )
# You need a blank line after fluid cav...two newline characters!
fid.write('*fluid cavity, name=FLUIDCAVITY, behavior=FLUID, refnode=INSTANCE.RP_CAV, surface=SURF_CAVITY\n\n')

###################
### INITIAL BCs ###
###################
fid.write('*boundary\n')
fid.write('** Top Ref Pt:  Only U3, UR3 permitted\n')
fid.write('INSTANCE.NS_RPTOP, 1, 2\n' +
          'INSTANCE.NS_RPTOP, 4, 5\n'
          )
# Specify boundary condition if halfring
if not fullring:
    fid.write('** AxialSymmetry Nodes:  Same as YSYMM\n')
    fid.write('INSTANCE.NS_AXIALSYMM, 2, \n' + 
              '**INSTANCE.NS_AXIALSYMM, 4, \n' +
              '**INSTANCE.NS_AXIALSYMM, 6, \n'
             )
# Bottom symplane conditions...order matters
fid.write('** Bottom surface nodes:  Same as ZSYMM\n')
fid.write('INSTANCE.NS_BOTTOMSURFACE, 3, 3 \n')
# Cavity ref node
fid.write('** Cavity reference node fully constrained\n')
fid.write('INSTANCE.RP_CAV, 1, 6\n')

###################
###### STEP #######
###################
fid.write('****************************************\n')
fid.write('*************** STEP *******************\n')
fid.write('****************************************\n')
if UAMP:
    fid.write('*amplitude, name=FORCE_AMP, Definition=User, Variables=1\n')

if constit in ['H8','h8','anis','ANIS']:
    numinc, freq = 10000, 5
else:
    numinc, freq = 1000, 1
fid.write('*step, name=STEP, nlgeom=yes, inc={}\n'.format(numinc))
if n.isnan(alpha):
    # Pure tension case, apply only U3
    # NOT TESTED
    fid.write('*static\n' +
              '0.005, 1., 1e-05, .005\n'
              )
    fid.write('**[1]Initial incr, [2]Total step, [3]Min incr, [4]Max incr\n')
    fid.write('*boundary\n' + 
              'INSTANCE.RP_TOP, 3, 3, 0.3\n'
              )
elif n.isclose(alpha, 0.5):
    # Just pressure.  Apply fluid flux only
    # Calculate volume of cavity
    vol = pi*R*R*zcoord
    if not fullring:
        vol*=.5
    fid.write('*static\n' +
              '0.00125, 1., 1e-06, .002\n'
              )
    fid.write('**[1]Initial incr, [2]Total step, [3]Min incr, [4]Max incr\n')
    fid.write('*fluid flux\n' + 
              'INSTANCE.RP_CAV, {:.3f}\n'.format(vol/5/2)
              )
else:
    if not UAMP:
    # Pressurizing and pulling.  Must use riks or UAMP
        # [1]Inital arc len, [2]total step, [3]minimum increm, [4]max increm (no max if blank), [5]Max LPF, [6]Node whose disp is monitored, [7]DOF, [8]Max Disp
        fid.write('*static, riks\n' +
                '0.0005, 1.0, 1e-15, , .0022, ASSEMBLY.RIKSMON, '+
                '{:.0f}, {:.6f}\n'.format(riks_DOF_num, riks_DOF_val)
                 )
        fid.write('**[1]Inital arc len, [2]total step, [3]minimum increm, [4]max increm'+
                    '(no max if blank),' +
                  ' [5]Max LPF, [6]Node whose disp is monitored, [7]DOF, [8]Max Disp\n')
        if not fullring:
            force*=.5
        fid.write('*cload\n' +
                      'INSTANCE.NS_RPTOP, 3, {:.5f}\n'.format(force) 
                  )
        fid.write('*boundary\n' + 
                  'INSTANCE.RP_CAV, 8, 8, {:.5f}\n'.format(press)
                  )
    elif UAMP:
        vol = pi*R*R*zcoord
        if not fullring:
            vol*=.5
        fid.write('*static\n'+
                  '0.000325, 1.0, 1e-06, 0.000325\n'
                  )
        fid.write('*fluid flux\n' + 
                'INSTANCE.RP_CAV, {:.3f}\n'.format(2*vol/5/2)
                  )
        if not fullring:
            force *= 0.5
        fid.write('*cload, op=MOD, amplitude=FORCE_AMP\n'
                  'INSTANCE.NS_RPTOP, 3, {:.6f}\n'.format(force/press)
                  )

# Solver controls
fid.write('** Increase number of increment cutbacks permitted\n')
fid.write('*controls, parameters=time incrementation\n')
fid.write(' , , , , , , , 10, , , , ,\n')           

# field output
fid.write('*output, field, frequency={}\n'.format(freq))
fid.write('** COORn must be called under history output, but COORD can be called in field\n')
fid.write('*node output, nset=INSTANCE.NS_DISPROT_LO\n' +   # disprot nodesets
          'U, UR\n'
          )
fid.write('*node output, nset=INSTANCE.NS_DISPROT_HI\n' +   # disprot nodesets
          'U, UR\n'
          )
fid.write('*node output, nset=INSTANCE.NS_ALLNODES\n' +
          'U, UR\n'
          )
fid.write('*node output, nset=INSTANCE.NS_RPTOP\n' +    # refpt node
          'U, UR, CF\n'
          )
fid.write('*element output, elset=INSTANCE.ES_LEPROF\n' + 
            'COORD\n'
          )
fid.write('*element output, elset=INSTANCE.ES_ALLELEMENTS, directions=YES\n' +    # sts, stn in element sets
           'S, PE, LE, P'
          )
if constit in ['H8','h8','anis','ANIS']:
    fid.write(', SDV1, SDV2\n')
else:
    fid.write('\n')

# History output.  Different for RIKS vs UAMP
if not UAMP:
    fid.write('*output, history, frequency=1\n')
    fid.write('*node output, nset=INSTANCE.RP_CAV\n' + 
              'CVOL, PCAV\n'
             )
elif UAMP:
    # Then we need separate histpory requests for PCAV and CVOL
    fid.write('*output, history, frequency=1, sensor, name=PRESS\n')
    fid.write('*node output, nset=INSTANCE.RP_CAV,global=NO\n' +
              'PCAV\n'
              )
    fid.write('*output, history, frequency=1\n')
    fid.write('*node output, nset=INSTANCE.RP_CAV\n' + 
              'CVOL\n'
             )
fid.write('** Note:  This node print will result in abaqus printing\n' + 
          '** Nearly 80000 lines worth of meaningless nodal-connectivity\n' + 
          '** information to the dat file.\n'
          )
fid.write('*node print, nset=INSTANCE.RP_CAV, summary=NO\n' +
          'CVOL\n'
          )

fid.write('*end step\n')
# end step
fid.close()

