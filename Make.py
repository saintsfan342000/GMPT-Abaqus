import numpy as n
pi = n.pi
import os
from sys import argv

def println():
    print('-----------------------------------------------------------')

if len(argv) != 7:
    println()
    print('Requires 6 command-line arguments.')
    print('[1] Experiment no. whose geometry and loading you\'re modeling.')
    print('[2] Number of elements thru-thickness in the test section')
    print('[3] Imperfection magnitude (as fraction, not percentage)')
    print('[4] Eccentricity.  Give a number, or give "auto" and the geometry will match that of the test specimen.')
    print('[5] Input file name')
    print('[6] Constitutive model: "vm", "h8", or "anis"')
    println()
else:
    argv, expt, num_el_fine_th, dt, eccen, inpname, constit = argv
    d = n.genfromtxt('ExptParams.dat', delimiter=', ', dtype=str)
    # [0]Expt no, [1]IDg, [2]Mean thickness, [3]Min thickness
    alpha, ID, tg, tmin = d[ d[:,0]==expt, 1:].ravel()   
    if eccen == 'auto':
        eccen = 1 - float(tmin)/float(tg)
        print('Eccen = {:.2f}'.format(eccen*100))
    
    inpname = inpname.split('.')[0] # In case I give an extension.

    Lg = '4.0'
    Ltop = '0.5'  # Length of thick section above radius/chamf
    ODtop = '1.9685'    # Radius of thick section    
    R = '0.5'    # Radius of chamf
    
    ### Nodes
    println()
    print('Generating nodes.')
    os.system('python Nodes.py {:s} {:s} {:s} {:s} {:s} {:s} {:s} {:s}'.format(
                Lg, Ltop, ODtop, ID, tg, R, num_el_fine_th, dt))
    print('Done.')

    ### Elements
    println()
    print('Generating elements.')
    os.system('python Elements.py')
    print('Done.')

    ### Sets
    println()
    print('Generating node and element sets.')
    os.system('python Sets.py {}'.format(Lg))
    print('Done.')

    ### Eccentricity
    println()
    print('Generating Eccentricity.')
    os.system('python Eccen.py {} {} {}'.format(eccen, Lg, R))
    print('Done.')

    ### Write Input
    println()
    print('Writing input file, {}.inp'.format(inpname))
    os.system('python WriteInput.py {:s} {:s} {:s} {:s} {:s} {:s} {:s}'.format(
              expt, inpname, constit, num_el_fine_th, 
              dt, str(eccen), ID ))
    print('Done!')
