import numpy as n
from numpy import pi
import matplotlib.pyplot as p
from sys import argv
import figfun as f
import glob
import os

'''
Compare the strain paths, sts-stns, and ur_profs
'''

expts = [11, 4, 2, 8, 3, 12]
jobno = [2, 1, 2, 6, 2, 1]
key = n.genfromtxt('ExptSummary.dat', delimiter=',')
key = key[ n.in1d(key[:,0],expts) ]
key = key[ n.argsort(key[:,4]) ]

repl = {'88':'11', '99':'8', '00':'8'}

for k, (exp,job) in enumerate(zip(expts, jobno)):

    st = str(exp)
    simpath = 'Jobs/{}/{}'.format(exp,job)
    
    exppath = 'Martin_Experiments/GM/PT/GMPT-{}_FS15SS5'.format(exp)
    if exp in [99, 0, 88]:
        exppath = exppath.replace(st,repl[st])
    
    for i in range(10):
        if not os.path.exists(exppath):
            exppath = '../' + exppath
        else:
            break
    else:
        raise IOError("Can't find the experiment directory!")


    calpath = 'Martin_Experiments/GM/Calibration/GMPT-{}'.format(exp)
    if exp in [99, 0, 88]:
        calpath = calpath.replace(st, repl[st])
    for i in range(10):
        if not os.path.exists(calpath):
            calpath = '../' + calpath
        else:
            break
    else:
        raise IOError("Can't find the calibration directory!")

    relpath = '.'
    for i in range(10):
        if not os.path.exists('{}/ExptSummary.dat'.format(relpath)):
            relpath = '../' + relpath
        else:
            key = n.genfromtxt('{}/ExptSummary.dat'.format(relpath), delimiter=',')
            break

    key = key[ key[:,0] == int(exp) ].ravel()
    Ro, to = key[6:8]


    # Simulation data
    cutoff = n.genfromtxt('{}/Cutoff.dat'.format(simpath),dtype=int)
    d = n.genfromtxt('{}/{}_NewResults.dat'.format(simpath, exp), delimiter=',')
    if cutoff > 0:
        d = d[:cutoff+1]
    #if exp == 8:
    #    d[:,[5,6,12]]*=1.8
    #d[:,3:5]+=1
    #[0] Force (kip), [1]Pressure (ksi), [2]Vol, [3]NomAxSts, [4]NomHoopSts, 
    # [5]d/Lg lo, [6]d/Lg Back, [7,8,9]S11,22,33, [10,11,12]LE11,22,33, 
    # [13]PEEQ(SDV1), [14]EqSts(SDV2)

    d[:,[5,6,10,11,12]]*=100
    simloc = n.argmax(d[:,1])
    u = n.genfromtxt('{}/{}_UR.dat'.format(simpath,exp), delimiter=',').T
    u = u[ u[:,0] <= 2]
    # Mirror so that it looks like the expt plot
    uneg = u.copy()
    uneg[:,0]*=-1
    u = n.vstack((u,uneg))
    u = u[ u[:,0].argsort() ]


    # Experiment data
    xd = n.genfromtxt('{}/Results.dat'.format(exppath), delimiter=',')
    # [0] Stage, [1,2,3]eps_x,q,xq(point avg), [4]eps_x(1"ext), [5]eps_q(BFcirc@mid)
    # [6]d/L, Lg=4.289991
    # d[:,1:3] = n.log(1+d[:,1:3])
    xd[:,1:]*=100
    xstn= n.genfromtxt('{}/WholeFieldAverage.dat'.format(exppath), delimiter=',', usecols=(1,2))
    # Expt nominal hoop stn in area close to ES_ANALZONE
    xstn = n.log(1+xstn)*100
    xex, xeq = xstn.T
    xst = n.genfromtxt('{}/STPF.dat'.format(exppath), delimiter=',', usecols=(range(6)))
    # [0]Stage, [1]Time, [2]Force(kip), [3]Pressure(ksi), [4]NomAxSts(ksi), [5]NomHoopSts(ksi)
    exploc = xst[:,3].argmax()
    xu = n.genfromtxt('{}/ur_profiles.dat'.format(exppath), delimiter=',', usecols=(0,3*exploc+2,-1))
    xu[:,1:]*=100
    # Limitload y-coord, ur/Ro

    alpha, blank = n.polyfit(d[:,4],d[:,3],1)
    annotstr = '{} || {:.2f}'.format(exp,alpha)

    if k == 0:
        p.style.use('mysty-quad')
        p.rcParams['axes.labelsize'] = 18
        p.rcParams['lines.markersize'] = 3
        pad = .9
        hgap = 1.25
        vgap = 2
        axwt = 10/3
        axht = 8/3
        W = 2*pad + hgap + 2*axwt
        H = 2*pad + vgap + 2*axht
        fig = p.figure(figsize=(W,H))
        # Top left
        ax1 = fig.add_axes([pad/W, (pad+vgap+axht)/H, axwt/W, axht/H])
        # Top right
        ax2 = fig.add_axes([(pad+hgap+axwt)/W, (pad+vgap+axht)/H, axwt/W, axht/H])
        # Bottom left
        ax3 = fig.add_axes([pad/W, (pad+vgap/3-1/3)/H, axwt/W, (axht+2/3)/H])
        # Bottom right
        ax4 = fig.add_axes([(pad+hgap+axwt)/W, .8*pad/H, axwt/W, 1.5*axht/H])

    # ax1: Ax sts vs ax stn.  
    # delta/L
    a, = ax1.plot(d[:,5],d[:,3], '--', label=annotstr)
    x,= ax1.plot(xd[:,4], xst[:,4], a.get_color())
    ax1.plot(d[simloc,5],d[simloc,3], 'rD')
    ax1.plot(xd[exploc,4], xst[exploc,4], 'rD')
    #The point field averages
    #ax1.plot(d[:,11],d[:,2], '--', color=simcolor)
    #ax1.plot(n.log(xd[:,1]*.01+1)*100, xst[:,4], '--', color=expcolor)
    # Labels and annotations
    if k == len(expts) - 1:
        ax1.set_xlabel('$\\delta/\\mathsf{L}$ (%)')
        ax1.set_ylabel('$\\Sigma_\\mathsf{x}$\n(ksi)')
        ax1.axis(ymin=0, xmin=-3)
        f.ezlegend(ax1, title='Exp||$\\eta\\prime$', loc='lower right')
        f.eztext(ax1, '$\\mathsf{L}_\\mathsf{g}$=1"', 'ul')
        f.myax(ax1)

    # ax2: Hoop sts vs Hoop stn
    # e_q
    ax2.plot(d[:,11],d[:,4], '--', label=annotstr)
    ax2.plot(xeq, xst[:,5], a.get_color())
    ax2.plot(d[simloc,11], d[simloc,4], 'rD')
    ax2.plot(xeq[exploc], xst[exploc,5], 'rD')
    if k == len(expts) - 1:
        ax2.axis(xmin=0, xmax=10)
        ax2.set_xlabel('e$_\\theta$ (%)')
        ax2.set_ylabel('$\\Sigma_\\theta$\n(ksi)')
        #f.ezlegend(ax2, loc='upper right')
        f.eztext(ax2, 'Avg. of Pts.', 'br')
        f.myax(ax2)

    #[0] Force (kip), [1]Pressure (ksi), [2]Vol, [3]NomAxSts, [4]NomHoopSts, 
    # [5]d/Lg lo, [6]d/Lg Back, [7,8,9]S11,22,33, [10,11,12]LE11,22,33, 
    # [13]PEEQ(SDV1), [14]EqSts(SDV2)
    
    #xex, xeq = n.genfromtxt('{}/CalData.dat'.format(calpath),
    #                        delimiter=',', usecols=(7,8), unpack=1)*100

    #d[:,12] -= (d[:,9]-0.32*d[:,8])/10036
    #d[:,11] -= (d[:,8]-0.32*d[:,9])/10036

    # ax3: Axial vs hoop stn
    ax3.plot(d[:,11], d[:,12], '--', label=annotstr)
    ax3.plot(xeq, xex, a.get_color())
    ax3.plot(d[simloc,11], d[simloc, 12], 'rD')
    ax3.plot(xeq[exploc], xex[exploc], 'rD')
    if k == len(expts) - 1:
        ax3.axis([0,10,-4,6])
        ax3.set_xlabel('e$_\\theta$ (%)')
        ax3.set_ylabel('e$_\\mathsf{x}$\n(%)')
        #f.ezlegend(ax3, loc='upper right')
        f.myax(ax3)

    # Ax4: Ur Profiles
    ax4.plot(100*u[:,(simloc+2)]/u[:,1], u[:,0]*2/4, '--', label=annotstr)
    ax4.plot(xu[:,1], 2*xu[:,0]/4, a.get_color())
    if k == len(expts) - 1:
        ax4.set_xlabel('u$_\\mathsf{r}$/R$_\\mathsf{o}$ (%)')
        ax4.set_ylabel('$\\frac{\\mathsf{2y}_\\mathsf{o}}{\\mathsf{L}_\\mathsf{g}}$')
        ax4.axis(xmin=0,ymin=-1,ymax=1)
        f.eztext(ax4, 'LL', 'br')
        f.myax(ax4, autoscale=.75, nudge=('left',.2,0))

fig.savefig('AllPlot.jpg', dpi=200, bbox_inches='tight')

p.show('all')
