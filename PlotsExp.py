import numpy as n
from numpy import pi
import matplotlib.pyplot as p
from sys import argv
import figfun as f
p.style.use('mysty-quad')
p.rcParams['axes.labelsize'] = 18
import glob
import os

'''
Compare the strain paths, sts-stns, and ur_profs
'''

if len(argv) != 2:
    print('Give name of job, which is prefix to the *_results.dat' + 
            'and *_UR.dat files.  Proper opertation is for your pwd'+
            ' to be the folder where the *dat files are.  If you call'+
            '"$ python ../../Plots.py <jobname>" it will work.')
    raise IndexError('')
    
name = argv[1]
exp = name.split('-')[-1]
exppath = 'Martin_Experiments/GM/PT/GMPT-{}_FS15SS5'.format(exp)

# Find the experimental relative path to experiments
for i in range(10):
    if not os.path.exists('{}'.format(exppath)):
        exppath = '../' + exppath
    else:
        break
else:
    raise IOError("Can't find the experiment directory!")

# Simulation data
d = n.genfromtxt('{}_results.dat'.format(name), delimiter=',')
# [0] Force (kip), [1]Pressure (ksi), [2]NomAxSts, [3]NomHoopSts, [4]d/Lg lo, 
# [5]d/Lg hi, [6,7,8]S11,22,33, [9,10,11]LE11,22,33, [12]Vol
d[:,[4,5,9,10,11]]*=100
u = n.genfromtxt('{}_UR.dat'.format(name), delimiter=',').T
# Mirror so that it looks like the expt plot
uneg = u.copy()
uneg[:,0]*=-1
u = n.vstack((u,uneg))
u = u[ u[:,0].argsort() ]
simloc = n.argmax(d[:,1])

# Experiment data
xd = n.genfromtxt('{}/Results.dat'.format(exppath), delimiter=',')
# [0] Stage, [1,2,3]eps_x,q,xq(point avg), [4]eps_x(1"ext), [5]eps_q(BFcirc@mid)
# [6]d/L, Lg=4.289991
xd[:,1:]*=100
xstn= n.genfromtxt('{}/WholeFieldAverage.dat'.format(exppath), delimiter=',', usecols=(1,2))
# Expt nominal hoop stn in area close to ES_ANALZONE
xstn = n.log(1+xstn)*100
xex, xeq = xstn.T
xst = n.genfromtxt('{}/STPF.dat'.format(exppath), delimiter=',', usecols=(range(6)))
# [0]Stage, [1]Time, [2]Force(kip), [3]Pressure(ksi), [4]NomAxSts(ksi), [5]NomHoopSts(ksi)
exploc = xst[:,3].argmax()
xu = n.genfromtxt('{}/ur_profiles.dat'.format(exppath), delimiter=',', usecols=(0,3*exploc+2))
xu[:,1]*=100
# Limitload y-coord, ur/Ro

alpha, blank = n.polyfit(d[:,3],d[:,2],1)
annotstr = '{}.  $\\alpha$ = {:.2f}'.format(name,alpha)


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
ax3 = fig.add_axes([pad/W, (pad+vgap/3)/H, axwt/W, axht/H])
# Bottom right
ax4 = fig.add_axes([(pad+hgap+axwt)/W, .8*pad/H, axwt/W, 1.5*axht/H])

fig.suptitle(annotstr)

# Fig1: Ax sts vs ax stn.  
# delta/L
a, = ax1.plot(d[:,4],d[:,2], label='Anal')
simcolor = a.get_color()
x,= ax1.plot(xd[:,4], xst[:,4], label='Exp')
expcolor = x.get_color()
ax1.plot(d[simloc,4],d[simloc,2], '^', mfc=simcolor, mec='k')
ax1.plot(xd[exploc,4], xst[exploc,4], '^', mfc=expcolor, mec='k')
#The point field averages
#ax1.plot(d[:,11],d[:,2], '--', color=simcolor)
#ax1.plot(n.log(xd[:,1]*.01+1)*100, xst[:,4], '--', color=expcolor)
# Labels and annotations
ax1.set_xlabel('$\\delta/\\mathsf{L}$')
ax1.set_ylabel('$\\Sigma_\\mathsf{x}$\n(ksi)')
f.ezlegend(ax1, loc='upper right')
f.eztext(ax1, '$\\mathsf{L}_\\mathsf{g}$=1"', 'ul')
f.myax(ax1)

# Fig2: Hoop sts vs Hoop stn
# delta/L
ax2.plot(d[:,10],d[:,3], label='Anal')
ax2.plot(xeq, xst[:,5], label='Exp')
ax2.plot(d[simloc,10], d[simloc,3], '^', mfc=simcolor, mec='k')
ax2.plot(xeq[exploc], xst[exploc,5], '^', mfc=expcolor, mec='k')
ax2.axis(xmin=0)
ax2.set_xlabel('e$_\\theta$')
ax2.set_ylabel('$\\Sigma_\\theta$\n(ksi)')
f.ezlegend(ax2, loc='upper right')
f.eztext(ax2, 'Avg. of Pts.', 'br')
f.myax(ax2)

# Fig3: Axial vs hoop stn
ax3.plot(d[:,10], d[:,11], label='Anal')
ax3.plot(xeq, xex, label='Exp')
ax3.plot(d[simloc,10], d[simloc, 11], '^', mfc=simcolor, mec='k')
ax3.plot(xeq[exploc], xex[exploc], '^', mfc=expcolor, mec='k')
ax3.axis(xmin=0)
ax3.set_xlabel('e$_\\theta$')
ax3.set_ylabel('e$_\\mathsf{x}$')
f.ezlegend(ax3, loc='upper right')
f.myax(ax3)


ax4.plot(100*u[:,(simloc+2)]/u[:,1], u[:,0]*2/4, label='Anal')
ax4.plot(xu[:,1], 2*xu[:,0]/4, label='Exp')
ax4.set_xlabel('u$_\\mathsf{r}$/R$_\\mathsf{o}$ (%)')
ax4.set_ylabel('$\\frac{\\mathsf{2y}_\\mathsf{o}}{\\mathsf{L}_\\mathsf{g}}$')
ax4.axis(xmin=0,ymin=-1,ymax=1)
f.eztext(ax4, annotstr, 'ur')
f.myax(ax4, autoscale=.75, nudge=('left',.2,0))

p.savefig('F5_ExpPlot.png', dpi=125, bbox_inches='tight')
