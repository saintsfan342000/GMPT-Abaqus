import numpy as n
from numpy import pi
from pandas import read_csv
import matplotlib.pyplot as p
from sys import argv
import figfun as f
p.style.use('mysty')
import glob
import os

size_factor = .8

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
key = n.genfromtxt('{}/../ExptSummary.dat'.format(exppath), delimiter=',')
key = key[ key[:,0] == int(exp) ].ravel()
Ro, to = key[6:8]

# Simulation data
d = n.genfromtxt('{}_results.dat'.format(name), delimiter=',')
# [0] Force (kip), [1]Pressure (ksi), [2]NomAxSts, [3]NomHoopSts, [4]d/Lg lo, 
# [5]d/Lg hi, [6,7,8]S11,22,33, [9,10,11]LE11,22,33, [12]Vol
d[:,[4,5,9,10,11]]*=100
simloc = n.argmax(d[:,1])
u = n.genfromtxt('{}_UR.dat'.format(name), delimiter=',').T
u = u[ u[:,0] <= 2]
# Mirror so that it looks like the expt plot
uneg = u.copy()
uneg[:,0]*=-1
u = n.vstack((u,uneg))
u = u[ u[:,0].argsort() ]
#u[:,2:]*=(1/u[:,1])
lep = n.genfromtxt('{}_LEprof.dat'.format(name), delimiter=',')
lepneg = lep.copy()
lepneg[:,0]*=-1
lep = n.vstack((lep,lepneg))
lep = lep[ lep[:,0].argsort() ]


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
xu = n.genfromtxt('{}/ur_profiles.dat'.format(exppath), delimiter=',', usecols=(0,3*exploc+2, -1))
xu[:,1:]*=100
# Limitload y-coord, ur/Ro
xlep = n.genfromtxt('{}/LEp_profiles.dat'.format(exppath), delimiter=',')


alpha, blank = n.polyfit(d[:,3],d[:,2],1)
annotstr = '{}.  $\\alpha$ = {:.2f}'.format(name,alpha)
   
p.style.use('mysty-12')
p.rcParams['font.size'] = 18*size_factor
p.rcParams['axes.labelsize'] = 22*size_factor

W,H = 18,10
fvec = n.array([W,H,W,H])

# ax1:  LEp Profile
x,y,w,h = 1.5,1.2,10,4
ax1_loc = n.array([x,y,w,h]/fvec)
# ax2:  Ur Profile
x,y,w,h = 12.75,2,3.5,6
ax2_loc = n.array([x,y,w,h]/fvec)
# ax3:  Hooop-sts / hoop-stn
x,y,w,h = 7.25,6.1,4,3
ax3_loc = n.array([x,y,w,h]/fvec)
# ax4:  Ax-sts / hoop-sts
x,y,w,h = 1.75,6.1,4,3
ax4_loc = n.array([x,y,w,h]/fvec)
# ax_sl:  Slider
x,y,w,h = 4,.1,10,.25
ax_sl_loc = n.array([x,y,w,h]/fvec)

W,H = map(lambda x: x*size_factor,(W,H))
fig = p.figure(figsize=(W,H))
for i in [1,2,3,4,'_sl']:
    exec('ax{} = fig.add_axes(ax{}_loc)'.format(i,i))

expcolor = 'C0'
simcolor = 'C1'

# ax1:  LEp Profile
# Exp at LL and failure
ax1.plot(xlep[:,0]/to, xlep[:,exploc+1], expcolor)
ax1.plot(xlep[:,0]/to, xlep[:,-1], expcolor)
# Initial Sim
line1, = ax1.plot(lep[:,0],lep[:,1], simcolor)
ax1.axis([-8,8,0,n.nanmax(xlep[:,1:])*1.05])
ax1.set_xlabel('s/t$_\\mathsf{o}$')
ax1.set_ylabel('e$_\\mathsf{e}$')
f.myax(ax1, autoscale='preserve')

# ax2:  Ur Profile
# Exp at LL and failure
ax2.plot(xu[:,1], 2*xu[:,0]/4, expcolor, label='Exp')
ax2.plot(xu[:,2], 2*xu[:,0]/4, expcolor)
# Initial sim
line2, = ax2.plot(100*u[:,(0+2)]/u[:,1], 2*u[:,0]/4, simcolor, label='Anal')
ax2.axis(xmin=0, ymin=-1, ymax=1)
ax2.set_xlabel('u$_\\mathsf{r}$/R$_\\mathsf{o}$ (%)')
ax2.set_ylabel('$\\frac{\\mathsf{2y}_\\mathsf{o}}{\\mathsf{L}_\\mathsf{g}}$')
f.ezlegend(ax2, loc='lower center', bbox_to_anchor=(0.5, 1.05), fontsize=20)
f.myax(ax2, nudge=('left',.3,.5), autoscale='preserve')

# ax3:  Hooop-sts / hoop-stn
# Exp. response and LL
ax3.plot(xeq, xst[:,5], expcolor)
ax3.plot(xeq[exploc], xst[exploc,5], '^', mfc=expcolor, mec='k')
# Sim response and LL
ax3.plot(d[:,10], d[:,3], simcolor)
ax3.plot(d[simloc,10], d[simloc,3], '^', mfc=simcolor, mec='k')
# Sim marker initial value
line3, = ax3.plot(d[0,10], d[0,3], 'o', mec=simcolor)
#if not exp in [1,4]:
#    ax3.axis(xmin=0,ymin=0)
ax3.axis(xmin=0)
ax3.set_xlabel('$\\bar{\\epsilon}_\\theta$ (%)')
ax3.set_ylabel('$\\sigma_\\theta$\n($\\mathsf{ksi}$)')
f.myax(ax3, autoscale='preserve')

# ax4:  Ax-sts / ax-stn
# Exp reponse and LL
ax4.plot(xd[:,4], xst[:,4], expcolor)
ax4.plot(xd[exploc,4], xst[exploc,4], '^', mfc=expcolor, mec='k')
# Sim response and LL
ax4.plot(d[:,4],d[:,2], simcolor)
ax4.plot(d[simloc,4],d[simloc,2], '^', mfc=simcolor, mec='k')
# Sim marker initial value
line4, = ax4.plot(d[0,4], d[0,2], 'o', mec=simcolor)
ax4.axis(ymax=1.1*xst[:,4].max())
#if not exp in [1,4]:
#    ax4.axis(xmin=0,ymin=0, ymax=1.05*xst[:,4].max())
ax4.axis(ymin=0)
ax4.set_xlabel('$\\bar{\\epsilon}_\\mathsf{x}$ (%)')
ax4.set_ylabel('$\\sigma_{\\mathsf{x}}$\n($\\mathsf{ksi}$)')
f.myax(ax4, autoscale='preserve')

fig.text(.5, .98, annotstr, ha='center', va='top', transform=fig.transFigure)
#fig.text(.5, .95, 'P = {:.0f} psi'.format(P[0]*1000), ha='center', va='top', transform=fig.transFigure)

from matplotlib.widgets import Slider, Button, RadioButtons
slider = Slider(ax=ax_sl, label='', valmin=0, valmax=d.shape[0]-1, valinit=0, valfmt='%.0f', facecolor=line1.get_color())

def update(val):
    i = int(val)
    line1.set_ydata(lep[:,1+i])
    line2.set_xdata(100*u[:,i+2]/u[:,1])
    line3.set_data(d[i,10], d[i,3])
    line4.set_data(d[i,4], d[i,2])
    fig.canvas.draw_idle()

slider.on_changed(update)
p.show()

