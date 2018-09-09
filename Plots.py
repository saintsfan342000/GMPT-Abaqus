import numpy as n
from numpy import pi
import matplotlib.pyplot as p
from sys import argv
import figfun as f
import glob
import os

if len(argv) != 2:
    print('Give name of job, which is prefix to the *_results.dat' + 
            'and *_UR.dat files.  Proper opertation is for your pwd'+
            ' to be the folder where the *dat files are.  If you call'+
            '"$ python ../../Plots.py <jobname>" it will work.')
    raise IndexError('')
    
name = argv[1]

d = n.genfromtxt('{}_NewResults.dat'.format(name), delimiter=',').T
#[0] Force (kip), [1]Pressure (ksi), [2]Vol, [3]NomAxSts, [4]NomHoopSts, 
# [5]d/Lg lo, [6]d/Lg Back, [7,8,9]S11,22,33, [10,11,12]LE11,22,33, 
# [13]PEEQ(SDV1), [14]EqSts(SDV2)
u = n.genfromtxt('{}_UR.dat'.format(name), delimiter=',')
# Mirror so that it looks like the expt plot
uneg = u.copy()
uneg[0]*=-1
u = n.hstack((u,uneg))
u = u[:, u[0].argsort() ]

loc = n.argmax(d[1])
#d = d[:,:loc+20]

alpha, blank = n.polyfit(d[4],d[3],1)
annotstr = '{}\n$\\alpha$ = {:.2f}'.format(name,alpha)

incs = n.linspace(0,loc,5)[1:].astype(int)
incs = n.hstack((incs,len(d[0])-1))

p.style.use('mysty-quad')
p.rcParams['axes.labelsize'] = 18
pad = .9
hgap = 1.25
vgap = 2
axwt = 10/3
axht = 8/3
W = 2*pad + hgap + 2*axwt
H = 2*pad + vgap + 2*axht
fig1 = p.figure(figsize=(W,H))
# Top left
ax1 = fig1.add_axes([pad/W, (pad+vgap+axht)/H, axwt/W, axht/H])
# Top right
ax3 = fig1.add_axes([(pad+hgap+axwt)/W, (pad+vgap+axht)/H, axwt/W, axht/H])
# Bottom left
ax2= fig1.add_axes([pad/W, (pad+vgap/3)/H, axwt/W, axht/H])
# Bottom right
ax4 = fig1.add_axes([(pad+hgap+axwt)/W, .8*pad/H, axwt/W, 1.5*axht/H])

# Ur profiles
p.style.use('mysty-sub')
p.rcParams['font.size'] = 18
p.rcParams['axes.labelsize'] = 22
colors = []
for j in incs:
    l, = ax4.plot(100*u[j+2]/.8393,2*u[0]/4)
    colors.append(l.get_color())
ax4.set_xlabel('u$_\\mathsf{r}$/R$_\\mathsf{o}$ (%)')
ax4.set_ylabel('$\\frac{\\mathsf{2y}_\\mathsf{o}}{\\mathsf{L}_\\mathsf{g}}$')
ax4.axis(xmin=0,ymin=-1,ymax=1)
f.eztext(ax4, annotstr, 'ur')
f.myax(ax4, autoscale=.65, nudge=('left',.3,-.4))

#[0] Force (kip), [1]Pressure (ksi), [2]Vol, [3]NomAxSts, [4]NomHoopSts, 
# [5]d/Lg lo, [6]d/Lg Back, [7,8,9]S11,22,33, [10,11,12]LE11,22,33, 
# [13]PEEQ(SDV1), [14]EqSts(SDV2)

#[0] Force (kip), [1]Pressure (ksi), [2]NomAxSts, [3]NomHoopSts, [4]d/Lg lo, 
# [5]d/Lg Back, [6,7,8]S11,22,33, [9,10,11]LE11,22,33, [12]Vol

# Ax vs hoop sts
ax1.plot(d[4], d[3], 'b',label='Nominal')
ax1.plot(d[8], d[9], 'r', label='True')
for q,j in enumerate(incs):
    ax1.plot(d[4, j],d[3,j],marker='o',mfc=colors[q],mec='b')
    ax1.plot(d[8, j],d[9,j],marker='o',mfc=colors[q],mec='r')
ax1.set_xlabel('$\\sigma_\\theta$ (ksi)')
ax1.set_ylabel('$\\sigma_\\mathsf{x}$\n(ksi)')
f.ezlegend(ax1, loc='lower right')
newstr= '{}\nP$_\\mathsf{{max}}$ = {:.0f} psi'.format(annotstr,d[1,-1]*1000)
f.eztext(ax1, newstr, 'tl')
f.myax(ax1)

# Axial vs hoop stn
ax3.plot(d[11], d[12], 'b')
for q,j in enumerate(incs):
    ax3.plot(d[11, j],d[12,j],marker='o',mfc=colors[q],mec=colors[q])  
ax3.set_xlabel('$\\mathsf{e}_\\theta$')
ax3.set_ylabel('$\\mathsf{e}_\\mathsf{x}$')
ax3.axis('equal')
f.eztext(ax3, annotstr, 'br')
f.myax(ax3)

# P vs dv/vo
v0 = d[2,0]
dv = (d[2] - v0)/v0
ax2.plot(dv, d[1]*1000)
for q,j in enumerate(incs):
    ax2.plot(dv[j], d[1,j]*1000, marker='o',mfc=colors[q],mec=colors[q])
ax2.set_xlabel('$\\Delta\\mathsf{V}/\\mathsf{V}_\\mathsf{o}$')
ax2.set_ylabel('P\n(psi)')
f.myax(ax2)


# Ax sts vs ax stn
fig2, ax21, ax22 = f.make21()
ax21.plot(d[5],d[3],'b', label='$\\delta/\\mathsf{L}$')
ax21.plot(d[12],d[3],'r', label='$e_\\mathsf{x}$')
for q,j in enumerate(incs):
        ax21.plot(d[5, j],d[3,j],marker='o',mfc=colors[q],mec='b') 
        ax21.plot(d[12, j],d[3,j],marker='o',mfc=colors[q],mec='r') 
[i.set_visible(False) for i in ax21.xaxis.get_ticklabels()[::2]]
ax21.set_xlabel('$\\epsilon_\\mathsf{x}$')
ax21.set_ylabel('$\\sigma_\\mathsf{x}$\n(ksi)')
ax21.axis(ymax=1.05*d[3].max())
f.ezlegend(ax21, loc='lower right', fontsize=16)
f.eztext(ax21, annotstr, 'ul')
f.myax(ax21)

ax22.plot(d[11],d[4], 'b')
for q,j in enumerate(incs):
    ax22.plot(d[11, j],d[4,j],marker='o',mfc=colors[q],mec=colors[q])
ax22.axis(ymax=1.05*d[4].max())
ax22.set_xlabel('$\\mathsf{e}_\\theta$')
ax22.set_ylabel('$\\sigma_\\theta$\n(ksi)')
f.myax(ax22)


fig1.savefig('F1.png',bbox_inches='tight')
fig2.savefig('F2.png',bbox_inches='tight')
n.savetxt('ProfileStages.dat',X=incs,fmt='%.0f',header='ProfileStages.  (genfromtxt reads this in as a 1D array)')

