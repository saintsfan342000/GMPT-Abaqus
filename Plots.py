import numpy as n
from numpy import pi
import matplotlib.pyplot as p
from sys import argv
import figfun as f
p.style.use('mysty')
import glob

name = glob.glob('*_results.dat')[0].split('_')[0]

d = n.genfromtxt('{}_results.dat'.format(name), delimiter=',').T
#[0] Force (kip), [1]Pressure (ksi), [2]NomAxSts, [3]NomHoopSts,' 
#[4]d/Lg lo, [5]d/Lg hi, [6,7,8]S11,22,33, 
# [9,10,11]LE11,22,33, [12]Vol'
u = n.genfromtxt('{}_UR.dat'.format(name.split('_')[0]),delimiter=',')

loc = n.argmax(d[1])
d = d[:,:loc+1]

Pmax = d[1,-1]
Pspace = n.linspace(0,Pmax,5)[1:]
incs = n.empty((len(Pspace)))
for q,j in enumerate(Pspace):
    incs[q] = n.nonzero(d[1]>=j)[0][0]
incs = incs.astype(int)

incs = n.linspace(0,loc,5)[1:].astype(int)

fig4 = p.figure()
ax4 = fig4.add_subplot(111)
colors = []
for j in incs:
    l, = ax4.plot(u[j+2],u[0])
    colors.append(l.get_color())
ax4.set_xlabel('U$_\\mathsf{r}$ (in)')
ax4.set_ylabel('Z$_\\mathsf{o}$\n(in)')
ax4.axis(ymax=2, xmin=0)
f.myax(ax4)

# Ax vs hoop sts
fig1 = p.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(d[3], d[2], 'b',label='Nominal')
ax1.plot(d[7], d[8], 'r', label='True')
for q,j in enumerate(incs):
    ax1.plot(d[3, j],d[2,j],marker='o',mfc=colors[q],mec='b')
    ax1.plot(d[7, j],d[8,j],marker='o',mfc=colors[q],mec='r')
ax1.set_xlabel('$\\sigma_\\theta$ (ksi)')
ax1.set_ylabel('$\\sigma_\\mathsf{x}$\n(ksi)')
f.ezlegend(ax1, loc='lower right')
f.myax(ax1)

# Ax sts vs ax stn
fig2, ax21, ax22 = f.make12()
ax21.plot(d[4],d[2],'b', label='$\\delta/\\mathsf{L}$')
ax21.plot(d[11],d[2],'r', label='$e_\\mathsf{x}$')
for q,j in enumerate(incs):
        ax21.plot(d[4, j],d[2,j],marker='o',mfc=colors[q],mec='b') 
        ax21.plot(d[11, j],d[2,j],marker='o',mfc=colors[q],mec='r') 
[i.set_visible(False) for i in ax21.xaxis.get_ticklabels()[::2]]
ax21.set_xlabel('$\\epsilon_\\mathsf{x}$')
ax21.set_ylabel('$\\sigma_\\mathsf{x}$\n(ksi)')
f.ezlegend(ax21, loc='lower right', fontsize=16)
f.myax(ax21)

ax22.plot(d[10],d[3], 'b')
for q,j in enumerate(incs):
    ax22.plot(d[10, j],d[3,j],marker='o',mfc=colors[q],mec=colors[q])
ax22.set_xlabel('$\\mathsf{e}_\\theta$')
ax22.set_ylabel('$\\sigma_\\theta$\n(ksi)')
f.myax(ax22)

# Axial vs hoop stn
fig3 = p.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(d[10], d[11], 'b')
for q,j in enumerate(incs):
    ax3.plot(d[10, j],d[11,j],marker='o',mfc=colors[q],mec=colors[q])  
ax3.set_xlabel('$\\mathsf{e}_\\theta$')
ax3.set_ylabel('$\\mathsf{e}_\\mathsf{x}$')
ax3.axis('equal')
f.myax(ax3)

fig1.savefig('F1.png',bbox_inches='tight')
fig2.savefig('F2.png',bbox_inches='tight')
fig3.savefig('F3.png',bbox_inches='tight')
fig4.savefig('F4.png',bbox_inches='tight')

p.show()

