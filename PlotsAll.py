import numpy as n
from numpy import pi
import matplotlib.pyplot as p
from sys import argv
import figfun as f
p.style.use('mysty')
import glob

jobs = ['Eccen0p4','Eccen3p8']
paths = ['4_GMPT1','1_GMPT1']
label = ['0.4%','3.8%']
lc = ['chocolate','blueviolet']
ls=['-','--']

for k,i in enumerate(jobs):

    d = n.genfromtxt('../{}_{}/{}_results.dat'.format(paths[k],i,i), delimiter=',').T

    #[0] Force (kip), [1]Pressure (ksi), [2]NomAxSts, [3]NomHoopSts,' 
    #[4]d/Lg lo, [5]d/Lg hi, [6,7,8]S11,22,33, 
    # [9,10,11]LE11,22,33, [12]Vol'

    loc = n.argmax(d[1])

    d = d[:,:loc+1]
    Pmax = d[1,-1]
    Pspace = n.linspace(0,Pmax,5)[1:]
    incs = n.empty((len(Pspace)))
    for q,j in enumerate(Pspace):
        incs[q] = n.nonzero(d[1]>=j)[0][0]
    incs = incs.astype(int)

    # Profiles
    if k == 0:
        fig4 = p.figure()
        ax4 = fig4.add_subplot(111)
        colors = []
    u = n.genfromtxt('../{}_{}/{}_UR.dat'.format(paths[k],i,i),delimiter=',')
    incs = n.linspace(0,loc,5)[1:].astype(int)
    if i == jobs[0]:
        for j in incs:
            l, = ax4.plot(u[j+2],u[0],label=label[k])
            colors.append(l.get_color())
    else:
        for q,j in enumerate(incs):
            l2, =  ax4.plot(u[j+2],u[0],label=label[k],
                     linestyle='--',color=colors[q])
    if i == jobs[-1]:
        ax4.set_xlabel('U$_\\mathsf{r}$ (in)')
        ax4.set_ylabel('Z$_\\mathsf{o}$\n(in)')
        ax4.axis(ymax=2, xmin=0)
        leg = ax4.legend([l,l2],[L.get_label() for L in [l,l2]],
                   loc='upper right')
        p.setp(leg.get_lines(), color='k')
        f.myax(ax4)


    # Ax vs hoop sts
    if k == 0:
        fig1 = p.figure()
        ax1 = fig1.add_subplot(111)
    m,b = n.polyfit(d[7],d[8],1)
    for q,j in enumerate(incs):
        ax1.plot(d[7, j],d[8,j],marker='o',mfc=colors[q],mec=colors[q])
    ax1.plot(d[7], d[8], label='{:s}, {:.2f}'.format(label[k],m), zorder=-10,color=lc[k],ls=ls[k])
    if i == jobs[-1]:  
        ax1.set_xlabel('$\\sigma_\\theta$ (ksi)')
        ax1.set_ylabel('$\\sigma_\\mathsf{x}$\n(ksi)')
        ax1.legend(loc='lower right')
        f.myax(ax1)

    # Ax sts vs ax stn (logarithmic)
    if k == 0:
        fig2, ax21, ax22 = f.make12()
    for q,j in enumerate(incs):
        ax21.plot(d[11, j],d[2,j],marker='o',mfc=colors[q],mec=colors[q])   
    ax21.plot(d[11],d[2],label=label[k],zorder=-10,color=lc[k],ls=ls[k])
    if i == jobs[-1]:
        ax21.set_xlabel('$\\mathsf{e}_\\mathsf{x}$')
        ax21.set_ylabel('$\\sigma_x$\n(ksi)')
        ax21.legend(loc='lower right')
        [z.set_visible(False) for z in ax21.xaxis.get_ticklabels()[::2]]
        f.myax(ax21)
    for q,j in enumerate(incs):
        ax22.plot(d[10, j],d[3,j],marker='o',mfc=colors[q],mec=colors[q])
    ax22.plot(d[10],d[3],label=label[k],zorder=-10,color=lc[k],ls=ls[k])

    if i == jobs[-1]:
        ax22.set_xlabel('$\\mathsf{e}_\\theta$')
        ax22.set_ylabel('$\\sigma_\\theta$\n(ksi)')
        ax22.legend(loc='lower right')
        f.myax(ax22)

    # Ax sts vs ax stn (nominal)
    if k == 0:
        fig5, ax51, ax52 = f.make12()
    for q,j in enumerate(incs):
        ax51.plot(d[4, j],d[2,j],marker='o',mfc=colors[q],mec=colors[q])   
    ax51.plot(d[4],d[2],label=label[k],zorder=-10,color=lc[k],ls=ls[k])
    if i == jobs[-1]:
        ax51.set_xlabel('$\\epsilon_\\mathsf{x}$')
        ax51.set_ylabel('$\\sigma_x$\n(ksi)')
        ax51.legend(loc='lower right')
        [z.set_visible(False) for z in ax51.xaxis.get_ticklabels()[::2]]
        f.myax(ax51)
    for q,j in enumerate(incs):
        ax52.plot(u[j+2, 0]/u[1,0],d[3,j],marker='o',mfc=colors[q],mec=colors[q])
    ax52.plot(u[2:loc+1+2,0]/u[1,0],d[3],label=label[k],zorder=-10,color=lc[k],ls=ls[k])
    if i == jobs[-1]:
        ax52.set_xlabel('$\\epsilon_\\theta$')
        ax52.set_ylabel('$\\sigma_\\theta$\n(ksi)')
        ax52.legend(loc='lower right')
        f.myax(ax52)
    # Axial vs hoop stn
    if k == 0:
        fig3 = p.figure()
        ax3 = fig3.add_subplot(111)
    for q,j in enumerate(incs):
        ax3.plot(d[10, j],d[11,j],marker='o',mfc=colors[q],mec=colors[q])      
    ax3.plot(d[10], d[11], label=label[k], zorder=-10,color=lc[k],ls=ls[k])

    if i == jobs[-1]:
        ax3.set_ylabel('$\\mathsf{e}_\\mathsf{x}$')
        ax3.set_xlabel('$\\mathsf{e}_\\theta$')
        ax3.legend(loc='lower right')
        ax3.axis('equal')
        f.myax(ax3)

    print(i, loc+1)

fig1.savefig('F1.png',bbox_inches='tight')
fig2.savefig('F2.png',bbox_inches='tight')
fig3.savefig('F3.png',bbox_inches='tight')
fig4.savefig('F4.png',bbox_inches='tight')
fig5.savefig('F5.png',bbox_inches='tight')

p.show()

