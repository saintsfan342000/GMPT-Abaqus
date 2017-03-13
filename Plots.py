import numpy as n
from numpy import pi
import matplotlib.pyplot as p
from sys import argv
import figfun as f
p.style.use('mysty')
import glob

name = glob.glob('*_results.dat')[0]

d = n.genfromtxt('{}'.format(name), delimiter=',').T
#[0] Force (kip), [1]Pressure (ksi), [2]NomAxSts, [3]NomHoopSts,' 
#[4]d/Lg lo, [5]d/Lg hi, [6,7,8]S11,22,33, 
# [9,10,11]LE11,22,33, [12]Vol'

loc = n.argmax(d[1])

d = d[:,:loc+1]

# Ax vs hoop sts
fig1 = p.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(d[3], d[2], label='Nominal')
ax1.plot(d[7], d[8], label='True')
ax1.set_xlabel('$\\sigma_\\theta$ (ksi)')
ax1.set_ylabel('$\\sigma_\\mathsf{x}$\n(ksi)')
ax1.legend(loc='lower right')
f.myax(ax1)

# Ax sts vs ax stn
fig2, ax21, ax22 = f.make12()
ax21.plot(d[4],d[2],label='$\\delta/\\mathsf{L}$')
ax21.plot(d[11],d[2],label='$e_\\mathsf{x}$')
[i.set_visible(False) for i in ax21.xaxis.get_ticklabels()[::2]]
ax21.set_xlabel('$\\epsilon_\\mathsf{x}$')
ax21.set_ylabel('$\\sigma_\\mathsf{x}$\n(ksi)')
ax21.legend(loc='lower right')
f.myax(ax21)

ax22.plot(d[10],d[3])
ax22.set_xlabel('$\\mathsf{e}_\\theta$')
ax22.set_ylabel('$\\sigma_\\theta$\n(ksi)')
f.myax(ax22)

# Axial vs hoop stn
fig3 = p.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(d[10], d[11])
ax3.set_xlabel('$\\mathsf{e}_\\theta$')
ax3.set_ylabel('$\\mathsf{e}_\\mathsf{x}$')
ax3.axis('equal')
f.myax(ax3)

fig1.savefig('F1.png',bbox_inches='tight')
fig2.savefig('F2.png',bbox_inches='tight')
fig3.savefig('F3.png',bbox_inches='tight')

p.show()

