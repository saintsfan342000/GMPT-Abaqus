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

save = 0
export = 1

try:
    case = argv[1]
except IndexError:
    case = 'BestSoFar'

try: 
    expdataplot = argv[2].upper()
except IndexError:
    expdataplot = 'Exp'

if case.upper() in ['BAD']:
    expts = [11, 4, 2, 8, 3, 12]
    jobno = [1 for i in range(len(expts))]
elif case.upper() in ['1PCT']:
    expts = [11, 4, 2, 8, 3, 12]
    jobno = [3, 3, 3, 7, 4, 3]
elif case.upper() in ['2PCT']:
    expts = [2,3,4,8,11,12]
    jobno = [4,5,4,9,4, 4]
elif case.upper() in ['H8']:
    expts = [11, 4, 2, 8, 3, 12, 10]
    jobno = [8888]*len(expts)
elif case.upper() in ['VM']:
    expts = [11, 4, 2, 8, 3, 12, 10]
    jobno = [2222]*len(expts)
else:
    case = 'BestSoFar'
    expts = [11, 4, 2, 8, 3, 12, 10]
    jobno = [2,  2, 2, 6, 3, 2,  5]

if export:
    import pandas as pd
    fname = 'GMPT-{}.xlsx'.format(case)
    fid = pd.ExcelWriter(fname)
    exphead = 'Expt AxSts, HoopSts, AxStn, HoopStn'.split(', ')
    simhead = 'Sim AxSts, HoopSts, AxStn, HoopStn'.split(', ')
    LL = n.empty((len(expts), 10))
    FA = n.empty((len(expts), 10))
    table_head = ('Expt, alpha, Rm, t, Ecc(%), SigX_LL, SigTH_LL, eX_LL, eQ_LL, ' + 
        'SigX_Fa, SigTh_Fa, e_Fa_Avg, eQ_Fa_Avg, eX_Fa_max, eQ_Fa_max').split()
    table = n.zeros((len(expts), len(table_head)))
    table = pd.DataFrame(table, index=expts, columns=table_head)
    
expts = n.array(expts)
jobno = n.array(jobno)   
    
key = n.genfromtxt('ExptSummary.dat', delimiter=',')
key = key[ n.in1d(key[:,0],expts) ]
key = key[ key[:,4].argsort() ]
for k,i in enumerate(expts):
    loc = n.nonzero(key[:,0]==i)[0][0]
    key[loc,1] = jobno[k]

expts, jobno = key[:,:2].astype(int).T
Radii, thicknesses, eccens = key[:, 6:].T
alphs = key[:,4]

repl = {'88':'11', '99':'8', '00':'8'}

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


for k, (exp,job, Rm, th, ecc, alp) in enumerate(zip(expts, 
                                    jobno, Radii, thicknesses, eccens, alphs)):

    print(case, expdataplot)
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

    # Simulation data
    cutoff = n.genfromtxt('{}/Cutoff.dat'.format(simpath),dtype=int, delimiter=',')[0]
    d = n.genfromtxt('{}/{}_NewResults.dat'.format(simpath, exp), delimiter=',')
    if cutoff > 0:
        d = d[:cutoff+1]
        print(exp, cutoff)
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
    
    # Expt nominal hoop stn in area close to ES_ANALZONE
    xstn= n.genfromtxt('{}/WholeFieldAverage.dat'.format(exppath), delimiter=',', usecols=(1,2))
    xstn = n.log(1+xstn)*100
    
    if expdataplot == 'CALDATA':
        #From Cal data!
        if k == 0: 
            case+='_CalData'
        # [0]Stage, [1]Wp, [2]SigX_Tru, [3]SigQ_True, [4]ex, [5]eq, [6]e3, 
        # [7]ep_x, [8]ep_q, [9]ep_3, [10]R_tru, [11]th_tru
        xstn = n.genfromtxt('{}/CalData.dat'.format(calpath), delimiter=',', usecols=(4,5))*100

    xex, xeq = xstn.T
    xst = n.genfromtxt('{}/STPF.dat'.format(exppath), delimiter=',', usecols=(range(6)))
    # [0]Stage, [1]Time, [2]Force(kip), [3]Pressure(ksi), [4]NomAxSts(ksi), [5]NomHoopSts(ksi)
    
    from scipy.signal import savgol_filter as sg
    winlen = xst.shape[0]//10
    if winlen%2 == 0: winlen+=1
    xst[:,5] = sg(xst[:,5], winlen, 1)
    xst[:,4] = sg(xst[:,4], winlen, 1)
    winlen//=2
    if winlen%2 == 0: winlen+=1
    xeq = sg(xeq, winlen, 1)
    xex = sg(xex, winlen, 1)
    xd[:,4] = sg(xd[:,4], winlen, 1) 
    
    exploc = xst[:,5].argmax()
    xu = n.genfromtxt('{}/ur_profiles.dat'.format(exppath), delimiter=',', usecols=(0,3*exploc+2,-1))
    xu[:,1:]*=100
    # Limitload y-coord, ur/Ro

    
    alpha, blank = n.polyfit(d[:,4],d[:,3],1)
    annotstr = '{}/{} || {:.2f}'.format(exp,job,alpha)

    # ax1: Ax sts vs ax stn.  
    # delta/L
    a, = ax1.plot(d[:,5],d[:,3], ':', label=annotstr)
    x,= ax1.plot(xd[:,4], xst[:,4], a.get_color(), alpha=0.5)
    ax1.plot(d[simloc,5],d[simloc,3], 'rD')
    ax1.plot(xd[exploc,4], xst[exploc,4], 'rD')
    #The point field averages
    #ax1.plot(d[:,11],d[:,2], ':', color=simcolor)
    #ax1.plot(n.log(xd[:,1]*.01+1)*100, xst[:,4], ':', color=expcolor)
    # Labels and annotations
    if k == len(expts) - 1:
        ax1.set_xlabel('$\\delta/\\mathsf{L}$ (%)')
        ax1.set_ylabel('$\\Sigma_\\mathsf{x}$\n(ksi)')
        ax1.axis(ymin=0, xmin=-3)
        f.ezlegend(ax1, title='Exp/#  ||  $\\eta\\prime$', loc='lower right')
        f.eztext(ax1, '$\\mathsf{L}_\\mathsf{g}$=1"', 'ul')
        f.myax(ax1)

        
    # ax2: Hoop sts vs Hoop stn
    # e_q
    ax2.plot(d[:,11],d[:,4], ':', label=annotstr)
    ax2.plot(xeq, xst[:,5], a.get_color(), alpha=0.5)
    ax2.plot(d[simloc,11], d[simloc,4], 'rD')
    ax2.plot(xeq[exploc], xst[exploc,5], 'rD')
    if k == len(expts) - 1:
        ax2.axis(xmin=0, xmax=10, ymin=0)
        ax2.set_xlabel('e$_\\theta$ (%)')
        ax2.set_ylabel('$\\Sigma_\\theta$\n(ksi)')
        #f.ezlegend(ax2, loc='upper right')
        f.eztext(ax2, 'Avg. of Pts.', 'br')
        f.myax(ax2)

    # ax3: Axial vs hoop stn
    ax3.plot(d[:,11], d[:,12], ':', label=annotstr)
    ax3.plot(xeq, xex, a.get_color(), alpha=0.5)
    ax3.plot(d[simloc,11], d[simloc, 12], 'rD')
    ax3.plot(xeq[exploc], xex[exploc], 'rD')
    if k == len(expts) - 1:
        ax3.axis([0,10,-4,6])
        ax3.set_xlabel('e$_\\theta$ (%)')
        ax3.set_ylabel('e$_\\mathsf{x}$\n(%)')
        #f.ezlegend(ax3, loc='upper right')
        f.myax(ax3)

    # Ax4: Ur Profiles
    ax4.plot(100*u[:,(simloc+2)]/u[:,1], u[:,0]*2/4, ':', label=annotstr)
    ax4.plot(xu[:,1], 2*xu[:,0]/4, a.get_color(), alpha=0.5)
    if k == len(expts) - 1:
        ax4.set_xlabel('u$_\\mathsf{r}$/R$_\\mathsf{o}$ (%)')
        ax4.set_ylabel('$\\frac{\\mathsf{2y}_\\mathsf{o}}{\\mathsf{L}_\\mathsf{g}}$')
        ax4.axis(xmin=0,ymin=-1,ymax=1)
        f.eztext(ax4, 'LL', 'br')
        f.myax(ax4, autoscale=.75, nudge=('left',.2,0))
    
    if export:
        expdata = n.c_[xst[:,4:6], xex, xeq]
        simdata = n.c_[d[:,3:5], d[:,12], d[:,11]]
        LL[k,:2] = exp, alpha
        LL[k,2:6] = expdata[exploc]
        LL[k,6:] = simdata[simloc]
        FA[k,:2] = exp, alpha
        FA[k, 2:6] = expdata[-1]
        FA[k, 6:] = simdata[-1]
        pd.DataFrame(expdata).to_excel(fid, sheet_name='GMPT-{} || {:.2f}'.format(exp,alpha), 
                    header=exphead, index=False)
        pd.DataFrame(simdata).to_excel(fid, sheet_name='GMPT-{} || {:.2f}'.format(exp,alpha), 
                    header=simhead, index=False, startcol=(expdata.shape[1]+2))

        xmax = n.genfromtxt('{}/MaxPt.dat'.format(exppath), delimiter=',')[-1]
        # [0]NEx [1]NEy [2]Gamma [3]F11-1 [4]F22-1 [5]atan(F12/F22) [6]epeq
        xmax = n.log(1+xmax)*100
        # [0]Expt, [1]alpha, [2]Rm, [3]t, [4]Ecc(%), [5]SigX_LL, [6]SigTH_LL, [7]eX_LL, [8]eQ_LL' + 
        # [9]SigX_Fa, [10]SigTh_Fa, [11]e_Fa_Avg, [12]eQ_Fa_Avg, [13]eX_Fa_max, [14]eQ_Fa_max'
        data = [exp, '{:g}'.format(alp), '{:.1f}\n{:.3f}'.format(Rm*25.4,Rm),
                '{:.1f}\n{:.3f}'.format(th*25.4,th), '{:.2f}'.format(ecc)]
        data.extend(['{:.0f}\n{:.1f}'.format(j*6.89655, j) for j in LL[k,2:4]])
        data.extend(['{:.2f}'.format(j) for j in LL[k,4:6]])
        data.extend(['{:.0f}\n{:.1f}'.format(j*6.89655, j) for j in FA[k,2:4]])
        data.extend(['{:.2f}'.format(j) for j in FA[k,4:6]])
        data.extend(['{:.2f}'.format(j) for j in xmax[:2][::-1]])
        table.loc[exp] = data

if export:
    
    pd.DataFrame(LL).to_excel(fid, sheet_name='LL', header=['Exp', 'Alpha']+exphead+simhead, index=False)
    pd.DataFrame(FA).to_excel(fid, sheet_name='Failure', header=['Exp', 'Alpha']+exphead+simhead, index=False)
    pd.DataFrame(table).to_excel(fid, sheet_name='Summary Table', index=False)
    fid.close()
              
fig.suptitle(case)

if save:
    fig.savefig('AllPlot_{}.jpg'.format(case), dpi=200, bbox_inches='tight')
    fig.savefig('AllPlot_{}.pdf'.format(case),bbox_inches='tight')

p.show('all')

if export:
    fig, ax1, ax2 = f.make12()
    l, = ax1.plot(expdata[:,3], expdata[:,1], alpha=0.5)
    ax1.plot(simdata[:,3], simdata[:,1], '--', color=l.get_color())
    ax1.plot(LL[-1,5], LL[-1, 3], 'rD')
    ax1.plot(LL[-1,5+4], LL[-1, 3+4], 'ro')
    l, = ax2.plot(expdata[:,2], expdata[:,0], alpha=0.5)
    ax2.plot(simdata[:,2], simdata[:,0], '--', color=l.get_color())
    ax2.plot(LL[-1,4], LL[-1, 2], 'rD')
    ax2.plot(LL[-1,4+4], LL[-1, 2+4], 'ro')    
    p.show()
