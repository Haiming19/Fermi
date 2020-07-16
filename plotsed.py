import numpy as np
import matplotlib.pyplot as plt

srcName="4FGL J0521.2+1637"  #4FGL name


#------------------------------------------------------------------------------------------
N0 = 7.46213653926e-14
N0_err = 2.17542004251e-15
scale = 1e-13
gamma = 2.18929599463
gamma_err = 0.0212860988071
E0 = 1563.118286
cov_gg = 0.000453098002425
#N0 = N0*scale
# = N0_err*scale


filename = 'bestfitted_global.txt'
f = np.load('sed.npz')

minEs=f['minEs']/1.e3
maxEs=f['maxEs']/1.e3
ecent=f['ecent']
fluxPts=f['fluxPts']*1.e3
fluxErrs=f['fluxErrs']*1.e3
tsbands=f['tsbands']
gammas=f['gammas']

ts,flux,fluxerr,emin,emax,ul,ecent = tsbands,fluxPts,fluxErrs,minEs,maxEs,fluxPts,ecent
slt = (fluxErrs != 0)

xlabel = 'Energy [GeV]'
ylabel = r'$\rm E^2{dN/dE}\ [erg\ cm^{-2}\ s^{-1}]$'

emid = np.sqrt(emin*emax)
edelt = emax - emin

#----------------------------
for i in range(len(ecent)):
    vfv=flux[i]*ecent[i]**2*0.0016022
    err=fluxerr[i]*ecent[i]**2*0.0016022
    UL=ul[i]*ecent[i]**2*0.0016022
    print ecent[i],emin[i],emax[i],vfv,err

#----------------------------
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.tick_params(top='on', right='on', which='both')

#esqf=[y*x**(2-g)*(g-1)/((e/1000.)**(1-g)-(E/1000.)**(1-g))*0.001602]
#esqfErr=[z*x**(2-g)*(g-1)/((e/1000.)**(1-g)-(E/1000.)**(1-g))*0.001602]

#plt.errorbar(emid[slt],flux[slt]*emid[slt]**(2-gamma)*(gamma-1)/(emin[slt]**(1-gamma)-emax[slt]**(1-gamma))*0.0016022, \
    #yerr=fluxerr[slt]*emid[slt]*0.0016022,xerr=(emid[slt]-emin[slt],emax[slt]-emid[slt]),fmt='k.')

plt.errorbar(emid[slt],flux[slt]*ecent[slt]**2*0.0016022, \
    yerr=fluxerr[slt]*ecent[slt]**2*0.0016022,xerr=(emid[slt]-emin[slt],emax[slt]-emid[slt]),fmt='k.')



if sum(~slt)>0:
    uul = ul[~slt]*ecent[~slt]**2*0.0016022
    plt.errorbar(emid[~slt],uul,xerr=(emid[~slt]-emin[~slt],emax[~slt]-emid[~slt]), \
        yerr=(-10**(np.log10(uul)-0.2)+uul,np.zeros_like(ul[~slt])),capsize=3, fmt='k.',uplims=1)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.xscale('log')
plt.yscale('log')
#------------------------------------------------------------------------------------------

f = lambda E,N0,E0,gamma: N0*(E/E0)**(-1*gamma)
ferr = lambda E,F,N0,N0err,E0,cov_gg: F*np.sqrt(N0err**2/N0**2 + ((np.log(E/E0))**2)*cov_gg)

E = np.logspace(2,5.477121254719663,100)
F = f(E,N0,E0,gamma)
Ferr = ferr(E,F,N0,N0_err,E0,cov_gg)

vfv=1.622e-6*E**2*F
vfverr1=1.6022e-6*E**2*(F+Ferr)
vfverr2=1.6022e-6*E**2*(F-Ferr)
E1 = E/1.e3
#plt.figure(figsize=(8,8))
#plt.xlabel('Energy [MeV]')
#plt.ylabel(r'E$^2$ dN/dE [MeV s$^{-1}$ cm$^{-2}$]')
plt.fill_between(E1,vfverr1, vfverr2, facecolor='grey', alpha=0.3)
#plt.loglog(E1,vfverr1)
#plt.loglog(E1,vfverr2)
plt.plot(E1,vfv,"red")

'''def pl(ee,p0,gm,e0):
    return p0*(ee/e0)**-gm

def logp(ee,n0,al,bt,eb):
    return n0*(ee/eb)**(-(al+bt*np.log(ee/eb)))

fitresult = eval(open(filename,'r').read())[srcName]
x = np.logspace(np.log10(minEs[0]),np.log10(maxEs[-1]),200)
sptype = fitresult["SpectrumType"]

from scipy.integrate import simps

if sptype == 'LogParabola':
    norm_  = np.float(fitresult["norm"].split('+/-')[0])
    alpha_ = np.float(fitresult["alpha"].split('+/-')[0])
    beta_  = np.float(fitresult["beta"].split('+/-')[0])
    Eb_    = np.float(fitresult["Eb"].split('+/-')[0])
    Flux_  = np.float(fitresult["Flux"].split('+/-')[0])
    func_ = lambda x: logp(x,norm_,alpha_,beta_,Eb_)
    PrefactorScale = Flux_/simps(func_(x*1.e3),x*1.e3)
    y = logp(x,norm_*PrefactorScale,alpha_,beta_,Eb_/1.e3)*1.e3 
elif sptype == 'PowerLaw':
    Prefactor_ = np.float(fitresult["Prefactor"].split('+/-')[0])*PrefactorScale
    Index_ = np.float(fitresult["Index"].split('+/-')[0])
    Scale_ = np.float(fitresult["Scale"].split('+/-')[0])
    Flux_  = np.float(fitresult["Flux"].split('+/-')[0])
    func_ = lambda x: pl(x,Prefactor_,Index_,Scale_)
    PrefactorScale = Flux_/simps(func_(x*1.e3),x*1.e3)
    y = pl(x,Prefactor_,Index_,Scale_/1.e3)*1.e3
plt.plot(x,y*x**2*0.0016021766208)'''

plt.ylim(5e-15,1e-11)

#font1 = {'size': 22} 
ax=plt.gca()
ax.spines['bottom'].set_linewidth(1.3)
ax.spines['left'].set_linewidth(1.3)
ax.spines['right'].set_linewidth(1.3)
ax.spines['top'].set_linewidth(1.3)
#plt.tick_params(labelsize=20)
plt.text(1e-1, 1.5e-14, "%s" %srcName, size = 10, alpha = 0.01)



plt.savefig(srcName.replace(' ','_')+'_sed0.eps')
plt.show()
