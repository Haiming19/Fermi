#!/usr/bin/env python
from UpperLimits import UpperLimits
#Import all of the needed modules
import pyLikelihood as pyLike
#import UnbinnedAnalysis as ua
import BinnedAnalysis as ba
import numpy as np
from os import path
from gt_apps import filter, maketime, expCube, expMap, diffResps, evtbin, srcMaps, gtexpcube2
#import matplotlib.pyplot as plt
#import pyfits
#from likeSED import *


sourcename="3C138"
#time_start = 410629417. #time(s)
#time_end = 426181417.1
#ra = 202.6356      # ra (degree)
#dec = 30.3886      # dec (degree)
#evclass = 128
#evtype = 3         #front(1)+black(2)
#rad = 10.          #ROI (degree) 10x10
#emin = 100.        #Emin unit in MeV
#emax = 300000.     #Emax unit in MeV
#zmax = 90          #maximum zenith angle (degree)
#irfs = "CALDB"
#SC = "SC.fits"     #spacecraft file 
#infile = "@events.list"    #event file 
#model = "./model.xml"  #xml file

#print "**Working on the %d th bin: range (%f,%f)" % (i,time[0],time[1])
##----------------------------------------------------------------------------------------------    
#print 'Step1:***Running gtselect***'

#filter['evclass'] = evclass
#filter['evtype'] = evtype
#filter['ra'] = ra
#filter['dec'] = dec
#filter['rad'] = rad
#filter['emin'] = emin
#filter['emax'] = emax
#filter['zmax'] = zmax
#filter['tmin'] = time_start
#filter['tmax'] = time_end
#filter['infile'] = infile
#filter['outfile'] = sourcename+'_filtered.fits'
#filter.run()
##----------------------------------------------------------------------------------------------
#print 'Step2:***Running gtmktime***'

#maketime['scfile'] = SC
#maketime['filter'] = '(DATA_QUAL==1)&&(LAT_CONFIG==1)'
#maketime['roicut'] = 'yes'
#maketime['evfile'] = sourcename+'_filtered.fits'
#maketime['outfile'] = sourcename+'_filtered_gti.fits'
#maketime.run()
            
##----------------------------------------------------------------------------------------------
#print 'Step3:***Running gtbin***'

#evtbin['algorithm'] = 'CCUBE'
#evtbin['evfile'] = sourename+'_filtered_gti.fits'
#evtbin['outfile'] = sourcename+'_CCUBE.fits'
#evtbin['scfile'] = 'NONE'
#evtbin['nxpix'] = 140
#evtbin['nypix'] = 140
#evtbin['binsz'] = 0.1            
#evtbin['coordsys'] = 'CEL'
#evtbin['xref'] = ra
#evtbin['yref'] = dec
#evtbin['axisrot'] = 0.0
#evtbin['proj'] = 'AIT'
#evtbin['ebinalg'] = 'LOG'
#evtbin['emin'] = emin
#evtbin['emax'] = emax
#evtbin['enumbins'] = 30            
#evtbin.run()
                       
##----------------------------------------------------------------------------------------------            
#print 'Step4:***Running gtltcube***'

#expCube['evfile'] = sourcename+'_filtered_gti.fits'
#expCube['scfile'] = SC
#expCube['outfile'] = sourcename+'_expCube.fits'
#expCube['dcostheta'] = 0.025
#expCube['binsz'] = 1
#expCube.run()
           
##----------------------------------------------------------------------------------------------  
#print 'Step5:***Running gtexpcube2***'

#gtexpcube2['infile'] = sourcename+'_expCube.fits'
#gtexpcube2['outfile'] = sourcename+'_BinnedExpMap.fits'
#gtexpcube2['cmap'] = 'none'
#gtexpcube2['irfs'] = 'P8R2_SOURCE_V6'
#gtexpcube2['nxpix'] = 400
#gtexpcube2['nypix'] = 400
#gtexpcube2['binsz'] = 0.1            
#gtexpcube2['coordsys'] = 'CEL'
#gtexpcube2['xref'] = ra
#gtexpcube2['yref'] = dec
#gtexpcube2['axisrot'] = 0.0
#gtexpcube2['proj'] = 'AIT'
#gtexpcube2['emin'] = emin
#gtexpcube2['emax'] = emax
#gtexpcube2['enumbins'] = 30            
#gtexpcube2.run()
            
##----------------------------------------------------------------------------------------------  
#print 'Step6:***Running gtsrcmaps***'
model = "./model.xml"
#srcMaps['expcube'] = sourcename+'_expCube.fits'
#srcMaps['cmap'] = sourcename+'_CCUBE.fits'
#srcMaps['bexpmap'] = sourcename+'_BinnedExpMap.fits' 
#srcMaps['srcmdl'] = model
#srcMaps['outfile'] = sourcename+'_srcMaps.fits' 
#srcMaps['irfs'] = 'CALDB'
#srcMaps.run()

#----------------------------------------------------------------------------------------------              
print '***Running likelihood analysis***'
#---------------------------------------------------------------------
#model = "./model.xml"

obs = ba.BinnedObs(srcMaps=sourcename+'_srcMaps.fits',
                expCube=sourcename+'_expCube.fits',
                binnedExpMap=sourcename+'_BinnedExpMap.fits',
                irfs='CALDB')


like = ba.BinnedAnalysis(obs,model,optimizer='MINUIT')
like.tol = 0.0001
likeObj = pyLike.Minuit(like.logLike)
print "#------------------------------------------------------#"
loglikelihood= like.fit(verbosity=0,covar=True,optObject=likeObj)
print loglikelihood

#print likeobj.getRetCode()


#like = ba.BinnedAnalysis(obs,model,optimizer='NewMINUIT')
#like.tol = 0.0001
#likeobj = pyLike.NewMinuit(like.logLike)
#print "#------------------------------------------------------#"
#print like.fit(verbosity=0,covar=True,optObject=likeobj)
#print "#------------------------------------------------------#"

#like.logLike.writeXml('Fittednone.xml')  

name="4FGL J0521.2+1637"

TS=like.Ts(name)


flux1=like.energyFlux(name,emin=100,emax=300000)
fluxerr1=like.energyFluxError(name,emin=100,emax=300000)
flux=like.flux(name,emin=100,emax=300000)
fluxerr=like.fluxError(name,emin=100,emax=300000)



print loglikelihood,TS,flux,fluxerr,flux1,fluxerr1

#f=open("TSresult_allnone.txt","w+")
#f1=open("sname.txt","w+")
#sourceDetails = {}
#for source in like.sourceNames():
    #sourceDetails[source] = like.Ts(source)
    #print >>f1,source
#for source,TS in sourceDetails.iteritems():
    #if (TS < 16):
        #like.deleteSource(source)
    #else:
    #print >>f,source, TS
#like.logLike.writeXml('modelfitted.xml')      
#like.writeXml('modelfitted2.xml')
like.logLike.writeXml('fitted.xml') 


#------------------
import matplotlib.pyplot as plt
N0 = like.model[name].funcs['Spectrum'].getParam('Prefactor').value()
N0_err = like.model[name].funcs['Spectrum'].getParam('Prefactor').error()
scale = 1e-13
N0 = N0*scale
N0_err = N0_err*scale
gamma = like.model[name].funcs['Spectrum'].getParam('Index').value()
gamma_err = like.model[name].funcs['Spectrum'].getParam('Index').error()
E0 = like.model[name].funcs['Spectrum'].getParam('Scale').value()

freeParValues = []
for sourcename in like.sourceNames():
    for element in like.freePars(sourcename):
        freeParValues.append(element.getValue())

g_index = freeParValues.index(like.freePars(name)[1].getValue())
cov_gg = like.covariance[g_index][g_index]

print N0,N0_err,gamma,gamma_err,E0,cov_gg,g_index,scale

f = lambda E,N0,E0,gamma: N0*(E/E0)**(-1*gamma)
ferr = lambda E,F,N0,N0err,E0,cov_gg: F*np.sqrt(N0err**2/N0**2 + ((np.log(E/E0))**2)*cov_gg)

E = np.logspace(2,5,100)
F = f(E,N0,E0,gamma)
Ferr = ferr(E,F,N0,N0_err,E0,cov_gg)

plt.figure(figsize=(8,8))
plt.xlabel('Energy [MeV]')
plt.ylabel(r'E$^2$ dN/dE [MeV s$^{-1}$ cm$^{-2}$]')
plt.loglog(E,1.6022e-6*E**2*(F+Ferr))
plt.loglog(E,1.6022e-6*E**2*(F-Ferr))
plt.plot(E,1.622e-6*E**2*F)
plt.show()


#-----------------------------------------

#np.savez('sed.npz',minEs=sed.likeIn.bins[0],maxEs=sed.likeIn.bins[1],ecent=sed.centers,fluxPts=sed.data[0],fluxErrs=sed.data[1],tsbands=sed.data[2],gammas=sed.data[3])


#------------------

#---------------------------------------
#ul=UpperLimits(like)
#ul[name].compute(emin=10000,emax=500000)
#UL=ul[name].results
#print UL

#f.close()
#f1.close()


#f1=open("ULresult_2d.txt","w+")
##model1="./model02.xml"
##if (like.TS(name) <9):
    ##like1 = ba.BinnedAnalysis(obs,model1,optimizer='MINUIT')
#ul=UpperLimits(like)
#ul['Geminga'].compute(emin=10000,emax=500000)
#UL=ul['Geminga'].results
#print UL
#print >>f1,UL

    #if (TS < 9):
        #print "Deleting..."
        #like.deleteSource(source)

#like.fit(verbosity=0,covar=True,optObject=likeobj)

#print likeobj.getRetCode()

    
#like.logLike.writeXml('modelfitted.xml')  
#f1.close()            
print '***Finished***'
