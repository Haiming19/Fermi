#!/usr/bin/env python
# coding: utf-8

#Import all of the needed modules
from UpperLimits import UpperLimits
from bdlikeSED import *
import pyLikelihood as pyLike
import BinnedAnalysis as ba
import numpy as np
import gt_apps as my_apps
from os import path
from gt_apps import filter, maketime, expCube, expMap, diffResps, evtbin, srcMaps, gtexpcube2
from make4FGLxml import *
#from  gtltcube_mp import *

sourcename="4C+55.17" 

sourcename1="4FGL J0957.6+5523"  #4FGL name

time_start = 239557417. #time(s)
time_end = 584235367.  # years
ra = 149.409     # ra (degree)
dec = 55.3827      # dec (degree)
evclass = 128
evtype = 3         #=front(1)+black(2)
rad = 15.          #ROI (degree) 30
emin = 100       #Emin unit in MeV
emax = 300000     #Emax unit in MeV
zmax = 90          #maximum zenith angle (degree)
irfs = "P8R3_SOURCE_V2"
SC = "/home/ganyingying/AGN/3FGL-J0957.6+5523/L19070722455808DBD12D42_SC00.fits"     #spacecraft file 
infile = "@events.list"    #event file 

#----------------------------------------------------------------------------------------------    
'''print 'Step1:***Running gtselect***'

filter['evclass'] = evclass
filter['evtype'] = evtype
filter['ra'] = ra
filter['dec'] = dec
filter['rad'] = rad
filter['emin'] = emin
filter['emax'] = emax
filter['zmax'] = zmax
filter['tmin'] = time_start
filter['tmax'] = time_end
filter['infile'] = infile
filter['outfile'] = sourcename+'_filtered_.fits'
filter.run()
###----------------------------------------------------------------------------------------------
print 'Step2:***Running gtmktime***'

maketime['scfile'] = SC   #spacecraft file 
maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
maketime['roicut'] = 'no'
maketime['evfile'] = sourcename+'_filtered_.fits'
maketime['outfile'] = sourcename+'_filtered_gti.fits'
maketime.run()

#####-----------------------Create a source model XML file 4FGL--------------------------------------------------
mymodel = srcList('gll_psc_v17.fit',sourcename+'_filtered_gti.fits','model.xml')  
mymodel.makeModel('gll_iem_v07.fits', 'gll_iem_v07', 'iso_P8R3_SOURCE_V2_v1.txt', 'iso_P8R3_SOURCE_V2_v1')
###----------------------------------------------------------------------------------------------
print 'Step3:***Running gtbin***'

evtbin['algorithm'] = 'CCUBE'
evtbin['evfile'] = sourcename+'_filtered_gti.fits'
evtbin['outfile'] = sourcename+'_CCUBE.fits'
evtbin['scfile'] = 'NONE'
evtbin['nxpix'] = 100
evtbin['nypix'] = 100
evtbin['binsz'] = 0.2            
evtbin['coordsys'] = 'CEL'
evtbin['xref'] = ra
evtbin['yref'] = dec
evtbin['axisrot'] = 0.0
evtbin['proj'] = 'AIT'
evtbin['ebinalg'] = 'LOG'
evtbin['emin'] = emin
evtbin['emax'] = emax
evtbin['enumbins'] = 37            
evtbin.run()
####---------------------------------------------------------------------------------------------- 
print 'Step4:***Running gtltcube***'
#---------RUN multi-threaded 10 CORES---------------------
#gtltcube_mp(bins=10, SCFile=SC, EVFile=sourcename+'_filtered_gti.fits', OutFile=sourcename+'_expCube.fits',SaveTemp=False,zmax=zmax)   
#---------10 CORES---------------------
expCube['zmax']=zmax
expCube['evfile'] = sourcename+'_filtered_gti.fits'
expCube['scfile'] = SC
expCube['outfile'] = sourcename+'_expCube.fits'
expCube['dcostheta'] = 0.025
expCube['binsz'] = 1
expCube.run()
           
#####----------------------------------------------------------------------------------------------  
print 'Step5:***Running gtexpcube2***'

gtexpcube2['infile'] = sourcename+'_expCube.fits'
gtexpcube2['outfile'] = sourcename+'_BinnedExpMap.fits'
gtexpcube2['cmap'] = 'none'
gtexpcube2['irfs'] = irfs
gtexpcube2['nxpix'] = 300
gtexpcube2['nypix'] = 300
gtexpcube2['binsz'] = 0.2           
gtexpcube2['coordsys'] = 'CEL'
gtexpcube2['xref'] = ra
gtexpcube2['yref'] = dec
gtexpcube2['axisrot'] = 0.0
gtexpcube2['proj'] = 'AIT'
gtexpcube2['emin'] = emin
gtexpcube2['emax'] = emax
gtexpcube2['enumbins'] = 37            
gtexpcube2.run()
            
###----------------------------------------------------------------------------------------------  
print 'Step6:***Running gtsrcmaps***'
model = "./model.xml"  #xml file 
srcMaps['expcube'] = sourcename+'_expCube.fits'
srcMaps['cmap'] = sourcename+'_CCUBE.fits'
srcMaps['bexpmap'] = sourcename+'_BinnedExpMap.fits' 
srcMaps['srcmdl'] = model
srcMaps['outfile'] = sourcename+'_srcMaps.fits' 
srcMaps['irfs'] = 'CALDB'
srcMaps.run()'''

#----------------------------------------------------------------------------------------------              
print('***Running likelihood analysis***')
#---------------------------------------------------------------------
model = "./model.xml"
#model = "./modelfitted2.xml"
obs = ba.BinnedObs(srcMaps=sourcename+'_srcMaps.fits',
                expCube=sourcename+'_expCube.fits',
                binnedExpMap=sourcename+'_BinnedExpMap.fits',
                irfs='CALDB')

#---------------第一次整体拟合目的在于删除一些弱源---------------
like = ba.BinnedAnalysis(obs,model,optimizer='MINUIT')
like.tol = 1e-8
likeObj = pyLike.Minuit(like.logLike)
like.fit(verbosity=0,covar=True,optObject=likeObj)
sourceDetails = {}
for source in like.sourceNames():
    sourceDetails[source] = like.Ts(source)

for source,TS in sourceDetails.iteritems():
    if (TS < 16):
        like.deleteSource(source)
        
like.logLike.writeXml('modelfitted.xml')      


model1="./modelfitted.xml"
#----------------第二次整体拟合---------------------------
like1 = ba.BinnedAnalysis(obs,model1,optimizer='MINUIT')
like1.tol = 1e-8
likeObj1 = pyLike.Minuit(like1.logLike)
like1.fit(verbosity=0,covar=True,optObject=likeObj1)
like1.logLike.writeXml('bestfitted_global.xml')

from mytools import printResult,sp2pl
printResult(like1,saveresultto='bestfitted_global.txt')
sp2pl('bestfitted_global.xml','sed_in.xml',sourcename1)

assert 1
#----------------产生SED---------------------------
like2 = ba.BinnedAnalysis(obs,'sed_in.xml',optimizer='MINUIT')
like2.tol = 1e-8
#inputs=bdlikeInput(like2,sourcename+'_filtered_gti.fits',sourcename+'_CCUBE.fits',"sc.fits",sourcename1)
inputs=bdlikeInput(like2,sourcename+'_filtered_gti.fits',sourcename+'_CCUBE.fits',"/home/ganyingying/AGN/3FGL-J0957.6+5523/L19070722455808DBD12D42_SC00.fits",sourcename1)
inputs.plotBins()
inputs.fullFit()
sed=bdlikeSED(inputs)
sed.getECent()
sed.fitBands()
#sed.Plot()
#-----------------------------------------

print('***Finished***')

np.savez('sed.npz',minEs=sed.likeIn.bins[0],maxEs=sed.likeIn.bins[1],ecent=sed.centers,fluxPts=sed.data[0],fluxErrs=sed.data[1],tsbands=sed.data[2],gammas=sed.data[3])

