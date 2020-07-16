#!/usr/bin/env python

#Import all of the needed modules

import numpy as np
from os import path
from gt_apps import filter, maketime, expCube, expMap, diffResps, evtbin, srcMaps, gtexpcube2
#from make3FGLxml import *
from make4FGLxml import *
#from make3FHLxml import *
from gtltcube_mp import *
sourcename="3C138"
time_start = 239557417. #time(s)
time_end = 584949858.  #8 years
ra = 80.2912    # ra (degree)
dec = 16.6395     # dec (degree)
evclass = 128
evtype = 3         #=front(1)+black(2)
rad = 10.          #ROI (degree) 30
emin = 100.        #Emin unit in MeV
emax = 300000.     #Emax unit in MeV
zmax = 90          #maximum zenith angle (degree)
irfs = "P8R3_SOURCE_V2"
SC = "SC.fits"     #spacecraft file 
infile = "@events.list"    #event file 

##----------------------------------------------------------------------------------------------    
print 'Step1:***Running gtselect***'

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
filter['outfile'] = sourcename+'_filtered.fits'
filter.run()
###----------------------------------------------------------------------------------------------
print 'Step2:***Running gtmktime***'

maketime['scfile'] = SC
maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
maketime['roicut'] = 'no'
maketime['evfile'] = sourcename+'_filtered.fits'
maketime['outfile'] = sourcename+'_filtered_gti.fits'
maketime.run()
#####-----------------------Create a source model XML file--------------------------------------------------
##mymodel = srcList('gll_psc_v16.fit',sourcename+'_filtered_gti.fits','model_3FGL.xml')  
##mymodel.makeModel('gll_iem_v06.fits', 'gll_iem_v06', 'iso_P8R3_SOURCE_V2.txt', 'iso_P8R3_SOURCE_V2')
#####-----------------------Create a source model XML file 4FGL--------------------------------------------------
mymodel = srcList('gll_psc_v19.fit',sourcename+'_filtered_gti.fits','model.xml')  
mymodel.makeModel('gll_iem_v07.fits', 'gll_iem_v07', 'iso_P8R3_SOURCE_V2_v1.txt', 'iso_P8R3_SOURCE_V2_v1')
###-----------------------Create a source model XML file 3FHL--------------------------------------------------
##mymodel = srcList('gll_psch_v13.fit',sourcename+'_filtered_gti.fits','model_3FHL.xml')  
##mymodel.makeModel('gll_iem_v07.fits', 'gll_iem_v07', 'iso_P8R3_SOURCE_V2_v1.txt', 'iso_P8R3_SOURCE_V2_v1')
###-------------------------------------------------------------------------                                  
###----------------------------------------------------------------------------------------------
print 'Step3:***Running gtbin***'

evtbin['algorithm'] = 'CCUBE'
evtbin['evfile'] = sourcename+'_filtered_gti.fits'
evtbin['outfile'] = sourcename+'_CCUBE.fits'
evtbin['scfile'] = 'NONE'
evtbin['nxpix'] = 140
evtbin['nypix'] = 140
evtbin['binsz'] = 0.1            
evtbin['coordsys'] = 'CEL'
evtbin['xref'] = ra
evtbin['yref'] = dec
evtbin['axisrot'] = 0.0
evtbin['proj'] = 'AIT'
evtbin['ebinalg'] = 'LOG'
evtbin['emin'] = emin
evtbin['emax'] = emax
evtbin['enumbins'] = 12            
evtbin.run()
###----------------------------------------------------------------------------------------------            
print 'Step4:***Running gtltcube***'
gtltcube_mp(bins=4, SCFile=SC, EVFile=sourcename+'_filtered_gti.fits', OutFile=sourcename+'_expCube.fits',SaveTemp=False,zmax=zmax) 
#expCube['zmax'] = zmax
#expCube['evfile'] = sourcename+'_filtered_gti.fits'
#expCube['scfile'] = SC
#expCube['outfile'] = sourcename+'_expCube.fits'
#expCube['dcostheta'] = 0.025
#expCube['binsz'] = 1
#expCube.run()
           
###----------------------------------------------------------------------------------------------  
print 'Step5:***Running gtexpcube2***'

gtexpcube2['infile'] = sourcename+'_expCube.fits'
gtexpcube2['outfile'] = sourcename+'_BinnedExpMap.fits'
gtexpcube2['cmap'] = 'none'
gtexpcube2['irfs'] = irfs
gtexpcube2['nxpix'] = 400
gtexpcube2['nypix'] = 400
gtexpcube2['binsz'] = 0.1           
gtexpcube2['coordsys'] = 'CEL'
gtexpcube2['xref'] = ra
gtexpcube2['yref'] = dec
gtexpcube2['axisrot'] = 0.0
gtexpcube2['proj'] = 'AIT'
gtexpcube2['emin'] = emin
gtexpcube2['emax'] = emax
gtexpcube2['enumbins'] = 12            
gtexpcube2.run()
            
##----------------------------------------------------------------------------------------------  
#print 'Step6:***Running gtsrcmaps***'
#model = "./model.xml"  #xml file 
#srcMaps['expcube'] = sourcename+'_expCube.fits'
#srcMaps['cmap'] = sourcename+'_CCUBE.fits'
#srcMaps['bexpmap'] = sourcename+'_BinnedExpMap.fits' 
#srcMaps['srcmdl'] = model
#srcMaps['outfile'] = sourcename+'_srcMaps.fits' 
#srcMaps['irfs'] = 'CALDB'
#srcMaps.run()
#---------------------------------------------------------------------------------------------- 

print "Finish!!!"
