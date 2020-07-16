#!/usr/bin/env python

#Import all of the needed modules

import numpy as np
from os import path
from gt_apps import filter, maketime, expCube, expMap, diffResps, evtbin, srcMaps, gtexpcube2

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
zmax = 105          #maximum zenith angle (degree)
irfs = "P8R3_SOURCE_V2"
SC = "SC.fits"     #spacecraft file 
infile = "@events.list"    #event file 


##----------------------------------------------------------------------------------------------
print 'Step3:***Running gtbin***'

evtbin['algorithm'] = 'CMAP'
evtbin['evfile'] = sourcename+'_filtered_gti.fits'
evtbin['outfile'] = sourcename+'_CMAP.fits'
evtbin['scfile'] = 'NONE'
evtbin['nxpix'] = 40
evtbin['nypix'] = 40
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

print "Finish!!!"
