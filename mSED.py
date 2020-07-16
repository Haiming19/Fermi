#!/usr/bin/env python
from UpperLimits import UpperLimits
#Import all of the needed modules
import pyLikelihood as pyLike
#import UnbinnedAnalysis as ua
import BinnedAnalysis as ba
import numpy as np
from os import path

from gt_apps import filter, maketime, expCube, expMap, diffResps, evtbin, srcMaps, gtexpcube2
from  gtltcube_mp import *

f1=open("SED_3c138_3.txt","w+") 
emin_select = '100'
emax_select = '300000'
binnum=4
bin_i = np.logspace(np.log10(float(emin_select)), np.log10(float(emax_select)), int(binnum) + 1)
emin_bin = bin_i[:-1]
emax_bin = bin_i[1:]
sourcename="3c138"


for i in range(0,4):
    
    print "The %dth bin" %i
    time_start = 239557417. #time(s)
    time_end = 584949858.  #8 years
    ra = 80.2912    # ra (degree)
    dec = 16.6395    # dec (degree)
    evclass = 128
    evtype = 3         #=front(1)+black(2)
    rad = 10.          #ROI (degree) 30
    emin = emin_bin[i]        #Emin unit in MeV
    emax = emax_bin[i]     #Emax unit in MeV
    zmax = 105          #maximum zenith angle (degree)
    irfs = "P8R3_SOURCE_V2"
    SC = "SC.fits"     #spacecraft file 
    infile = "@events.list"    #event file 
    
    #----------------------------------------------------------------------------------------------    
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
    #filter['outfile'] = sourcename+'_'+str(i)+'_filtered_.fits'
    #filter.run()
    ####----------------------------------------------------------------------------------------------
    #print 'Step2:***Running gtmktime***'
    
    #maketime['scfile'] = SC
    #maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
    #maketime['roicut'] = 'no'
    #maketime['evfile'] = sourcename+'_'+str(i)+'_filtered_.fits'
    #maketime['outfile'] = sourcename+'_'+str(i)+'_filtered_gti.fits'
    #maketime.run()
    ######-----------------------Create a source model XML file--------------------------------------------------
    ###mymodel = srcList('gll_psc_v16.fit',sourcename+'_filtered_gti.fits','model_3FGL.xml')  
    ###mymodel.makeModel('gll_iem_v06.fits', 'gll_iem_v06', 'iso_P8R3_SOURCE_V2.txt', 'iso_P8R3_SOURCE_V2')
    ######-----------------------Create a source model XML file 4FGL--------------------------------------------------
    ###mymodel = srcList('gll_psc_v17.fit',sourcename+'_filtered_gti.fits','model.xml')  
    ###mymodel.makeModel('gll_iem_v07.fits', 'gll_iem_v07', 'iso_P8R3_SOURCE_V2_v1.txt', 'iso_P8R3_SOURCE_V2_v1')
    ####-----------------------Create a source model XML file 3FHL--------------------------------------------------
    ###mymodel = srcList('gll_psch_v13.fit',sourcename+'_filtered_gti.fits','model_3FHL.xml')  
    ###mymodel.makeModel('gll_iem_v07.fits', 'gll_iem_v07', 'iso_P8R3_SOURCE_V2_v1.txt', 'iso_P8R3_SOURCE_V2_v1')
    ####-------------------------------------------------------------------------                                  
    ####----------------------------------------------------------------------------------------------
    #print 'Step3:***Running gtbin***'
    
    #evtbin['algorithm'] = 'CCUBE'
    #evtbin['evfile'] = sourcename+'_'+str(i)+'_filtered_gti.fits'
    #evtbin['outfile'] = sourcename+'_'+str(i)+'_CCUBE.fits'
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
    #evtbin['enumbins'] = 12            
    #evtbin.run()
    #####---------------------------------------------------------------------------------------------- 
    ###"python gtltcube_mp.py 2 SC "+name+"_binned_gti.fits "+name+"_ltcube.fits --zmax "+zmax
    ##gtltcube_mp(bins=2, SCFile=SC, EVFile=sourcename+'_'+str(i)+'_filtered_gti.fits', OutFile=sourcename+'_'+str(i)+'_expCube.fits', SaveTemp, zmax)
    #print 'Step4:***Running gtltcube***'
    #gtltcube_mp(bins=4, SCFile=SC, EVFile=sourcename+'_'+str(i)+'_filtered_gti.fits', OutFile=sourcename+'_'+str(i)+'_expCube.fits',SaveTemp=False,zmax=zmax)   
    ###expCube['zmax']=zmax
    ###expCube['evfile'] = sourcename+'_'+str(i)+'_filtered_gti.fits'
    ###expCube['scfile'] = SC
    ###expCube['outfile'] = sourcename+'_'+str(i)+'_expCube.fits'
    ###expCube['dcostheta'] = 0.025
    ###expCube['binsz'] = 1
    ###expCube.run()
               
    ######----------------------------------------------------------------------------------------------  
    #print 'Step5:***Running gtexpcube2***'
    
    #gtexpcube2['infile'] = sourcename+'_'+str(i)+'_expCube.fits'
    #gtexpcube2['outfile'] = sourcename+'_'+str(i)+'_BinnedExpMap.fits'
    #gtexpcube2['cmap'] = 'none'
    #gtexpcube2['irfs'] = irfs
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
    #gtexpcube2['enumbins'] = 12            
    #gtexpcube2.run()
                
    #####----------------------------------------------------------------------------------------------  
    print 'Step6:***Running gtsrcmaps***'
    model = "./model.xml"  #xml file 
    srcMaps['expcube'] = sourcename+'_'+str(i)+'_expCube.fits'
    srcMaps['cmap'] = sourcename+'_'+str(i)+'_CCUBE.fits'
    srcMaps['bexpmap'] = sourcename+'_'+str(i)+'_BinnedExpMap.fits' 
    srcMaps['srcmdl'] = model
    srcMaps['outfile'] = sourcename+'_'+str(i)+'_srcMaps.fits' 
    srcMaps['irfs'] = 'CALDB'
    srcMaps.run()
    
    ##----------------------------------------------------------------------------------------------              
    print '***Running likelihood analysis***'
    #---------------------------------------------------------------------
    #model = "./model1.xml"
    #model = "./modelfitted2.xml"
    obs = ba.BinnedObs(srcMaps=sourcename+'_'+str(i)+'_srcMaps.fits',
                    expCube=sourcename+'_'+str(i)+'_expCube.fits',
                    binnedExpMap=sourcename+'_'+str(i)+'_BinnedExpMap.fits',
                    irfs='CALDB')
    
    
    like = ba.BinnedAnalysis(obs,model,optimizer='MINUIT')
    like.tol = 1e-8
    likeObj = pyLike.Minuit(like.logLike)
    print "#------------------------------------------------------#"
    loglikelihood= like.fit(verbosity=0,covar=True,optObject=likeObj)
    print loglikelihood
    like.logLike.writeXml('modelfitted'+'_'+str(i)+'.xml')
    
    sourcename1="4FGL J0521.2+1637"
    
    #-------------------------------------------------
    ts=like.Ts(sourcename1)
    print ts
    if ts<9:
        model1="./model2.xml" 
        like1 = ba.BinnedAnalysis(obs,model1,optimizer='MINUIT')
        like1.tol = 1e-8
        likeObj1 = pyLike.Minuit(like1.logLike)
        like1.fit(verbosity=0,covar=True,optObject=likeObj1)
        ul=UpperLimits(like1)
        upperlimits=ul[sourcename1].compute(emin=emin,emax=emax)#delta=4.8 means 3sigma
        flux = upperlimits[0]
        fluxerr = 0.0
        ecent=10**((np.log10(emax)+np.log10(emin))*0.5)
        Corr=1./(emax-emin)
        vfv=1.6e-6*flux*Corr*ecent**2
        vfverr=0.0
        Gamma=2.0
        GammaErr=0.0
    else:
        flux=like.flux(sourcename1,emin=emin ,emax=emax)
        fluxerr=like.fluxError(sourcename1,emin=emin ,emax=emax)
        Corr=1./(emax-emin) 
        ecent=10**((np.log10(emax)+np.log10(emin))*0.5)
        vfv=1.6e-6*flux*Corr*ecent**2
        vfverr=1.6e-6*fluxerr*Corr*ecent**2
        #Gamma = like.model[sourcename1].funcs['Spectrum'].getParam('Index').value()
        #GammaErr = like.model[sourcename1].funcs['Spectrum'].getParam('Index').error()
    print >>f1,ecent,ecent-emin,emax-ecent,vfv,vfverr,ts,flux,fluxerr
    print ecent,ecent-emin,emax-ecent,vfv,vfverr,ts,flux,fluxerr
    #-----------------------------------------
    
    #------------------------------
f1.close()
print '***Finished***'
