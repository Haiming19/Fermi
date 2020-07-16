#!/usr/bin/env python

#Import all of the needed modules
import pyLikelihood as pyLike
#import UnbinnedAnalysis as ua
import BinnedAnalysis as ba
import numpy as np
from os import path
from gt_apps import filter, maketime, expCube, expMap, diffResps, evtbin, srcMaps, gtexpcube2
from UpperLimits import UpperLimits

class LightCurve():

    """A really simple light curve generator.  Not really robust.  If
    you want to really make a light curve, you shouldn't use this.
    This is just an example of what can be done.  The class has an
    init function and two other functions.  If you run the
    generate_files function then the runLikelihood function you can
    get some data you can plot."""

    def __init__(self,time_start = 239557417., time_end = 584949858., delta_t =86400*180):  #bin=1d. 7d  

        self.times = np.arange(time_start,time_end,delta_t)

        self.ra = 80.2912     #best fitted ra
        self.dec = 16.6395      #best fitted dec
        self.evclass = 128
        self.evtype = 3
        self.rad = 10.
        self.emin = 100.
        self.emax = 300000.
        self.zmax = 105
        self.irfs = "P8R3_SOURCE_V2"
        self.SC = "SC.fits"
        self.source = "4FGL J0521.2+1637"
        self.model = "./lc.xml"
        self.infile = "@events.list"
        self.model1 = "./lc2.xml"

        #Set up some arrays to store the results
        self.IntenergyFlux = np.zeros_like(self.times[:-1])
        self.IntenergyFluxErr = np.zeros_like(self.times[:-1])
        self.IntFlux = np.zeros_like(self.times[:-1])
        self.IntFluxErr = np.zeros_like(self.times[:-1])
        self.Gamma = np.zeros_like(self.times[:-1])
        self.GammaErr = np.zeros_like(self.times[:-1])
        self.TS = np.zeros_like(self.times[:-1])
        self.retCode = np.zeros_like(self.times[:-1])
        self.npred = np.zeros_like(self.times[:-1])

    def dolightcurve(self):

        """This function makes all the needed files like the livetime
        cube and the exposure map."""
        f1=open("result_ul.txt","w+")
        f=open("result.txt","w+") 
        #f3=open("result3.txt","w+")         
        #Loop through the start (times[:-1])/stop(times[1:]) times
        for i,time in enumerate(zip(self.times[:-1],self.times[1:])):     
            
            print "**Working on the %d th bin: range (%f,%f)" % (i,time[0],time[1])
#----------------------------------------------------------------------------------------------    
            print 'Step1:***Running gtselect***'
            filter['evclass'] = self.evclass
            filter['evtype'] = self.evtype
            filter['ra'] = self.ra
            filter['dec'] = self.dec
            filter['rad'] = self.rad
            filter['emin'] = self.emin
            filter['emax'] = self.emax
            filter['zmax'] = self.zmax
            filter['tmin'] = time[0]
            filter['tmax'] = time[1]
            filter['infile'] = self.infile
            filter['outfile'] = 'LC_filtered_binned_'+str(i)+'.fits'
            filter.run()
        
#----------------------------------------------------------------------------------------------
            print 'Step2:***Running gtmktime***'
            maketime['scfile'] = self.SC
            maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
            maketime['roicut'] = 'no'
            maketime['evfile'] = 'LC_filtered_binned_'+str(i)+'.fits'
            maketime['outfile'] = 'LC_filtered_gti_binned_'+str(i)+'.fits'
            maketime.run()
            
#----------------------------------------------------------------------------------------------
            print 'Step3:***Running gtbin***'
            evtbin['algorithm'] = 'CCUBE'
            evtbin['evfile'] = 'LC_filtered_gti_binned_'+str(i)+'.fits'
            evtbin['outfile'] = 'LC_CCUBE_'+str(i)+'.fits'
            evtbin['scfile'] = 'NONE'
            evtbin['nxpix'] = 140
            evtbin['nypix'] = 140
            evtbin['binsz'] = 0.1            
            evtbin['coordsys'] = 'CEL'
            evtbin['xref'] = self.ra
            evtbin['yref'] = self.dec
            evtbin['axisrot'] = 0.0
            evtbin['proj'] = 'TAN'
            evtbin['ebinalg'] = 'LOG'
            evtbin['emin'] = self.emin
            evtbin['emax'] = self.emax
            evtbin['enumbins'] = 20            
            evtbin.run()
                       
#----------------------------------------------------------------------------------------------            
            print 'Step4:***Running gtltcube***'
            #gtltcube_mp(bins=6, SCFile=self.SC, EVFile='LC_filtered_gti_binned_'+str(i)+'.fits', OutFile='LC_expCube_binned_'+str(i)+'.fits',SaveTemp=False,zmax=self.zmax) 

            expCube['zmax'] = self.zmax
            expCube['evfile'] = 'LC_filtered_gti_binned_'+str(i)+'.fits'
            expCube['scfile'] = self.SC
            expCube['outfile'] = 'LC_expCube_binned_'+str(i)+'.fits'
            expCube['dcostheta'] = 0.025
            expCube['binsz'] = 1
            expCube.run()
           
#----------------------------------------------------------------------------------------------  
            print 'Step5:***Running gtexpcube2***'
            gtexpcube2['infile'] = 'LC_expCube_binned_'+str(i)+'.fits'
            gtexpcube2['outfile'] = 'LC_BinnedExpMap_'+str(i)+'.fits'
            gtexpcube2['cmap'] = 'none'
            gtexpcube2['irfs'] = self.irfs
            gtexpcube2['nxpix'] = 400
            gtexpcube2['nypix'] = 400
            gtexpcube2['binsz'] = 0.1            
            gtexpcube2['coordsys'] = 'CEL'
            gtexpcube2['xref'] = self.ra
            gtexpcube2['yref'] = self.dec
            gtexpcube2['axisrot'] = 0.0
            gtexpcube2['proj'] = 'TAN'
            gtexpcube2['emin'] = self.emin
            gtexpcube2['emax'] = self.emax
            gtexpcube2['enumbins'] = 20            
            gtexpcube2.run()
            
#----------------------------------------------------------------------------------------------  
            print 'Step6:***Running gtsrcmaps***'
            srcMaps['expcube'] = 'LC_expCube_binned_'+str(i)+'.fits'
            srcMaps['cmap'] = 'LC_CCUBE_'+str(i)+'.fits'
            srcMaps['bexpmap'] = 'LC_BinnedExpMap_'+str(i)+'.fits' 
            srcMaps['srcmdl'] = self.model
            srcMaps['outfile'] = 'LC_srcMaps_'+str(i)+'.fits' 
            srcMaps['irfs'] = 'CALDB'
            srcMaps.run()
            
#----------------------------------------------------------------------------------------------              
            print '***Running likelihood analysis***'
            print i
#---------------------------------------------------------------------
            obs = ba.BinnedObs(srcMaps='LC_srcMaps_'+str(i)+'.fits',
                            expCube='LC_expCube_binned_'+str(i)+'.fits',
                            binnedExpMap='LC_BinnedExpMap_'+str(i)+'.fits',
                            irfs='CALDB')
            like = ba.BinnedAnalysis(obs,self.model,optimizer='MINUIT')
            #like1 = ba.BinnedAnalysis(obs,self.model1,optimizer='MINUIT')
            #---------------------------------------------------------------------
            #obs = ua.UnbinnedObs('LC_filtered_gti_bin'+str(i)+'.fits',
                                 #self.SC,
                                 #expMap='LC_expMap_bin'+str(i)+'.fits',
                                 #expCube='LC_expCube_bin'+str(i)+'.fits',
                                 #irfs=self.irfs)
            #like = ua.UnbinnedAnalysis(obs,self.model,optimizer='NewMinuit')
            #---------------------------------------------------------------------
            
            like.tol = 1e-8
            likeObj = pyLike.Minuit(like.logLike)
            loglikelihood=like.fit(verbosity=0,covar=True,optObject=likeObj)
     #       like.fit(verbosity=0,covar=True)  
     
            like.logLike.writeXml(str(i)+'.xml')  
            

            self.IntenergyFlux[i] = like.energyFlux(self.source,emin=self.emin,emax=self.emax)
            self.IntenergyFluxErr[i] = like.energyFluxError(self.source,emin=self.emin,emax=self.emax)
            self.IntFlux[i] = like.flux(self.source,emin=self.emin,emax=self.emax)
            self.Gamma[i] = like.model[self.source].funcs['Spectrum'].getParam('Index').value()
            self.IntFluxErr[i] = like.fluxError(self.source,emin=self.emin,emax=self.emax)
            self.GammaErr[i] = like.model[self.source].funcs['Spectrum'].getParam('Index').error()
            self.TS[i] = like.Ts(self.source)
            self.retCode[i] = likeObj.getRetCode()
            self.npred[i] = like.NpredValue(self.source)
            ###zhm
            print i,time[0],time[1],self.IntFlux[i],self.IntFluxErr[i],self.Gamma[i],self.GammaErr[i],self.TS[i],self.IntenergyFlux[i],self.IntenergyFluxErr[i],loglikelihood
            print >>f,i,time[0],time[1],self.IntFlux[i],self.IntFluxErr[i],self.Gamma[i],self.GammaErr[i],self.TS[i],self.IntenergyFlux[i],self.IntenergyFluxErr[i],loglikelihood
            if (self.TS[i] <9):
                like1 = ba.BinnedAnalysis(obs,self.model1,optimizer='MINUIT')
                ul=UpperLimits(like1)
                UL=ul[self.source].compute(emin=self.emin,emax=self.emax)
                print >>f1,i,time[0],time[1],UL[0]
                print i,time[0],time[1],UL[0]
            
           # source="4FGL J0515.8+1527"
            #IntenergyFlux = like.energyFlux(source,emin=self.emin,emax=self.emax)
           # IntenergyFluxErr = like.energyFluxError(source,emin=self.emin,emax=self.emax)
            #IntFlux = like.flux(source,emin=self.emin,emax=self.emax)
            #Gamma = like.model[source].funcs['Spectrum'].getParam('Index').value()
            #IntFluxErr = like.fluxError(source,emin=self.emin,emax=self.emax)
            #GammaErr = like.model[source].funcs['Spectrum'].getParam('Index').error()
            #TS = like.Ts(source)
            #print >>f3,i,time[0],time[1],IntFlux,IntFluxErr,Gamma,GammaErr,TS
            #Clean up at the end of this iteration by
            #deleting some objects and incrementing the
            #bin start and stop times as well as the bin number
            
            ###zhm
            del likeObj
            del like
            del obs  
        #f3.close()         
        f1.close()
        f.close()

#if __name__ == '__main__':
    #lc = LightCurve()
    #lc.dolightcurve()
