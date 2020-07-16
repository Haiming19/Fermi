import xml.dom.minidom as minidom

def parameter_element(free, name, maximum, minimum, scale, value):
    """Create an XML document parameter description element"""
    impl = minidom.getDOMImplementation()
    xmldoc_out = impl.createDocument(None,None,None)
    parameter = xmldoc_out.createElement('parameter')
    parameter.setAttribute('free', str(free))
    parameter.setAttribute('name', str(name))
    parameter.setAttribute('max', str(maximum))
    parameter.setAttribute('min', str(minimum))
    parameter.setAttribute('scale', str(scale))
    parameter.setAttribute('value', str(value))
    return parameter


def sp2pl(infile,outfile,name,PrefactorScale=1.e-13,index=2.2,scaleValue=1000.0,freeIndex=True):
    xmldoc = minidom.parse(infile)
    if 0: # test
        xmlStr=xmldoc.toprettyxml(' ').splitlines(True)
        outStr=filter(lambda xmlStr: len(xmlStr) and not xmlStr.isspace(),xmlStr)
        outxml=open('test0.xml','w')
        outxml.write(''.join(outStr))
        outxml.close()

    new_sp= xmldoc.createElement('spectrum')
    new_sp.setAttribute('apply_edisp','false')
    new_sp.setAttribute('type',"PowerLaw")
    new_sp.appendChild(parameter_element("1", "Prefactor", "1e+04", "1e-04", "%s"%str(PrefactorScale),"1.0"))
    new_sp.appendChild(parameter_element("%s"%int(freeIndex), "Index", "3.5", "0.0", "-1.0", "%s"%str(index)))
    new_sp.appendChild(parameter_element("0", "Scale", "500000.0", "30.0", "1.0", "%s"%str(scaleValue)))

    name1 = name
    slist = xmldoc.getElementsByTagName('source')
    for s in slist:
        srcName = s.getAttribute('name')
        if srcName != name1: continue
        old_sp  = s.getElementsByTagName('spectrum')[0]
        s.replaceChild(new_sp,old_sp)
        break    
    
    xmlStr=xmldoc.toprettyxml(' ').splitlines(True)
    outStr=filter(lambda xmlStr: len(xmlStr) and not xmlStr.isspace(),xmlStr)
    outxml=open(outfile,'w')
    outxml.write(''.join(outStr))
    outxml.close()


import numpy as np
import pprint
def printResult(like,saveresultto,extraInfo2Print={},summedLike=None):
    e1 = like.energies.min()
    e2 = like.energies.max()
    if summedLike != None: like = summedLike
    if like.covar_is_current is False:
        print("covar is not calculated yet")
        return 0

    freeParNames = []
    index = 0
    dicttot = {"logLikelihood": like.logLike.value()}
    for source in like.sourceNames():
        if like.model[source].src.fixedSpectrum(): continue
        TsValue = like.Ts(source)
        dict_ = {}
        #print source,like.Ts(source)

        dict_['TS value'] = np.str(TsValue)
        flux_ = like.flux(source,emin=e1,emax=e2,energyFlux=False)
        fluxErr_ = like.fluxError(source,emin=e1,emax=e2,energyFlux=False)
        dict_['Flux'] = "%.5e +/- %.5e" %(flux_,fluxErr_)

        for pname in like.model[source].funcs['Spectrum'].paramNames:
            error_ = like.model[source].funcs['Spectrum'].params[pname].error()
            value_ = like.model[source].funcs['Spectrum'].params[pname].value()
            dict_[pname] = "%f +/- %f" %(value_,error_)
        
            bounds_ = like.model[source].funcs['Spectrum'].params[pname].getBounds()
            if ((value_-bounds_[0]) < 1.e-2) or ((bounds_[1]-value_) < 1.e-2):
                print("Warning:  " +pname+ "  of  " +source+ "  is close to initial value!")
                print(value_,bounds_[0],bounds_[1])
        dict_['SpectrumType'] = like.model[source].funcs['Spectrum'].func.genericName()
        dicttot[source] = dict_

        for p in like.freePars(source):
            pname = p.getName()
            isNorm = 0
            if pname == like.model[source].src.spectrum().normPar().getName(): isNorm = 1
            freeParNames.append((index,source,pname,isNorm,TsValue))
            index += 1

    cov = np.array(like.covariance)
    cov_diag = np.diag(cov)
    ab = np.repeat(cov_diag,len(cov_diag)).reshape((len(cov_diag),len(cov_diag))).T
    cor = np.array(like.covariance)/np.sqrt(ab*ab.T)
    cor = (cor - np.eye(len(cov_diag)))

    badParameters = []
    for i,source,pname,isnorm,ts in freeParNames:
        #if source == srcName: continue
        hindex = i
        Vindexs = np.where(np.abs(cor[i][:i+1]) > 0.9)[0]
        if Vindexs.size == 0: continue
        for vindex in Vindexs:
            indexv_,sourcev_,pnamev_,isNormv_,TsValuev_ = freeParNames[vindex]
            indexh_,sourceh_,pnameh_,isNormh_,TsValueh_ = freeParNames[hindex]
            #fixPar(hindex,vindex):
            if (TsValuev_ < TsValueh_):
                #like2.freeze(like2.par_index(sourcev_,pnamev_))
                badParameters.append((sourcev_,pnamev_))
            if (TsValueh_ < TsValuev_):
                #like2.freeze(like2.par_index(sourceh_,pnameh_))
                badParameters.append((sourceh_,pnameh_))
            if (TsValueh_ == TsValuev_):
                #if isNormv_: like2.freeze(like2.par_index(sourceh_,pnameh_))
                #if isNormh_: like2.freeze(like2.par_index(sourcev_,pnamev_))
                if isNormv_: badParameters.append((sourceh_,pnameh_))
                if isNormh_: badParameters.append((sourcev_,pnamev_))

    f = open(saveresultto,'w')
    pprint.pprint(dicttot,stream=f)
    print>>f, "#------------------------------The elegant separation line----------------------"
    for index,source,pname,isNorm,TsValue in freeParNames:
        print>>f, "#%d , %s , %s , %.3f" %(index,source,pname,TsValue)
    print>>f, "#------------------------------The elegant separation line----------------------"
    for i in range(len(cor)):
        tmpstr = "#%3d" %(i)
        for k in range(len(cor)):
            tmpstr = tmpstr + " %.3f " %(cor[i][k])
        print>>f, tmpstr
    print>>f, "#------------------------------The elegant separation line----------------------"
    for key in extraInfo2Print.keys():
        print>>f, "#",key,extraInfo2Print[key]
    f.close()

    return badParameters




if __name__=='__main__':
    sp2pl('m_0.xml','test.xml','4FGL J1441.4-1934')

