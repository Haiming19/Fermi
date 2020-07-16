#!/usr/bin/env python

#Import all of the needed modules

import numpy as np
from os import path
from scipy import stats


#data=np.loadtxt("result160.txt",unpack=True)
#print data[10]
#data2=np.loadtxt("result162.txt",unpack=True)
#print data2[10]

#data=np.loadtxt("result.txt",unpack=True)
#print data[10]
#data2=np.loadtxt("result2.txt",unpack=True)
#print data2[10]


TS_var=45.82
#for i in range(len(data[10])):
    #TS_var=TS_var+2*(-data[10][i]+data2[10][i])
    
#print TS_var

#ndof=len(data[10])
#print ndof
ndof=22
#TS_var=42.75

prop= stats.chi2.sf(TS_var,ndof-1)
print prop
sig = stats.norm.isf(prop)

print'%s sigma'%sig
