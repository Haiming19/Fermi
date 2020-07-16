import numpy as np
import scipy.ndimage as ndimage
import sys
import astropy.io.fits as pyfits
import matplotlib
import re
from math import pi, cos, sin

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.wcs as wcs
def get_circle(ra,dec,size,ts_bin,coor_obj,circle_pixel=0.05):
    ### plot the circle 
    c_ra_pixel,c_dec_pixel = coor_obj.wcs_world2pix(ra,dec,1)
    theta = np.arange(0, 2*np.pi, circle_pixel)  
    r_x = c_ra_pixel-0.5 + size * (1/ts_bin) * np.cos(theta)
    r_y = c_dec_pixel-0.5 + size * (1/ts_bin) * np.sin(theta)
    return (r_x,r_y)

 
name = sys.argv[1]
residHDU = pyfits.open(name+".fits")
wcsobj = wcs.WCS(name+".fits")
### Select J2000 position. 
### Decide coordinate range.
ts_bin = np.abs(residHDU[0].header['CDELT2'])
ra = residHDU[0].header['CRVAL1']
dec = residHDU[0].header['CRVAL2']
left = 0
right = residHDU[0].header['NAXIS1']
bottom = residHDU[0].header['NAXIS2']
top = 0
deltar = residHDU[0].header['CDELT2']

### Plot map.
data_set = residHDU[0].data
plt.figure(figsize=(24,8))
plt.subplot(1,2,1)
#plt.xlabel(residHDU[0].header['CTYPE1'])
#plt.ylabel(residHDU[0].header['CTYPE2'])
plt.xlabel("R.A. [degrees]")
plt.ylabel("Decl. [degrees]")
data_set = ndimage.gaussian_filter(data_set,1)#0.5
#plt.imshow(data_set,aspect='auto',cmap=plt.cm.jet,extent=[left,right,bottom,top],vmin=1.0)
#plt.imshow(data_set,aspect='auto',cmap=plt.cm.gnuplot2,extent=[left,right,bottom,top],vmin=4.0,vmax=25)
#plt.imshow(data_set,aspect='auto',cmap=plt.cm.gnuplot2,extent=[left,right,bottom,top],vmin=-0.1,vmax=0.3)
plt.imshow(data_set,aspect='auto',cmap=plt.cm.gnuplot2,extent=[left,right,bottom,top],vmin=4.0,vmax=20)
cbar=plt.colorbar()
cbar.set_label('TS value')
#cbar.set_label('Counts/cm2/s/MeV/sr')
#cbar.set_label('Counts/pixel')

### Plot the 4FGL source point source.
f = open(name+'.reg')
source_FGL = f.readlines()
f.close()
rs = r'J2000;point\(|,|\)\s\#|text={|}\n'
for i in np.arange(np.size(source_FGL))[3:]:
    position_FGL_ra = re.split(rs,source_FGL[i])[1]
    position_FGL_dec = re.split(rs,source_FGL[i])[2]
    position_FGL_ra_pixel,position_FGL_dec_pixel= wcsobj.wcs_world2pix(float(position_FGL_ra),float(position_FGL_dec),1)  
    position_FGL_name = re.split(rs,source_FGL[i])[-2] 
    if left<position_FGL_ra_pixel<right and bottom>position_FGL_dec_pixel>top:
        plt.plot(position_FGL_ra_pixel,position_FGL_dec_pixel,color='c',marker='+',markersize=4,ls='none')
        plt.text(position_FGL_ra_pixel+0.2,position_FGL_dec_pixel-0.5,position_FGL_name,color='c',fontsize=8)

#80.314,16.631
#plot 68 95 -----------------------------------------------------------------
lon_center=80.314 #x-position of the center
lat_center=16.631 #y-position of the center
aaxis=0.0464/0.05       #radius on the x-axis,0.9/0.025=36;
baxis=0.0387/0.05       #radius on the y-axis,0.5/0.025=20;
rad=pi/180
t_rot=-73.98*rad      #rotation angle in rad;
c_ra_pixel,c_dec_pixel = wcsobj.wcs_world2pix(lon_center,lat_center,1)

#theta2 = np.linspace(0, 2*np.pi, 180)
#r_x2 = c_ra_pixel2+aaxis*np.cos(theta2)
#r_y2 = c_dec_pixel2+baxis*np.sin(theta2)
##plt.grid(color='lightgray',linestyle='--')
#plt.plot(r_x2,r_y2,color='k',marker='.',markersize=3,ls='none')
#plt.annotate('RSGC1',xy=(c_ra_pixel2-2, c_dec_pixel2+3))
#
theta = np.linspace(0, 2*np.pi, 50)
Ell = np.array([aaxis*np.cos(theta),baxis*np.sin(theta)])  
R_rot = np.array([[cos(t_rot) , -sin(t_rot)],[sin(t_rot) , cos(t_rot)]]) #2-D rotation matrix;
Ell_rot = np.zeros((2,Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
#plt.plot( c_ra_pixel2+Ell[0,:],c_dec_pixel2+Ell[1,:] )     #initial ellipse
plt.plot( c_ra_pixel+Ell_rot[0,:] , c_dec_pixel+Ell_rot[1,:],'g',linewidth=0.7)    #rotated ellipse
#plt.annotate('HESS J',xy=(c_ra_pixel+4, c_dec_pixel),color='g', size=20)
#--------------------------------------------------------------------------

aaxis=0.0752/0.05       #radius on the x-axis,0.9/0.025=36;
baxis=0.0627/0.05       #radius on the y-axis,0.5/0.025=20;
rad=pi/180
t_rot=-73.98*rad      #rotation angle in rad;
c_ra_pixel,c_dec_pixel = wcsobj.wcs_world2pix(lon_center,lat_center,1)

#theta2 = np.linspace(0, 2*np.pi, 180)
#r_x2 = c_ra_pixel2+aaxis*np.cos(theta2)
#r_y2 = c_dec_pixel2+baxis*np.sin(theta2)
##plt.grid(color='lightgray',linestyle='--')
#plt.plot(r_x2,r_y2,color='k',marker='.',markersize=3,ls='none')
#plt.annotate('RSGC1',xy=(c_ra_pixel2-2, c_dec_pixel2+3))
#
theta = np.linspace(0, 2*np.pi, 50)
Ell = np.array([aaxis*np.cos(theta),baxis*np.sin(theta)])  
R_rot = np.array([[cos(t_rot) , -sin(t_rot)],[sin(t_rot) , cos(t_rot)]]) #2-D rotation matrix;
Ell_rot = np.zeros((2,Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
#plt.plot( c_ra_pixel2+Ell[0,:],c_dec_pixel2+Ell[1,:] )     #initial ellipse
plt.plot( c_ra_pixel+Ell_rot[0,:] , c_dec_pixel+Ell_rot[1,:],'g',linewidth=0.7)    #rotated ellipse
#plt.annotate('HESS J',xy=(c_ra_pixel+4, c_dec_pixel),color='g', size=20)
#-------------------------------------------------------------------------
position_name="3C 138"
lon_center1=80.2912 #x-position of the center
lat_center1=16.6394 #y-position of the center
c_ra_pixel1,c_dec_pixel1 = wcsobj.wcs_world2pix(lon_center1,lat_center1,1)
plt.plot(c_ra_pixel1,c_dec_pixel1,color='r',marker='o',markersize=4,ls='none')
plt.text(c_ra_pixel1+0.4,c_dec_pixel1+0.4,position_name,color='r',fontsize=8)

##-------------------------------------------------------------------------
plt.xlim(left,right)
plt.ylim(top,bottom)

### Set the left, right, bottom, top coordinates.

ax = plt.gca()
xx = np.linspace(left+1,right-1,5)
xxnamebottom = map(lambda x:str("%.1f"%(x)),wcsobj.wcs_pix2world(xx,bottom,1)[0])
yy = np.linspace(top+1,bottom-1,5)
yynameleft = map(lambda x:str("%.1f"%(x)),wcsobj.wcs_pix2world(left,yy,1)[1])

xx_float = map(lambda x:float("%.4f"%(x)),wcsobj.wcs_pix2world(xx,bottom,1)[0])
yy_float = map(lambda x:float("%.4f"%(x)),wcsobj.wcs_pix2world(left,yy,1)[1])
#for xx_each in xx_float:
    #dec_sample = np.linspace(yy_float[0]-ts_bin*5,yy_float[-1]+ts_bin*5,30)
    #each_pixel = wcsobj.wcs_world2pix(xx_each,dec_sample,1)
    #plt.plot(each_pixel[0],each_pixel[1],ls='-',color='darkblue',linewidth=0.5)
#for yy_each in yy_float:
    #ra_sample = np.linspace(xx_float[0]+ts_bin*5,xx_float[-1]-ts_bin*5,30)
    #each_pixel = wcsobj.wcs_world2pix(ra_sample,yy_each,1)
    #plt.plot(each_pixel[0],each_pixel[1],ls='-',color='darkblue',linewidth=0.5)

plt.setp(ax,xticks=xx,xticklabels=xxnamebottom,yticks=yy,yticklabels=yynameleft)
plt.savefig(name+'.eps',bbox_inches='tight')
plt.close()
print(name,'ok')

