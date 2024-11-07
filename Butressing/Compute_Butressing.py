#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 15:16:51 2018

@author: cmosbeux
"""

#%%
#Import some usfull library
#--------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util import numpy_support as VN
#from  matplotlib_toolkits.matplotlibot3d import Axes3D
from os import listdir
from os.path import isfile, isdir, join
import matplotlib
import sys
import re
import os
sys.path.append('/Users/cmosbeux/Documents/Developper/MyModules')
sys.path.append('/Users/cmosbeux/Documents/Developper/Python/My_python_scripts')
sys.path.append('/Users/cmosbeux/Documents/Developper/My_python_scripts/Graph_Data')
from mymodule import MidpointNormalize, colorbar
from lecture_file import lecture_contour

from drawcoast import *
import cartopy.crs as ccrs
from cartopy._crs import CRS, Geocentric, Geodetic, Globe, PROJ4_RELEASE
from SouthProj import *
import cartopy.io.shapereader as shpreader
import matplotlib.colors as mcolors
#from matplotlib.colors import colorConverter
import matplotlib.cm as cmp
from matplotlib.colors import LightSource, LogNorm, Normalize
from colorbar_terrain import rvb, rvb2
import math

import warnings
warnings.filterwarnings("ignore")

#%%------------------------------------------------------------------------
'''Read variable(s) in VTU files'''
#--------------------------------------------------------------------------

print '#########################################'
print '#	        VTU treatment 	        #'
print '#########################################'

input_directory='../output@comet/mesh2D_HR/'
onlyfiles = [ f for f  in listdir(input_directory) if isfile(join(input_directory+'/'+f)) ]

#What variables do you want to plot/use?
coords = []
list_of_var = ['sigma 1', 'sigma 2', 'sigma 3', 'sigma 4',  'ssavelocity_0', 'h', 'eta', 'vx', 'vy', 'gmask', 'zb' ] 
v = {}
for i in list_of_var:
    v[i]=[]


for f in onlyfiles:  
    fname = f
    fext= f.split('.')[1]
    if ('BMB_Initial_reg000_HR_' not in f) or (fext!='vtu'):
        continue
    fname = re.sub('\.vtu$', '', fname)
    print 'vtu inspected : ', f
    #Read the source file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(input_directory+'/'+f)
    reader.Update()
    output=reader.GetOutput()
    PointData=output.GetPointData()
    coords.extend(VN.vtk_to_numpy(output.GetPoints().GetData()))
    for i in list_of_var:
        array=VN.vtk_to_numpy(PointData.GetArray(i))
        v[i].extend(array)


#change lists to array
def arr(i):
    i=np.asarray(i)
    return i

#coords
coord_x=arr(coords).T[0]
coord_y=arr(coords).T[1]
           
#flow direction
ssax=arr(v['ssavelocity_0']).T[0]
ssay=arr(v['ssavelocity_0']).T[1]
v_norm=(ssax**2+ssay**2)**0.5
nx=ssax/v_norm
ny=ssay/v_norm
n=np.array([nx,ny])
ori1 = [math.acos(i) for i in (nx)]  
ori2 = [math.acos(i) for i in (ny)]  

#vertically integrated pressure from ocean
g=9.81
ri=917.
rw=1028.
h=arr(v['h'])
N0= 0.5*ri*(1-ri/rw)*g*h * 1e-6

p=0.5*ri*g*h * 1e-6
sxx=arr(v['sigma 1']) #2*v['viscosity']*v['strainrate 1']
syy=arr(v['sigma 2']) #2*v['viscosity']*v['strainrate 2']+0.5*ri*g*v['h']*1e-6
sxy=arr(v['sigma 4']) #+abs(0.5*ri*g*v['h']*1e-6) #2*v['viscosity']*v['strainrate 4']
szz=arr(v['sigma 3']) #2*v['viscosity']*v['strainrate 3']+0.5*ri*g*v['h']*1e-6

#n.Sigma.n 
sig = (nx*(sxx*nx+sxy*ny)+ny*(sxy*nx+syy*ny))
sigma1 = (sxx+syy)/2 + ((sxx-syy)**2/4+ sxy**2)**0.5
sigma2 = (sxx+syy)/2 - ((sxx-syy)**2/4+ sxy**2)**0.5

eta=arr(v['eta'])
sigma=sigma1
#s=sig/(N0)
#xxx=np.arange(0,len(N0),1)
#plt.scatter(xxx,-N0, s=3, marker='.', linewidth=0)
#plt.scatter(xxx,-sig, s=3, marker='.', linewidth=0)
#plt.scatter(xxx,s, s=3, marker='.', linewidth=0)
#plt.ylim(0,3)

K=1-(sigma/N0)

#%%
plt.figure()
cmap=cmp.get_cmap('jet', 20)
#plt.scatter(coord_x, coord_y, c=K, cmap=cmap, norm=Normalize(vmin=-7.26, vmax=-7.24), s=6, marker='.', linewidth=0)
plt.scatter(coord_x, coord_y, c=K, cmap=cmap, norm=Normalize(vmin=0.0, vmax=2), s=6, marker='.', linewidth=0)

#plt.scatter(coord_x,coord_y, c=var_plt, cmap=cmap, norm=Normalize(vmin=-0.005, vmax=0.000), s=6, marker='.', linewidth=0)
#plt.scatter(coord_x,coord_y, c=var_plt, cmap=cmap, norm=LogNorm(vmin=0.01, vmax=1), s=15, marker='.', linewidth=0)
plt.colorbar()

#%%
#plt.scatter(np.arange(0,len(sig),1),abs(sig/N0),marker='.', s=1) 
#plt.ylim(0,20)
#plt.grid()

