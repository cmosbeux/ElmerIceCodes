#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 11:37:31 2020

@author: cmosbeux

This scrit creates a box from a catchement (with an surrounding extra width) or
with 4 given corners locations

You can check your catchment box in "catchment_box.png"
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import sys


#%%----------------------------------------------------------------------------
#Do you want to use a catchement?
catchement = False

#If not, precise your box dimensions (x0,x1) and (y0,y1)
x = [0,100e3]
y = [0,50e3]


#%%----------------------------------------------------------------------------

def make_box(x,y, eps = 2.5e4):
    
    x = np.array(x)
    y = np.array(y)
    
    xmin, xmax, ymin, ymax = x.min()-eps, x.max()+eps, y.min()-eps, y.max()+eps
    
    x = np.arange(xmin, xmax+1, 5000.)
    y = np.arange(ymin, ymax+1, 5000.)     
    
    clow = np.array([x[::-1], np.ones_like(x)*ymin])
    cup = np.array([x, np.ones_like(x)*ymax])
    cl = np.array([np.ones_like(y)*xmin, y])
    cr = np.array([np.ones_like(y)*xmax, y[::-1]])
    return clow,cup,cl,cr

if catchement:
    try:
        x,y = np.loadtxt('../Catchments/catchment.txt').T
    except OSError:
        print('no catchment file found. ')
    clow,cup,cl,cr = make_box(x,y)
else:
    clow,cup,cl,cr = make_box(x,y, eps= 0.0)

plt.plot(x,y)
plt.plot(clow[0], clow[1], c= 'k')
plt.plot(cup[0], cup[1], c= 'k')
plt.plot(cl[0], cl[1], c= 'k')
plt.plot(cr[0], cr[1], c= 'k')
plt.savefig('catchment_box.png')
plt.close()

directory = 'boundaries'
if not os.path.isdir(directory):
    os.mkdir(directory)
    
for name, array in zip(['c1','c2', 'c3', 'c4'], [cup, cr, clow, cl]):
    np.savetxt(directory+'/'+name+'.txt', array.transpose())


sys.stdout.write('box domain created...\n\n')
