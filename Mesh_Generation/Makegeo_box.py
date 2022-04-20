"""Created on Thu Sep 14 11:53:45 2017

Python Script to create .geo files adapted from the matlab version of
F. Gillet-Chaulet, O. Gagliardini, May 2011

@author: cmosbeux
"""

import numpy as np
import os
import sys

#%% 
#Read inputs
sys.stdout.write('\nReading file parameters...\n\n')
f = open('param.txt', 'r') 
txt = f.readlines()
for line in txt:
    if 'lc1' in line:
        lc = line.split('=')[-1].rstrip()
        sys.stdout.write('resolution (km) = %s\n' % lc)
        lc = int(float(lc))
    elif 'name' in line:
        fname = line.split('=')[-1].rstrip()
        sys.stdout.write('mesh name = %s\n' % fname)
f.close()

fname = fname+'_%dkm' % (lc/1e3)
with open('mesh_name.txt', 'w') as f:
    f.write(fname)

#%%
fname=fname+'.geo' 

A1=np.loadtxt('boundaries/c1.txt').transpose()
A2=np.loadtxt('boundaries/c2.txt' ).transpose()
A3=np.loadtxt('boundaries/c3.txt'  ).transpose()
A4=np.loadtxt('boundaries/c4.txt' ).transpose()

array_A1=np.arange(0,np.size(A1[0]),1)
array_A2=np.arange(0,np.size(A2[0]),1)
array_A3=np.arange(0,np.size(A3[0]),1)
array_A4=np.arange(0,np.size(A4[0]),1)

long_A1=len(array_A1)
long_A2=len(array_A2)
long_A3=len(array_A3)
long_A4=len(array_A4)

A=np.concatenate((A1.transpose(),A2.transpose(),A3.transpose(),A4.transpose())).transpose()

with open(fname,  'w') as fid1:
    fid1.write('Mesh.Algorithm=5; \n')
    As=np.size(A[0])
    n=0
    leng=np.arange(0,As,1)
    for ii in leng:
        n=n+1
        fid1.write(''.join(('Point(', str(n), ')', '={', '%14.7e'%A[0][ii],',', '%14.7e'%A[1][ii], ',', '0.0',',','%0d'%lc,'};\n')))
    
    
    #spline for surface
    fid1.write('Spline(1)={')
    n_spline=0
    fid1.write(','.join( '%0d'%ii for ii in np.arange(1, long_A1+1)))
    fid1.write('}; \n')
    
    #spline for front
    fid1.write('Spline(2)={')
    n_spline=0
    fid1.write(','.join( '%0d'%ii for ii in np.arange(long_A1,long_A1+long_A2+1)))
    fid1.write('}; \n')
    
    #spline for base
    fid1.write('Spline(3)={')
    n_spline=0
    fid1.write(','.join( '%0d'%ii for ii in np.arange(long_A1+long_A2,long_A1+long_A2+long_A3+1)))
    fid1.write('}; \n')
    
    #spline for input
    fid1.write('Spline(4)={')
    n_spline=0
    fid1.write(','.join( '%0d'%ii for ii in np.arange(long_A1+long_A2+long_A3,long_A1+long_A2+long_A3+long_A4+1)))
    fid1.write(',1}; \n')
    
    #Splines make part of a line loop
    fid1.write('Line Loop(5)={1,2,3,4}; \n')
    #Line Loop define
    fid1.write('Plane Surface(6)={5}; \n')
    
    fid1.write('Physical Line(7) = {1}; \n')
    fid1.write('Physical Line(8) = {2}; \n')
    fid1.write('Physical Line(9) = {3}; \n')
    fid1.write('Physical Line(10) = {4}; \n')
    
    fid1.write('Physical Surface(11)={6}; \n')

      
fid1.close()

sys.stdout.write('\nGeo file created...\n\n')

#%%
#run gmsh
# gmsh_command = "gmsh -2 %s" % (fname)
# os.system(gmsh_command)

