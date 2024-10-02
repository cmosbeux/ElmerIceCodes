"""Created on Thu Sep 14 11:53:45 2017

Python Script to create .geo files adapted from the matlab version of
F. Gillet-Chaulet, O. Gagliardini, May 2011

@author: cmosbeux
"""

import numpy as np

def read_value_from_file(filename, key='lc1'):
    """Read the value for the specified key (e.g., 'lc1') from a text file.
    The file can have a header, and the key-value pair can be on any line.
    """
    with open(filename, 'r') as file:
        for line in file:
            # Strip leading/trailing spaces and newlines
            line = line.strip()
            # Check if the line starts with the key (e.g., 'lc1=')
            if line.startswith(f'{key}='):
                # Extract and return the value after the '=' symbol
                try:
                    value = float(line.split('=')[1])
                    return value
                except ValueError:
                    raise ValueError(f"Invalid value for {key} in file.")
    
    # If the key is not found, raise an error
    raise ValueError(f"Key '{key}' not found in file.")



print('BUILD GEO FILE')
print('--------------')

# Read lc1 from a text file
lc1_file = 'OPTIONS'
lc1 = read_value_from_file(lc1_file)
print(f'lc1 value loaded from {lc1_file}: {lc1}')


print('BUILD GEO FILE')
print('--------------')

fname='Contour.geo'
print('loading contour files...')
Contour_dir = '../Contours/'
A1=np.loadtxt(Contour_dir+'Contour_2.txt')
A2=np.loadtxt(Contour_dir+'Contour_1.txt')
print('done.')
#%%
print('building geo files...')
A1 = A1.transpose()
A2 = A2.transpose()

array_A1=np.arange(0,np.size(A1[0]),1)
array_A2=np.arange(0,np.size(A2[0]),1)

long_A1=len(array_A1)
long_A2=len(array_A2)

A=np.concatenate((A1.transpose(),A2.transpose())).transpose()

with open(fname,  'w') as fid1:
    fid1.write('Mesh.Algorithm=5; \n')
    As=np.size(A[0])
    n=0
    leng=np.arange(0,As,1)
    for ii in leng:
        n=n+1
        fid1.write(''.join(('Point(', str(n), ')', '={', '%14.7e'%A[0][ii],',', '%14.7e'%A[1][ii], ',', '0.0',',','%0d'%lc1,'};\n')))

    #spline for land
    fid1.write('Spline(1)={')
    n_spline=0
    fid1.write(','.join( '%0d'%ii for ii in np.arange(1, long_A1+1)))
    fid1.write('}; \n')

    #spline for ice front
    fid1.write('Spline(2)={')
    n_spline=0
    fid1.write(','.join( '%0d'%ii for ii in np.arange(long_A1,long_A1+long_A2+1)))
    fid1.write(',1}; \n')

    #Splines make part of a line loop
    fid1.write('Line Loop(5)={1,2}; \n')

    #Line Loop define
    fid1.write('Plane Surface(6)={5}; \n')

    fid1.write('Physical Line(7) = {1}; \n')
    fid1.write('Physical Line(8) = {2}; \n')

    fid1.write('Physical Surface(9)={6}; \n')

fid1.close()

print('done.')
print('\n')

