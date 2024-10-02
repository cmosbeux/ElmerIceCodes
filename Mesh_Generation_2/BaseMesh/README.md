This folder contains all the files and scripts to create a 2D mesh using GMSH. 

Steps:
	1) Make sure that gmsh is installed  on your machine
	2) You have to put your contour 1 and 2 (your two boundaries (inland and front)) in the '../contours' directory
		contour_1.txt and contour_2.txt should contain two columns for (x,y) coordinates
		with a tab separation between columns
	3) You can choose your mesh resolution in the OPTIONS File by defining a value for lc1 in meters
	
	3) Run script_maillage.sh
	4) If you have Elmer on your machine, you can directly convert your msh file to an ElmerFormat
