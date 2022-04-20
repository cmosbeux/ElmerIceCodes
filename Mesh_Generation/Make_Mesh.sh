#Make some cleaning

#Run python scripts
python make_box_domain.py
python Makegeo_box.py

#Run gmsh and ElmerGrid
name=$(cat mesh_name.txt)
gmsh -2 "$name.geo"

#ElmerGrid 14 2 "$name.msh" -autoclean
echo $name
rm "$name.geo" "mesh_name.txt"
