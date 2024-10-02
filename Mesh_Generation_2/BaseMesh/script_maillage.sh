#!/bin/bash

# Run Python script to create .geo file
python MakeGeo.py

# Generate mesh using Gmsh
gmsh -2 Contour.geo

# Rename the mesh file
mv Contour.msh Mesh.msh

# Output file and potential conversion
echo ' '
echo 'Your output file is "Mesh.msh". You can convert it to Elmer format by running:'
echo 'ElmerGrid 14 2 Mesh.msh -autoclean'
echo '———'

# Ask if they want to convert the mesh to Elmer format
read -p "Do you want to convert to Elmer format now? (Y/N): " answer

# Convert to uppercase to handle lowercase 'y'
answer=$(echo "$answer" | tr '[:lower:]' '[:upper:]')

# Check if the answer is 'Y'
if [ "$answer" == "Y" ]; then
    echo "Converting to Elmer format..."
    # Run ElmerGrid with the autoclean option
    ElmerGrid 14 2 Mesh.msh -autoclean
    echo "Conversion complete!"
else
    echo "Conversion aborted. You can run it manually later with:"
    echo "ElmerGrid 14 2 Mesh.msh -autoclean"
fi

