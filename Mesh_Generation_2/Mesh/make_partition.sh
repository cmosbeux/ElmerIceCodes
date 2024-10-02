#!/bin/bash

n=24

ElmerGrid 2 2 MESH_3/ -autoclean -metis "$n"  

mkdir mesh_"$n"

mv MESH_3/partitioning."$n" mesh_"$n"/.
