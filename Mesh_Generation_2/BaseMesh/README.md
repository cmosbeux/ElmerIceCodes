# 2D Mesh Generation with GMSH

This folder contains all the files and scripts necessary to create a 2D mesh using **GMSH**.

## Prerequisites

- Ensure that **GMSH** is installed on your machine.
- If you plan to convert the generated mesh to Elmer format, make sure **Elmer** is installed as well.

## File Structure

- `contour_1.txt` and `contour_2.txt` should be placed in the `../contours` directory. These files represent your two boundaries (inland and front).
  - Each file should contain two columns representing the (x, y) coordinates, with tab separation between the columns.
  
## Mesh Resolution

You can specify the mesh resolution in the `OPTIONS` file by defining a value for `lc1` in meters. The resolution controls the granularity of the mesh.

## Steps to Generate the Mesh

1. **Install GMSH**: Ensure that GMSH is installed and accessible via your terminal or command line.
2. **Prepare Contour Files**:
    - Place your boundary contours (inland and front) in the `../contours` directory.
    - Ensure the files are named `contour_1.txt` and `contour_2.txt`.
3. **Set Mesh Resolution**: Define the desired mesh resolution in the `OPTIONS` file by assigning a value to `lc1`.
4. **Run the Mesh Script**:
    - Execute the following command to generate the mesh:
      ```bash
      ./script_maillage.sh
      ```
5. **Convert to Elmer Format (Optional)**:
    - If Elmer is installed on your machine, you can directly convert the generated `.msh` file to Elmer format using the following command:
      ```bash
      ElmerGrid 14 2 Mesh.msh -autoclean
      ```

## Additional Information

- **GMSH Documentation**: For more details about GMSH, refer to the [GMSH official documentation](http://gmsh.info/doc/texinfo/gmsh.html).
- **Elmer Documentation**: Learn more about Elmer by visiting the [Elmer official documentation](https://www.csc.fi/web/elmer).

