# ElmerIceCodes - MELT

### Authors

Cyrille Mosbeubx, Benoit Urruty, Justine Caillet, Fabien Gillet-Chaulet

## Intro 

This repository contains the MELT module for the Elmer/Ice software, which is a finite element software for modeling multiphysical problems related to ice.

The MELT module provides functionality for parameterizing ocean-induced melt processes. It includes various parameterization options. THe following parameterizations are implemented:

- [x] PICO
  - [ ] we should avoid calling the `USF_CondDistance.F90` in the BodyForce. The USFs should ve coded in the `DistanceModule.F90`
- [ ] Linear
- [ ] Quadratic 

## Contents
- `Melt.F90`: Main subroutine for calling different melt parameterizations.

- `PICO.F90`: Implementation of the PICO parameterization for melt processes. This module relies on the other module:
  - `DistanceModule.F90`: Module for computing distances to grounding lines and ice fronts.
- Other supporting files and modules such as the `USF_CondDistance.F90` and the `USF_CondFront.F90`. 

## Usage

1. Compile the code with the given Makefile that runs the following command: `elmerf90 DistanceModule.F90 PICO.F90 Melt.F90 -o Melt $(LFLAGS) $(IFLAGS)`. The flags should be adpated to your machine.

2. The solver Melt should be called in you `.sif` and applied on a 2D body (the basal boundary in 3D or the bulk in 2D). For PICO, we can give:

```f90
Solver 5
  Equation = "MELT"
  Variable = -nooutput dummy
  Procedure = "./MELT/Melt" "MeltSolver"

  data file = File "Path-to-DATA_boxmodelv2-2.nc"

  Bottom Surface Variable Name = String "Surface Name"
  Grounding Line Melt = Logical False
  PanAntarctic = Logical False

  Parameterisation Name = string "PICO"
  DistGL Name = string "distGL"
  DistIF Name = string "distIF"

  Exported variable 1 = -dofs 1 distGL
  Exported variable 2 = -dofs 1 FrontMask
  Exported variable 3 = -dofs 1 distIF
  Exported Variable 4 = -dofs 1 -elem Boxes
  Exported Variable 5 = -dofs 1 -elem Melt
End
```
The solver works both elmentaly and nodaly. Nodal variables can be forced with `Nodal Melt = Logical True`

```f90
  Exported Variable 4 = -dofs 1 -elem Boxes
  Exported Variable 5 = -dofs 1 -elem Melt
End
```

The front boundary should be signaled by, e.g.:

```f90
Boundary Condition 1
 Name = "calving_front"
 Target Boundaries(1) = 1

 !Needed for IceFrontLocation
 Ice Front = Logical True
End
```

> :warning: **Warning**: For now, the solver still requests some conditions in the bodyforce to run PICO:

```f90
  !Applied to basal boundary
  distGL = Real 0.0  ! used in Solver DistanceSolver1
  distGL Condition = Variable GroundedMask
    Real procedure "PICO/USF_CondDistance" "CondDistance"
  distIF = Real 0.0
  distIF Condition = Variable FrontMask
    Real procedure "PICO/USF_CondFront" "CondFront"
```


## Contributing

Contributions to this project are welcome. If you encounter any issues or have suggestions for improvements, please open an issue on the GitHub repository.

## License

This code is licensed under the GNU General Public License, version 2 or later. See the `LICENSE` file for more details.

## Contact

For any inquiries or questions, please contact:
- Author: Cyrille Mosbeux
- Email: cyrille.mosbeux@univ-grenoble-alpes.fr
