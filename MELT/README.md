# ElmerIceCodes - MELT

### Authors

Cyrille Mosbeubx, Benoit Urruty, Justile Caillet, Fabien Gillet-Chaulet

## Intro 

This repository contains the MELT module for the Elmer/Ice software, which is a finite element software for modeling multiphysical problems related to ice and snow.

The MELT module provides functionality for parameterizing melt processes in ice and snow. It includes various parameterization options and solvers for computing distances to grounding lines and ice fronts.

## Contents

- `DistanceModule.F90`: Module for computing distances to grounding lines and ice fronts.
- `PICO.F90`: Implementation of the PICO parameterization for melt processes.
- `Melt.F90`: Main subroutine for calling different melt parameterizations.
- Other supporting files and modules.

## Usage

1. Compile the code with the given Makefile that runs the following command: `elmerf90 DistanceModule.F90 PICO.F90 Melt.F90 -o Melt $(LFLAGS) $(IFLAGS)`. The flags should be adpated to your machine.

2. The solver Melt should be called in you `.sif`. For PICO, we can give:

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
  Exported Variable 4 = -dofs 1 Boxes
  Exported Variable 5 = -dofs 1 Melt
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
- Email: cmosbeux@univ-grenoble-alpes.fr
