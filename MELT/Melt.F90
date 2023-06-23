!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Mosbeux
! *  Email:   cmosbeux@univ-grenoble-alpes.fr
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 31 May 2023
! *
! *****************************************************************************/

  
!-------------------------------------------------------------------------------
!> Subroutine for calling different melt parameterizations 
!   * PICO
!   * Quadratic 
!------------------------------------------------------------------------------
SUBROUTINE MeltSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Distance
  USE PICO

  IMPLICIT NONE
  
  TYPE(Model_t)  :: Model
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t), TARGET :: Solver
  
  REAL(KIND=dp) :: dt
  LOGICAL ::  TransientSimulation
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='Melt', DistGLName, DistIFName, ParamName


  Params => GetSolverParams()
  ParamName = ListGetString(Params, 'Parameterisation Name', UnFoundFatal=.TRUE.)

  CALL INFO(Trim(SolverName),'Starting Melt Parameterisation.', Level=5 )
  
  SELECT CASE (TRIM(ParamName))

    CASE ('pico')
      CALL INFO(Trim(SolverName),'Choice : PICO', Level=5 )
      !------------------------------------
      !1 - Compute distance to GL and IF
      !------------------------------------
      ! Grounding line
      DistGLName = ListGetString(Params, 'DistGL Name', UnFoundFatal=.TRUE.)
      CALL DistanceSolver2( Model,Solver,dt,TransientSimulation, DistGLName )

      ! Ice Front
      CALL IceFrontMask( Model,Solver,dt,TransientSimulation )

      DistIFName = ListGetString(Params, 'DistIF Name', UnFoundFatal=.TRUE.)
      CALL DistanceSolver2( Model,Solver,dt,TransientSimulation, DistIFName )
     
      !---------------------------------------
      !2 - Run the BoxModel Solver (PICO.F90)
      !---------------------------------------
      CALL boxmodel_solver(Model,Solver,dt,TransientSimulation)

    CASE DEFAULT
        CALL FATAL(SolverName,'Parameterisation does not exist.')
  END SELECT  

END SUBROUTINE MeltSolver