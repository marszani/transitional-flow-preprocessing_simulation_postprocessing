Minor Code Modification Documentation

This document summarizes the minor modifications for the Langtry-Menter SST transitional turbulence model in ADflow, detailing the affected modules and key subroutines with code excerpts. The major model implementation is documented discursively elsewhere.

⸻

Subroutine turbSolveDDADI in module TurbAPI

Handling for langtryMenterSST added:

select case (turbModel)
...
  case (komegaWilcox, komegaModified, menterSST, langtryMenterSST, ktau)
    call unsteadyTurbSpectral(itu1, itu2)  ! Compute unsteady turbulence spectral terms

  case (menterSST, langtryMenterSST)
    call SST_block_residuals(.false.)     ! Disable SST residual blocks for initial solve
    call SSTSolve                          ! Call SST solver for turbulence variables
end select


⸻

Subroutine turbResidual in module TurbAPI

Included langtryMenterSST case:

select case (turbModel)
...
  case (menterSST, langtryMenterSST)
    call SST_block_residuals(.True.)  ! Enable SST residual blocks to compute turbulence residuals
end select


⸻

Subroutine computeEddyViscosity in module turbUtils

Added langtryMenterSST to call SST eddy viscosity calculation:

select case (turbModel)
...
  case (menterSST, langtryMenterSST)
    call SSTEddyViscosity(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd) ! Compute eddy viscosity for SST
  ! other cases...
end select


⸻

Subroutine setBCVarNamesTurb in module BCData

Added BC variable names for transition variables:

select case (turbModel)
...
  case (langtryMenterSST)
    bcVarNames(offset + 1) = cgnsTurbK            ! Turbulent kinetic energy
    bcVarNames(offset + 2) = cgnsTurbOmega        ! Specific dissipation rate
    bcVarNames(offset + 3) = cgnsTransitionGamma  ! Transition onset parameter
    bcVarNames(offset + 4) = cgnsTransitionReThetat ! Transition Reynolds number parameter
end select


⸻

Logical Function setBCVarTurb in module BCData

Freestream initialization for transition variables:

select case (turbModel)
...
  case (langtryMenterSST)
    ref(itu1) = pRef / rhoRef      ! Initialize turbulent kinetic energy reference
    ref(itu2) = ref(itu1) / nuRef  ! Initialize omega reference from k and viscosity
    ref(iTransition1) = 0          
    ref(iTransition2) = 0        
end select


⸻

Main residual loop for RANS equations in module masterRoutines

if (equations == RANSEquations) then
...
  select case (turbModel)
    case (menterSST, langtryMenterSST)
      call SST_block_residuals(.True.)  ! Include SST residuals for turbulence model solve
  end select
end if


⸻

Subroutine volSolNames in module outputMod

Added volume solution variable names for Langtry-Menter SST model:

case (langtryMenterSST)
...
  solNames(itu1) = cgnsTurbK                ! Name for turbulent kinetic energy solution variable
  solNames(itu2) = cgnsTurbOmega            ! Name for omega solution variable
  solNames(iTransition1) = cgnsTransitionGamma    ! Name for transition gamma solution variable
  solNames(iTransition2) = cgnsTransitionReThetat ! Name for transition ReTheta solution variable

...

if (volWriteResTurb) then

  select case (turbModel)
  ...
  case (langtryMenterSST)
      nn = nn + 1
      solNames(nn) = cgnsResK

      nn = nn + 1
      solNames(nn) = cgnsResOmega

      nn = nn + 1
      solNames(nn) = cgnsTransitionGamma

      nn = nn + 1
      solNames(nn) = cgnsTransitionReThetat
      


⸻

Subroutine isoSurfNames in module outputMod

Added iso-surface solution variable names for Langtry-Menter SST model:

...
if (isoWriteTurb) then

  select case (turbModel)
  ...
  case (langtryMenterSST)
      nn = nn + 1
      solNames(nn) = cgnsResK

      nn = nn + 1
      solNames(nn) = cgnsResOmega

      nn = nn + 1
      solNames(nn) = cgnsTransitionGamma

      nn = nn + 1
      solNames(nn) = cgnsTransitionReThetat
      
 ...

if (isoWriteResTurb) then

  select case (turbModel)          
  ...
  case (langtryMenterSST)
      nn = nn + 1
      solNames(nn) = cgnsResK

      nn = nn + 1
      solNames(nn) = cgnsResOmega

      nn = nn + 1
      solNames(nn) = cgnsTransitionGamma

      nn = nn + 1
      solNames(nn) = cgnsTransitionReThetat


⸻

Subroutine writeCGNSHeader in module outputMod

Added support for Langtry-Menter SST in CGNS metadata output:

turbulentTest: if (equations == RANSEquations) then
  ...
  case (menterSST, langtryMenterSST)
    call writeCGNSMenterSSTInfo(cgnsInd, base)


⸻

Subroutine blockResCore in module blockette

Enabled residual block computation for Langtry-Menter SST in RANS mode:

if (equations == RANSEquations .and. turbRes) then
  ...
  case (menterSST, langtryMenterSST)
    call SST_block_residuals(.True.)


⸻

Subroutine referenceState in module initializeFlow

Freestream turbulent variables initialization:

select case (turbModel)
....
  case (komegaWilcox, komegaModified, menterSST, langtryMenterSST)
    wInf(itu1) = 1.5_realType * uInf2 * turbIntensityInf**2  ! Calculate freestream turbulent kinetic energy
    wInf(itu2) = wInf(itu1) / (eddyVisInfRatio * nuInf)      ! Calculate freestream omega from k
end select


⸻

Subroutine checkMonitor in module inputParamRoutines

Added monitors for new turbulence variables:

...
case (langtryMenterSST)
  nMon = nMon + 4; nMonSum = nMonSum + 4
  monNames(nMon - 3) = cgnsL2ResK        ! Monitor residual of k
  monNames(nMon - 2) = cgnsL2ResOmega    ! Monitor residual of omega
  monNames(nMon - 1) = cgnsL2ResGamma    ! Monitor residual of transition gamma
  monNames(nMon)     = cgnsL2ResRethetat ! Monitor residual of transition ReTheta


⸻

Subroutine setEquationParameters in module inputParamRoutines

...
case (langtryMenterSST)
  nw = 9                ! Number of variables including transition scalars
  nt2 = 9               ! Total turbulence variable count
  iTransition1 = 8      ! Index for transition gamma variable
  iTransition2 = 9      ! Index for transition ReTheta variable
  kPresent = .true.     ! Turbulence kinetic energy present
  eddyModel = .true.    ! Eddy viscosity model active
  transitionModel = GammaRetheta  ! Use Gamma-ReTheta transition model


⸻

Subroutine readTurbvar in module variableReading

Added langtryMenterSST to support reading turbulence variables using readTurbKwType:

...
case (menterSST, langtryMenterSST) ! this is a hack and should be implemented properly
  call readTurbKwType(nTypeMismatch)


⸻

transitionModel in module inputPhysics

Added transitionModel variable to support transition model:

! transitionModel      Which transition Model to use
....
integer(kind=intType) :: turbModel, cpModel, turbProd, transitionModel


⸻

gammaretheta in module constants

Added gammaretheta as a transition model:

integer(kind=intType), parameter :: &
    noTransitionModel = 0, & 
    gammaretheta = 1


⸻

C Interface (libadflowmodule.c)

...
{"langtrymentersst",0,{{-1}},NPY_INT},
...
static void f2py_setup_constants(char *maxstringlen,..... *langtrymentersst,char *v2f,char .....
...
f2py_constants_def[i_f2py++].data = langtrymentersst;
...
"Fortran 90/95 modules:\n""  constants --- maxstringlen,maxcgnsnamelen,..,langtrymentersst,....."


⸻

Python Interface (pyADflow.py)

"turbulencemodel": {
    ...
    "sa": self.adflow.constants.spalartallmaras,        # Spalart-Allmaras model
    "menter sst": self.adflow.constants.mentersst,      # Menter SST model
    "langtry menter sst": self.adflow.constants.langtrymentersst,  # Langtry-Menter SST model
    ...
},
"eddyVisInfRatio": [float, 0.009],     # Default freestream eddy viscosity ratio, but can be changed in run script 
"turbIntensityInf": [float, 0.001],    # Default freestream turbulence intensity, but can be changed in run script 


