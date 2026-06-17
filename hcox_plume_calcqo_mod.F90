!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_plume_mod.F90
!
! !DESCRIPTION: Module HCOX\_Plume\_Mod implements a thermodynamic plume rise
! model for wildfire emissions in HEMCO. It reads 2D fire emissions (CO),
! fire energy (Qo), and burned area (BurnAreaTot) from HEMCO, and then
! distributes them using the Anderson et al (2011) scheme.
!\\
!\\
! Required HEMCO_Config.rc base fields:
!   Fire_CO          - 2D CO emission rate [kg/m2/s]  (NetCDF var: ECO)
!   Fire_TFC         - 2D fire energy      [kg/m2]    (NetCDF var: ??)
!   Fire_BurnAreaTot - 2D burned area      [km2]      (NetCDF var: BurnAreaTot)
!   Fire_GrowthArea  - 2D new area burned  [km2]      (NetCDF var: ???) 
! 
!  Other things that have been changed:
!       - added init run and final to hcox_driver_mod.F90
!       - hcox_state_mod.F90
!       - add to CmakeList.txt
! !INTERFACE:
!
MODULE HCOX_Plume_Mod
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCOX_TOOLS_MOD
  USE HCOX_State_MOD, ONLY : Ext_State
  USE HCO_State_MOD,  ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_Plume_Run
  PUBLIC :: HCOX_Plume_Init
  PUBLIC :: HCOX_Plume_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DistributeEmissionsDensity
  PRIVATE :: DistributeEmissionsUniform
  PRIVATE :: DistributeEmissionsGaussian
!
! !REVISION HISTORY:
!  10 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE PARAMETERS:
!
  ! Physical constants for the plume rise algorithm
  REAL(sp), PARAMETER :: GRAV_CF   = 9.80025_sp  ! gravity [m/s^2]
  REAL(sp), PARAMETER :: RGASD_CF  = 287.05_sp   ! dry air gas constant [J/kg/K]
  REAL(sp), PARAMETER :: CPD_CF    = 1004.0_sp   ! dry air specific heat [J/kg/K]
  REAL(sp), PARAMETER :: LD_CF     = -0.0098_sp  ! dry adiabatic lapse rate [K/m]
  REAL(sp), PARAMETER :: H         = 18000000_sp ! Heat of Combustion for dry wood [J/kg]
  REAL(sp), PARAMETER :: Fire_Eff_Def = 0.2_sp   ! Default effieciency of fire
  INTEGER,  PARAMETER :: IMAX_CF   = 15000       ! max plume height and max iteration [m]
  INTEGER,  PARAMETER :: FOREST_ID(18) = [3, 4, 5, 6, 21, 22, 23, 24, 25, 26, 27, 32, 33, 43, 54, 60, 61, 62]
  INTEGER,  PARAMETER :: CROP_ID(8) = [29, 30, 31, 35, 36, 37, 38, 41]
  INTEGER,  PARAMETER :: GRASS_ID(10) = [2, 10, 41, 42, 46, 47, 51, 52, 59, 64]

  ! These are the first ID sets I generated, but there are only 73 Ids in the ROOT directory, so i guess it is the 1992 version
  ! instead of the 2001 one??
  !INTEGER,  PARAMETER :: FOREST_ID(20) = [3, 4, 5, 6, 21, 22, 23, 24, 25, 26, 27, 43, 54, 60, 61, 62, 77, 78, 95, 96]
  !INTEGER,  PARAMETER :: CROP_ID(11) = [30, 31, 35, 36, 37, 38, 39, 76, 92, 93, 94]
  !INTEGER,  PARAMETER :: GRASS_ID(10) = [2, 10, 41, 42, 46, 47, 51, 52, 59, 64]
!
! !MODULE VARIABLES:
!
  ! MyInst is the extension-specific derived type holding all per-instance state.
  TYPE :: MyInst
    INTEGER                         :: Instance
    INTEGER                         :: ExtNr         = -1
    INTEGER                         :: nSpc          =  0
    INTEGER,  ALLOCATABLE           :: SpcIDs(:)
    REAL(sp), ALLOCATABLE           :: SpcScl(:)
    CHARACTER(LEN=31), ALLOCATABLE  :: SpcNames(:)
    CHARACTER(LEN=61), ALLOCATABLE  :: SpcScalFldNme(:)
    TYPE(MyInst), POINTER           :: NextInst => NULL()

    !=================================================================
    ! Module specific variables/arraysdata pointers below
    ! Make pointers to the 3 inputs we need to read
    !=================================================================
    REAL(hp), POINTER :: Fire_CO(:,:)          => NULL()  ! CO [kg/m2/s]
    REAL(hp), POINTER :: Fire_GrowthArea(:,:)  => NULL()  ! Area burned by fuel  [km2]
    REAL(hp), POINTER :: Fire_TFC(:,:)         => NULL()  ! Total fire consumed [kg/m2]
    REAL(hp), POINTER :: Fire_BurnAreaTot(:,:) => NULL()  ! burned area [km2]
    REAL(hp), POINTER :: LANDTYPE(:,:,:)       => NULL()  ! land type
    REAL(sp)          :: Fire_Eff = 0.2_sp                ! fire efficiency default value
    INTEGER           :: DistMethod = 1
    LOGICAL           :: FirstRun = .TRUE.                ! Flag for reading landtype values
    
  END TYPE MyInst

  ! Pointer to all instances
  TYPE(MyInst), POINTER :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Plume_Run
!
! !DESCRIPTION: For each timestep, reads 2D fire fields and computes 3D
! wildfire CO emissions by distributing surface fluxes vertically up to
! the Anderson et al. (2011) thermodynamic plume height.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Plume_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER       :: ExtState    ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! HEMCO state
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  2024 - J. Landers - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MyInst), POINTER :: Inst => NULL()
    CHARACTER(LEN=255)    :: MSG, LOC                              ! hemco things
    CHARACTER(LEN=12)     :: LandName                              ! CHAR to read Landtypes
    INTEGER               :: T                                     ! INT to read Landtypes
    INTEGER               :: I, J, L, N, JJ                        ! loop variables
    INTEGER               :: NX, NY, NZ                            ! number of grid cells
    INTEGER               :: ref1, ref2, ref3, ref4, ref5          ! Levels correspoinding to lapse rates
    INTEGER               :: ii                                    ! plume height top and N layers
    REAL(sp)              :: FuelMult                              ! store values of land type coverage
    REAL(sp)              :: QPlume_ij, BurnArea_ij, CO_ij         ! fire values at ij
    REAL(sp)              :: GrowthArea_ij, TFC_ij                 ! fire values at ij
    REAL(sp)              :: Le, Le1, Le2, Le3, Le4                ! different lapse rate values
    REAL(sp)              :: tsurf, ps, area                       ! met/area at each fire
    REAL(sp)              :: dz, ddz, q, Qatm, mplume, work        ! compute variables
    REAL(sp), ALLOCATABLE :: Ph(:)                                 ! pressure profile [hPa]
    REAL(sp), ALLOCATABLE :: Zh(:)                                 ! altitude profile [m]
    REAL(sp), ALLOCATABLE :: Rho(:)                                ! air density profile [kg/m3]
    REAL(sp), ALLOCATABLE :: EmisCol(:)                            ! Single column of emissions
    REAL(sp), ALLOCATABLE :: Emis3D(:,:,:)                         ! 3D CO emissions [kg/m2/s]
    REAL(sp), ALLOCATABLE :: ZPlume2D(:,:)                         ! save plume top height [m]
    REAL(sp), ALLOCATABLE :: save_QPlume(:,:), save_Landtype(:,:)  !diagn

    !=================================================================
    ! HCOX_Plume_Run begins here!
    !=================================================================
    LOC = 'HCOX_Plume_Run (hcox_plume_mod.F90)'

    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get pointer to this instance
    CALL InstGet( ExtState%Plume, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find Plume instance Nr. ', ExtState%Plume
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Pull grid dimensions from HcoState
    NX = HcoState%NX
    NY = HcoState%NY
    NZ = HcoState%NZ

    ! Allocate fire field arrays on first call. We use HCO_EvalFld to read
    ! netcdf data for emissions, Qo and Burned area. HCO_EvalFld derives the
    ! grid size from the input array, so we need to allocate the arrays so it gets
    ! the right size
    IF ( .NOT. ASSOCIATED(Inst%Fire_CO) ) THEN
       ALLOCATE( Inst%Fire_CO(NX,NY),          &
                 Inst%Fire_GrowthArea(NX,NY),  &
                 Inst%Fire_TFC(NX,NY),         &
                 Inst%Fire_BurnAreaTot(NX,NY) )
       Inst%Fire_CO          = 0.0_hp
       Inst%Fire_GrowthArea  = 0.0_hp
       Inst%Fire_TFC         = 0.0_hp
       Inst%Fire_BurnAreaTot = 0.0_hp
    ENDIF

    ! Read current 2D fire fields from HEMCO data infrastructure
    ! Use HCO_EvalFld to read Fire_CO from HcoState, and write it to our current instance variable
    CALL HCO_EvalFld( HcoState, 'Fire_CO',          Inst%Fire_CO,          RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot get Fire_CO field', RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'Fire_GrowthArea',          Inst%Fire_GrowthArea,          RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot get Fire_GrowthArea field', RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'Fire_BurnAreaTot', Inst%Fire_BurnAreaTot, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot get Fire_BurnAreaTot field', RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'Fire_TFC', Inst%Fire_TFC, RC )
    IF (RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot get Fire_TFC field', RC, THISLOC=LOC )
       RETURN
    ENDIF

    IF ( Inst%FirstRun ) THEN
        ALLOCATE( Inst%LANDTYPE(NX, NY, 73) )
        Inst%LANDTYPE = 0.0_hp
        DO T = 0, 72
            WRITE(LandName, '(A8,I2.2)') 'LANDTYPE', T
            CALL HCO_EvalFld( HcoState, TRIM(LandName), Inst%LANDTYPE(:,:,T+1), RC )
            IF ( RC /= HCO_SUCCESS ) THEN
                    WRITE(MSG,*) 'Cannot get land type field: ', TRIM(LandName)
                    CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                    RETURN
            ENDIF
        ENDDO
        Inst%FirstRun = .FALSE.
    ENDIF


    ! Allocate arrays some arrays
    ALLOCATE( Ph(NZ), Zh(NZ), Rho(NZ), EmisCol(NZ), Emis3D(NX,NY,NZ), ZPlume2D(NX,NY) )
    ALLOCATE( save_QPlume(NX,NY), save_Landtype(NX,NY) )
    Emis3D   = 0.0_sp
    ZPlume2D = 0.0_sp

    !=================================================================
    ! Start looping through boxes. 
    ! L=1 = surface, L=NZ = top
    !=================================================================
    DO J = 1, NY
    DO I = 1, NX
       ! Preliminary stuff setting up values before calculations

       ! Read fire parameters parameters for this box
       TFC_ij   = REAL( Inst%Fire_TFC(I,J),          sp )  ! Total Fuel consumed [kg/m2]
       GrowthArea_ij = REAL( Inst%Fire_GrowthArea(I,J),  sp )  ! Fire area burned area [km2]
       BurnArea_ij = REAL( Inst%Fire_BurnAreaTot(I,J), sp ) !Total burned area for fire plume [km2]
       CO_ij   = REAL( Inst%Fire_CO(I,J),           sp )  ! CO flux [kg/m2/s]

       FuelMult = 1.0_sp
       DO JJ  = 1, SIZE(CROP_ID)
           IF (Inst%LANDTYPE(I,J, CROP_ID(JJ) + 1 ) > 0.0_sp ) THEN
               FuelMult = 0.01_sp
               EXIT
           ENDIF
       ENDDO
       DO II = 1, SIZE(GRASS_ID)
           IF (Inst%LANDTYPE(I,J,GRASS_ID(II) + 1 ) > 0.0_sp ) THEN
               FuelMult = 0.05_sp
               EXIT
           ENDIF
       ENDDO
       DO II = 1, SIZE(FOREST_ID)
           IF ( Inst%LANDTYPE(I,J,FOREST_ID(II) + 1 ) > 0.0_sp ) THEN
               FuelMult = 2.0_sp
               EXIT
           ENDIF
       ENDDO
       save_Landtype(I,J) = FuelMult
           
       ! Energy is Total fuel consumed [kg/m2] * burned area [km2] * Heat of combustion [J/kg] * Fire efficiency * area conversion
       ! [m2/km2]
       QPlume_ij = TFC_ij * GrowthArea_ij * H * Inst%Fire_Eff * FuelMult * 1.0e6_sp
       save_QPlume(I,J) = QPlume_ij

       ! Skip columns with no fire emissions or fire energy ??? WRONG?? FIX???
       IF ( QPlume_ij <= 0.0_sp .OR. CO_ij <= 0.0_sp .OR. BurnArea_ij <= 0.0_sp ) CYCLE

       ! Build Pressure Array as midpoint of edge pressures
       ! Airdensity just from AIRDEN
       DO L = 1, NZ
          Ph(L) = ( REAL( HcoState%Grid%PEDGE%Val(I,J,L), sp ) + &     
                  REAL( HcoState%Grid%PEDGE%Val(I,J,L+1), sp ) ) / 2.0_sp
          Rho(L) = REAL( ExtState%AIRDEN%Arr%Val(I,J,L), sp)
       ENDDO

       ! Zh is altitude at the BOTTOM of the grid box
       Zh(1) = 0.0_sp
       DO L = 2, NZ
          Zh(L) = Zh(L-1) + REAL(HcoState%Grid%BXHEIGHT_M%Val(I,J,L-1), sp)
       ENDDO

       ! Find the vertical layer for 800, 700, 500 and 200 hPa to calculate lapse rates later
       ref1 = 1   ! surface Should delete ref1  and/or switch ref2->ref1, ref3->ref2 etc.
       ref2 = NZ  ! defaults in case threshold not found in column
       ref3 = NZ
       ref4 = NZ
       ref5 = NZ

       DO L = 2, NZ
          IF ( Ph(L) <= 85000.0_sp ) THEN; ref2 = L; EXIT; ENDIF
       ENDDO
       DO L = 3, NZ
          IF ( Ph(L) <= 70000.0_sp ) THEN; ref3 = L; EXIT; ENDIF
       ENDDO
       DO L = 4, NZ
          IF ( Ph(L) <= 50000.0_sp ) THEN; ref4 = L; EXIT; ENDIF
       ENDDO
       DO L = 5, NZ
          IF ( Ph(L) <= 20000.0_sp ) THEN; ref5 = L; EXIT; ENDIF
       ENDDO

       ! Environmental lapse rates at reference levels
       ! Le = (T_upper - T_surface) / (Z_upper - Z_surface)
       ! Zh is the altitude at the bottom of the box. I assume TK is the
       ! midpoint temperature, so we add 0.5 of box height of refX to find the
       ! altitude at the midpoint (technically T2M is the temp at 2m but who
       ! cares
       Le1 = 0.0_sp; Le2 = 0.0_sp; Le3 = 0.0_sp; Le4 = 0.0_sp
       Le1 = ( REAL(ExtState%TK%Arr%Val(I,J,ref2), sp)           &
             - REAL(ExtState%T2M%Arr%Val(I,J), sp) )          &
           / ( Zh(ref2) + 0.5_sp*REAL(HcoState%Grid%BXHEIGHT_M%Val(I,J,ref2), sp) )   ! surface to 850 hPa
       Le2 = ( REAL(ExtState%TK%Arr%Val(I,J,ref3), sp)           &
             - REAL(ExtState%T2M%Arr%Val(I,J), sp) )          &
           / ( Zh(ref3) + 0.5*REAL(HcoState%Grid%BXHEIGHT_M%Val(I,J,ref3), sp ) )   ! surface to 700 hPa
       Le3 = ( REAL(ExtState%TK%Arr%Val(I,J,ref4), sp)           &
             - REAL(ExtState%T2M%Arr%Val(I,J), sp) )          &
           / ( Zh(ref4) + 0.5*REAL(HcoState%Grid%BXHEIGHT_M%Val(I,J,ref4), sp) )  
       Le4 = ( REAL(ExtState%TK%Arr%Val(I,J,ref5), sp)           &
             - REAL(ExtState%T2M%Arr%Val(I,J), sp) )          &
           / ( Zh(ref5) + 0.5*REAL(HcoState%Grid%BXHEIGHT_M%Val(I,J,ref5), sp) )   ! surface to 200 hPa


       ! Preamble calculations done, now start the Anderson algorithm

       ! Initalize initial variables
       tsurf  = REAL( ExtState%T2M%Arr%Val(I,J), sp )     ! surface-layer T [K]
       ps     = ExtState%PSC2_WET%Arr%Val(I,J) * 100.0_sp  ! hPa -> Pa
       area   = BurnArea_ij * 1.0e6_sp                     ! km2 -> m2

       dz     = 1000.0_sp   ! initial height guess [m]
       ddz    = 1000.0_sp   ! bisection step [m]
       Le     = Le1         ! Start at lowest lapse rate
       ii     = 0
       mplume = 0.0_sp

       DO WHILE ( (ii < IMAX_CF) .AND. (ABS(ddz) >= 1.0_sp) .AND. (dz <= IMAX_CF) )

          ! Energy per unit mass to heat plume to dry adiabat (Anderson eqn 8)
          q = -0.5_sp * CPD_CF * LD_CF * dz &
            * LOG( 1.0_sp + dz * (Le - LD_CF) / tsurf )


          ! Column mass of cylindrical plume (Anderson eqn 16)
          IF ( Le /= 0.0_sp ) THEN
             work = ( 1.0_sp + Le * dz / tsurf )**( -GRAV_CF / Le / RGASD_CF )
          ELSE
             ! Isothermal: average limiting cases Le +/- 0.5
             work = 0.5_sp * (                                                      &
                    ( 1.0_sp + Le * dz / tsurf )**( -GRAV_CF / (Le+0.5_sp) / RGASD_CF ) + &
                    ( 1.0_sp + Le * dz / tsurf )**( -GRAV_CF / (Le-0.5_sp) / RGASD_CF ) )
          ENDIF
          mplume = ps / GRAV_CF * area * ( 1.0_sp - work )

          ! Total energy to heat column to dry adiabat (Anderson eqn 17)
          Qatm = q * mplume

          ! Iterate dz
          IF ( Qatm > QPlume_ij ) THEN
             ddz = 0.5_sp * ddz
             dz  = dz - ddz
          ELSE
             dz  = dz + ddz
          ENDIF

          ! Select lapse rate segment appropriate for current dz
          IF      ( dz > 2000.0_sp .AND. dz <= 4000.0_sp ) THEN
             Le = Le2
          ELSE IF ( dz > 4000.0_sp .AND. dz <= 7000.0_sp ) THEN
             Le = Le3
          ELSE IF ( dz > 7000.0_sp ) THEN
             Le = Le4
          ELSE
             Le = Le1
          ENDIF

          ii = ii + 1
       ENDDO

       ! Set negative plume height to 0
       IF ( dz < 0.0_sp ) dz = 0.0_sp

       ! Save plume top height for diagnostic output
       ZPlume2D(I,J) = dz 

       ! Plume height calculation done
       ! Distribute Emissions
       IF ( Inst%DistMethod == 1 ) THEN
           CALL DistributeEmissionsUniform( dz, Zh, CO_ij, EmisCol)
       ELSE IF (Inst%DistMethod == 2 ) THEN
           CALL DistributeEmissionsDensity( dz, Zh, Rho, CO_ij, EmisCol )
       ELSE IF (Inst%DistMethod == 3 ) THEN
           CALL DistributeEmissionsGaussian(dz, Zh, CO_ij, EmisCol )
       ENDIF

       Emis3D(I,J,:) = EmisCol

    ENDDO  ! I
    ENDDO  ! J

    ! Contribute 3D emissions to HEMCO flux arrays for each species
    DO N = 1, Inst%nSpc
       CALL HCO_EmisAdd( HcoState, Emis3D, Inst%SpcIDs(N), RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'HCO_EmisAdd error: plume rise ' 
          CALL HCO_ERROR( MSG, RC )
          ! Deallocate arrays?
          RETURN
       ENDIF
    ENDDO

    ! Write diagnosed plume top heights to the ZPlume_Calc diagnostic
    CALL Diagn_Update( HcoState, cName='ZPlume', Array2D=ZPlume2D, RC=RC )
    CALL Diagn_Update( HcoState, cName='save_QPlume', Array2D=save_QPlume, RC=RC )
    CALL Diagn_Update( HcoState, cName='save_Landtype', Array2D=save_Landtype, RC=RC )
    IF (RC /= HCO_SUCCESS ) THEN
       MSG = 'Diagn_Update error: plume rise'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Cleanup
    DEALLOCATE( Ph, Zh, Rho, EmisCol, Emis3D, ZPlume2D )
    DEALLOCATE( save_QPlume, save_Landtype )
    Inst => NULL()

    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_Plume_Run
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Plume_Init
!
! !DESCRIPTION: Initialize thermodynamic plume rise HEMCO extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Plume_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,  ONLY : GetExtNr
    USE HCO_ExtList_Mod,  ONLY : GetExtOpt
    USE HCO_ExtList_Mod,  ONLY : GetExtSpcVal
    USE HCO_STATE_MOD,    ONLY : HCO_GetExtHcoID
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! HEMCO state
    INTEGER,          INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller    - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MyInst), POINTER :: Inst => NULL()
    INTEGER               :: ExtNr, N
    LOGICAL               :: Found
    CHARACTER(LEN=255)    :: MSG, LOC
    CHARACTER(LEN=255)    :: DistMethodStr   ! Holds DistMethod from HEMCO_Config.rc

    !=================================================================
    ! HCOX_Plume_Init begins here!
    !=================================================================
    LOC = 'HCOX_Plume_Init (hcox_plume_mod.F90)'

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Create instance and link to ExtState%Plume
    CALL InstCreate( ExtNr, ExtState%Plume, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot create Plume instance', RC )
       RETURN
    ENDIF

    ! Get HEMCO species IDs for species declared in HEMCO_Config.rc under this extension
    CALL HCO_GetExtHcoID( HcoState, ExtNr, Inst%SpcIDs, Inst%SpcNames, Inst%nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'ERROR 2', RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( Inst%nSpc == 0 ) THEN
       CALL HCO_ERROR( 'No Plume species specified in HEMCO_Config.rc', RC )
       RETURN
    ENDIF

    ! Per-species scaling factors (default 1.0; can override with Scaling_<SpcName>)
    CALL GetExtSpcVal( HcoState%Config, ExtNr, Inst%nSpc, &
                       Inst%SpcNames, 'Scaling', 1.0_sp, Inst%SpcScl, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        WRITE(MSG, '(A,F6.4)') 'Fire Efficiency set as: ', Inst%Fire_Eff 
        CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
        RETURN
    ENDIF

    ! Per-species scale field names (optional)
    CALL GetExtSpcVal( HcoState%Config, ExtNr, Inst%nSpc, &
                       Inst%SpcNames, 'ScaleField', HCOX_NOSCALE, Inst%SpcScalFldNme, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'ERROR 4', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Read Fire_Efficiency from HEMCO_Config (defaults to 0.2(
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Fire_Efficiency', &
                    OptValSp=Inst%Fire_Eff, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Error reading Fire Efficiency option ', RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( .NOT. Found ) Inst%Fire_Eff = Fire_Eff_Def !If Fire_Eff not in HEMCO_Confic.rc set to default value

    ! Read Distribution_Method from HEMCO_Config, to set CO distribution method
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Distribution_Method', &
                     OptValChar=DistMethodStr, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Error reading Distribution Method option', RC, THISLOC=LOC )
       RETURN
    ENDIF

    IF ( .NOT. FOUND ) THEN
        Inst%DistMethod = 1 ! Default to uniform distribution, if no option specified
        DistMethodStr = 'Uniform'
    ELSE
        ! Set Inst%DistMethod from DistMethodStr
        IF ( TRIM(DistMethodStr) == 'Uniform' ) THEN
            Inst%DistMethod = 1
        ELSE IF ( TRIM(DistMethodStr) == 'Density' ) THEN
            Inst%DistMethod = 2
        ELSE IF ( TRIM(DistMethodStr) == 'Gaussian' ) THEN
            Inst%DistMethod = 3
        ELSE
           CALL HCO_ERROR('Unknown Distribution_Method: '//TRIM(DistMethodStr), RC )
           RETURN
        ENDIF
    ENDIF
    
    ! Activate required meteorological fields in ExtState
    ExtState%TK%DoUse       = .TRUE.    ! 3D temperature [K]
    ExtState%AIRDEN%DoUse   = .TRUE.    ! 3D dry air density [kg/m3]
    ExtState%PSC2_WET%DoUse = .TRUE.    ! 2D surface pressure [hPa]
    ExtState%T2M%DoUse      = .TRUE.    ! 2D surface temperature [K]
    ! Register the plume top height diagnostic. Created here in code (not in
    ! HEMCO_Diagn.rc) so we can set AutoFill=0, preventing HEMCO from
    ! auto-populating it with CO emission data. OutOper='Mean' disables
    ! automatic unit conversion (data in metres, not kg/m2/s).
    CALL Diagn_Create( HcoState,                                            &
                       cName     = 'ZPlume',                           &
                       long_name = 'thermodynamic_plume_injection_height',   &
                       ExtNr     = Inst%ExtNr,                              &
                       SpaceDim  = 2,                                        &
                       OutUnit   = 'm',                                      &
                       OutOper   = 'Mean',                                   &
                       AutoFill  = 0,                                        &
                       COL       = HcoState%Diagn%HcoDiagnIDDefault,         &
                       RC        = RC                                        )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot create ZPlume diagnostic', RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL Diagn_Create( HcoState, &
                       cName = 'save_QPlume', &
                       ExtNr = Inst%ExtNr, &
                       SpaceDim = 2, &
                       OutUnit = 'J', &
                       OutOper = 'Mean', &
                       AutoFill = 0, &
                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
                       RC = RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot create save_QPlume diagnostic', RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL Diagn_Create( HcoState, &
                       cName = 'save_Landtype', &
                       ExtNr = Inst%ExtNr, &
                       SpaceDim = 2, &
                       OutUnit = 'm', &
                       OutOper = 'Mean', &
                       AutoFill = 0, &
                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
                       RC = RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot create save_Landtype diagnostic', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Log extension activation
    IF ( HcoState%amIRoot ) THEN
       msg = 'Using HEMCO extension: Plume (thermodynamic plume rise)'
       IF ( HcoState%Config%doVerbose ) THEN
          CALL HCO_Msg( msg, sep1='-', LUN=HcoState%Config%hcoLogLUN )
       ELSE
          CALL HCO_Msg( msg, LUN=HcoState%Config%hcoLogLUN )
       ENDIF
       MSG = ' - Species (Name, HcoID, Scaling):'
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
       DO N = 1, Inst%nSpc
          WRITE(MSG,*) TRIM(Inst%SpcNames(N)), ', ', Inst%SpcIDs(N), ', ', Inst%SpcScl(N)
          CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
       ENDDO
       WRITE(MSG,'(A,F6.4)')  'Fire Efficiency set as: ', Inst%Fire_Eff
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
       WRITE(MSG, '(A,A)' ) 'Distribution Method set to: ', TRIM(DistMethodStr)
       CALL HCO_MSG( MSG, LUN=HcoState%Config%hcoLogLUN )
    ENDIF

    ! Cleanup
    Inst => NULL()
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_Plume_Init
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Plume_Final
!
! !DESCRIPTION: Finalizes the thermodynamic plume rise extension and removes
! the instance from the instance list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Plume_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER :: ExtState
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! HCOX_Plume_Final begins here!
    !=================================================================
    CALL InstRemove( ExtState%Plume )

  END SUBROUTINE HCOX_Plume_Final
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DistributeEmissionsDensity
!
! !DESCRIPTION: Calculates the vertical column of CO emissions to add to Hco
! flux based on the density of each level:
! CO_ij * Rho(L)/Total_Rho
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DistributeEmissionsDensity( dz, Zh, Rho, CO_ij, EmisCol )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    REAL(sp),                   INTENT(IN)          :: dz, CO_ij
    REAL(sp),                   INTENT(IN)          :: Zh(:), Rho(:)
!
! !INPUT/OUTPUT PARAMETERS: 

!
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),                   INTENT(OUT)         :: EmisCol(:)

! !REVISION HISTORY:
!  16 Jun 2026 - J. Landers    - Initial Version
!  See https://github.com/JBLanders/GEOS-Chem-Plume-Rise for version history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:!
!
    INTEGER                  :: L, NZ
    INTEGER                  :: zplm_nlay,  zplm_topL
    REAL(sp)                 :: weight_f_sum
    REAL(sp), ALLOCATABLE    :: Weight(:)

    !========================================================================
    ! DistributeEmissionsDensity begins here!
    !========================================================================
    NZ = SIZE(Zh)
    ALLOCATE( Weight(NZ) )

    Weight    = 0.0_sp
    EmisCol   = 0.0_sp
    zplm_topL = 1     ! default: surface only

    ! Find level corresponding to plume top
    DO L = 1, NZ
       IF ( Zh(L) <= dz ) THEN
          zplm_topL = L
       ELSE
          EXIT
       ENDIF
    ENDDO
    zplm_nlay = zplm_topL   ! number of layers in plume

    ! Total weight is sum of air density from surf to plume top
    weight_f_sum = SUM( Rho(1:zplm_topL) )
    DO L = 1, zplm_topL
       IF ( weight_f_sum > 0.0_sp ) THEN
          Weight(L) = Rho(L) / weight_f_sum !Weight fraction is air density at lev / tot dens
       ELSE
          Weight(L) = 1.0_sp / REAL( zplm_nlay, sp )  ! uniform fallback
       ENDIF
    ENDDO

    ! Distribute CO emissions vertically: Emis3D = CO_sfc * weight(L)
    DO L = 1, zplm_topL
       EmisCol(L) = CO_ij * Weight(L)
    ENDDO

    DEALLOCATE( Weight )

 END SUBROUTINE DistributeEmissionsDensity
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DistributeEmissionsUniform
!
! !DESCRIPTION: Calculates the vertical column of CO emissions to add to Hco
! flux based on a uniform distribution. This should be identical to using 
! Zplume scale factor from the CFFEPS files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DistributeEmissionsUniform( dz, Zh, CO_ij, EmisCol )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    REAL(sp),                   INTENT(IN)          :: dz, CO_ij
    REAL(sp),                   INTENT(IN)          :: Zh(:)
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),                   INTENT(OUT)         :: EmisCol(:)
!
! !REVISION HISTORY:
!  16 Jun 2026 - J. Landers    - Initial Version
!  See https://github.com/JBLanders/GEOS-Chem-Plume-Rise for version history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:!
!
    INTEGER                  :: L, NZ, zplm_topL

    !========================================================================
    ! DistributeEmissionsUniform begins here!
    !========================================================================
    NZ = SIZE(Zh)
    EmisCol   = 0.0_sp
    zplm_topL = 1     ! default: surface only

    ! Find level corresponding to plume top
    DO L = 1, NZ
       IF ( Zh(L) <= dz ) THEN
          zplm_topL = L
       ELSE
          EXIT
       ENDIF
    ENDDO

   !Calculate Column of Species
   DO L = 1, zplm_topL
       EmisCol(L) = CO_ij / REAL(zplm_topL, sp)
   ENDDO

 END SUBROUTINE DistributeEmissionsUniform
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DistributeEmissionsGaussian
!
! !DESCRIPTION: Calculates the vertical column of CO emissions to add to Hco
! flux based on a Gaussian distribution. The plume is centred at the midpoint
! with the plume top and bottom at the 2-sigma point
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DistributeEmissionsGaussian( dz, Zh, CO_ij, EmisCol )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    REAL(sp),                   INTENT(IN)          :: dz, CO_ij
    REAL(sp),                   INTENT(IN)          :: Zh(:)
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),                   INTENT(OUT)         :: EmisCol(:)
!
! !REVISION HISTORY:
!  16 Jun 2026 - J. Landers    - Initial Version
!  See https://github.com/JBLanders/GEOS-Chem-Plume-Rise for version history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:!
!
    INTEGER                  :: L, NZ, zplm_topL
    REAL(sp)                 :: mu, sigma
    REAL(sp)                 :: weight_f_sum
    REAL(sp), ALLOCATABLE    :: Weight(:)

    !========================================================================
    ! DistributeEmissionsGaussian begins here!
    !========================================================================
    NZ = SIZE(Zh)
    ALLOCATE (Weight(NZ) )

    EmisCol      = 0.0_sp
    Weight       = 0.0_sp
    weight_f_sum = 0.0_sp
    zplm_topL = 1     ! default: surface only

    ! Find level corresponding to plume top
    DO L = 1, NZ
       IF ( Zh(L) <= dz ) THEN
          zplm_topL = L
       ELSE
          EXIT
       ENDIF
    ENDDO

    ! Calculate Gaussian parameters
    mu    = dz / 2.0_sp
    sigma = dz / 4.0_sp
   
    ! Guard against small dz causing divide by 0 error 
    IF (sigma <= 0.0_sp) sigma = 1.0_sp

    DO L = 1, zplm_topL
        Weight(L) = EXP( -0.5_sp * ( (Zh(L) - mu) / sigma ) **2 )
    ENDDO

    weight_f_sum = SUM( Weight(1:zplm_topL) )

    DO L = 1, zplm_topL
        EmisCol(L) = CO_ij * Weight(L) / weight_f_sum
    ENDDO

    DEALLOCATE( Weight )

 END SUBROUTINE DistributeEmissionsGaussian
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Returns a pointer to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED(Inst) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Creates a new extension instance and appends it to AllInst.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)    :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT) :: Instance
    TYPE(MyInst),  POINTER       :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER :: TmpInst => NULL()
    INTEGER               :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================
    Inst => NULL()

    ! Count existing instances
    TmpInst => AllInst
    nnInst  =  0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Allocate and initialize new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Prepend to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    Instance = Inst%Instance
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Removes an instance from the instance list and deallocates it.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER :: Instance
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER               :: RC
    TYPE(MyInst), POINTER :: PrevInst
    TYPE(MyInst), POINTER :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================
    PrevInst => NULL()
    Inst     => NULL()

    CALL InstGet( Instance, Inst, RC, PrevInst=PrevInst )

    IF ( ASSOCIATED(Inst) ) THEN

       ! Deallocate standard allocatable arrays
       IF ( ALLOCATED(Inst%SpcIDs)        ) DEALLOCATE(Inst%SpcIDs)
       IF ( ALLOCATED(Inst%SpcScl)        ) DEALLOCATE(Inst%SpcScl)
       IF ( ALLOCATED(Inst%SpcNames)      ) DEALLOCATE(Inst%SpcNames)
       IF ( ALLOCATED(Inst%SpcScalFldNme) ) DEALLOCATE(Inst%SpcScalFldNme)

       ! Deallocate fire field pointers
       IF ( ASSOCIATED(Inst%Fire_CO)          ) DEALLOCATE(Inst%Fire_CO)
       IF ( ASSOCIATED(Inst%Fire_GrowthArea)  ) DEALLOCATE(Inst%Fire_GrowthArea)
       IF ( ASSOCIATED(Inst%Fire_BurnAreaTot) ) DEALLOCATE(Inst%Fire_BurnAreaTot)
       IF ( ASSOCIATED(Inst%Fire_TFC)         ) DEALLOCATE(Inst%Fire_TFC)
       IF ( ASSOCIATED(Inst%LANDTYPE)         ) DEALLOCATE(Inst%LANDTYPE)
       ! Pop instance off the linked list
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
    ENDIF

    PrevInst => NULL()
    Inst     => NULL()

  END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_Plume_Mod
