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
!   Fire_QO          - 2D fire energy      [J]         (NetCDF var: Qo)
!   Fire_BurnAreaTot - 2D burned area      [km2]       (NetCDF var: BurnAreaTot)
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
  INTEGER,  PARAMETER :: IMAX_CF   = 15000       ! max plume height and max iteration [m]
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
    REAL(hp), POINTER :: Fire_QO(:,:)          => NULL()  ! fire energy [J]
    REAL(hp), POINTER :: Fire_BurnAreaTot(:,:) => NULL()  ! burned area [km2]

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
    CHARACTER(LEN=255)    :: MSG, LOC                           ! hemco things
    INTEGER               :: I, J, L, N                         ! loop variables
    INTEGER               :: NX, NY, NZ                         ! number of grid cells
    INTEGER               :: ref1, ref2, ref3, ref4, ref5       ! Levels correspoinding to lapse rates
    INTEGER               :: ii, zplm_topL, zplm_nlay           ! plume height top and N layers
    REAL(sp)              :: Qo_ij, Area_ij, CO_ij              ! fire values at ij
    REAL(sp)              :: Le, Le1, Le2, Le3, Le4             ! different lapse rate values
    REAL(sp)              :: tsurf, ps, area                    ! met/area at each fire
    REAL(sp)              :: dz, ddz, q, Qatm, mplume, work     ! compute variables
    REAL(sp)              :: weight_f_sum                       ! total weight to distribute emissions
    REAL(sp), ALLOCATABLE :: Ph(:)                              ! pressure profile [hPa]
    REAL(sp), ALLOCATABLE :: Zh(:)                              ! altitude profile [m]
    REAL(sp), ALLOCATABLE :: Rho(:)                             ! air density profile [kg/m3]
    REAL(sp), ALLOCATABLE :: Weight(:)                          ! vertical distribution weights
    REAL(sp), ALLOCATABLE :: Emis3D(:,:,:)                      ! 3D CO emissions [kg/m2/s]
    REAL(sp), ALLOCATABLE :: ZPlume2D(:,:)                      ! save plume top height [m]
!    REAL(sp), ALLOCATABLE :: save_Le1(:,:), save_Le2(:,:), save_Le3(:,:), save_Le4(:,:)
    REAL(sp), ALLOCATABLE :: save_ref1(:,:), save_ref2(:,:), save_ref3(:,:), save_ref4(:,:), save_ref5(:,:)
    REAL(sp), ALLOCATABLE :: save_Pres(:,:,:)

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
                 Inst%Fire_QO(NX,NY),          &
                 Inst%Fire_BurnAreaTot(NX,NY) )
       Inst%Fire_CO          = 0.0_hp
       Inst%Fire_QO          = 0.0_hp
       Inst%Fire_BurnAreaTot = 0.0_hp
    ENDIF

    ! Read current 2D fire fields from HEMCO data infrastructure
    ! Use HCO_EvalFld to read Fire_CO from HcoState, and write it to our current instance variable
    CALL HCO_EvalFld( HcoState, 'Fire_CO',          Inst%Fire_CO,          RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot get Fire_CO field', RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'Fire_QO',          Inst%Fire_QO,          RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot get Fire_QO field', RC, THISLOC=LOC )
       RETURN
    ENDIF

    CALL HCO_EvalFld( HcoState, 'Fire_BurnAreaTot', Inst%Fire_BurnAreaTot, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot get Fire_BurnAreaTot field', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Allocate arrays some arrays
    ALLOCATE( Ph(NZ), Zh(NZ), Rho(NZ), Weight(NZ), Emis3D(NX,NY,NZ), ZPlume2D(NX,NY) )
!    ALLOCATE( save_Le1(NX,NY), save_Le2(NX,NY), save_Le3(NX,NY), save_Le4(NX,NY) )
    ALLOCATE( save_ref1(NX,NY), save_ref2(NX,NY), save_ref3(NX,NY), save_ref4(NX,NY), save_ref5(NX,NY) )
    ALLOCATE( save_Pres(NX,NY,NZ) )
    Emis3D   = 0.0_sp
    ZPlume2D = 0.0_sp
!    save_Le1 = 0.0_sp
!    save_Le2 = 0.0_sp
!    save_Le3 = 0.0_sp
!    save_Le4 = 0.0_sp
    save_ref1 = 0.0_sp
    save_ref2 = 0.0_sp
    save_ref3 = 0.0_sp
    save_ref4 = 0.0_sp
    save_ref4 = 0.0_sp
    save_Pres = 0.0_sp
    !=================================================================
    ! Start looping through boxes. 
    ! L=1 = surface, L=NZ = top
    !=================================================================
    DO J = 1, NY
    DO I = 1, NX
       ! Preliminary stuff setting up values before calculations

       ! Read fire parameters parameters for this box
       Qo_ij   = REAL( Inst%Fire_QO(I,J),          sp )  ! fire energy [J]
       Area_ij = REAL( Inst%Fire_BurnAreaTot(I,J),  sp )  ! burned area [km2]
       CO_ij   = REAL( Inst%Fire_CO(I,J),           sp )  ! CO flux [kg/m2/s]

       ! Skip columns with no fire emissions or fire energy ??? WRONG?? FIX???
       IF ( Qo_ij <= 0.0_sp .OR. CO_ij <= 0.0_sp ) CYCLE

       ! Build Pressure, air density and altitude array from ExtState variables
       !DO L = 1, NZ
       !   Ph(L)  = REAL( ExtState%AIRDEN%Arr%Val(I,J,L), sp ) &
       !          * RGASD_CF                                      &
       !          * REAL( ExtState%TK%Arr%Val(I,J,L),     sp )   &
       !          / 100.0_sp                   ! [Pa] -> [hPa]
       !   save_Pres(I,J,L) = Ph(L)
       !   Rho(L) = REAL( ExtState%AIRDEN%Arr%Val(I,J,L), sp )   ! [kg/m3]
       !ENDDO
       DO L = 1, NZ
          Ph(L) = ( REAL( HcoState%Grid%PEDGE%Val(I,J,L), sp ) + &     
                  REAL( HcoState%Grid%PEDGE%Val(I,J,L+1), sp ) ) / 2.0_sp
          save_Pres(I,J,L) = Ph(L)
          Rho(L) = REAL( ExtState%AIRDEN%Arr%Val(I,J,L), sp)
       ENDDO

       ! Altitude via hypsometric equation; set Z=0 at surface (L=1)
!       Zh(1) = 0.0_sp
!       DO L = 2, NZ
!          Zh(L) = Zh(L-1)                                                  &
!                + ( RGASD_CF / GRAV_CF )                                    &
!                * 0.5_sp * ( REAL(ExtState%TK%Arr%Val(I,J,L-1), sp)        &
!                           + REAL(ExtState%TK%Arr%Val(I,J,L),   sp) )      &
!                * LOG( Ph(L-1) / Ph(L) )
!       ENDDO

       ! Build Alititude using BXH Values
       Zh(1) = 0.0_sp
       DO L = 2, NZ
          Zh(L) = Zh(L-1) + REAL(HcoState%Grid%BXHEIGHT_M%Val(I,J,L-1), sp)
       ENDDO

       ! Find the vertical layer for 800, 700, 500 and 200 hPa to calculate lapse rates later
       ref1 = 1   ! surface
       ref2 = NZ  ! defaults in case threshold not found in column
       ref3 = NZ
       ref4 = NZ
       ref5 = NZ

       DO L = 1, NZ
          IF ( Ph(L) <= 85000.0_sp ) THEN; ref2 = L; EXIT; ENDIF
       ENDDO
       DO L = 1, NZ
          IF ( Ph(L) <= 70000.0_sp ) THEN; ref3 = L; EXIT; ENDIF
       ENDDO
       DO L = 1, NZ
          IF ( Ph(L) <= 50000.0_sp ) THEN; ref4 = L; EXIT; ENDIF
       ENDDO
       DO L = 1, NZ
          IF ( Ph(L) <= 20000.0_sp ) THEN; ref5 = L; EXIT; ENDIF
       ENDDO

!       IF (ref1 == ref2) THEN
!           ref2 = ref2 + 1
!       ENDIF
       save_ref1(I,J) = ref1
       save_ref2(I,J) = ref2
       save_ref3(I,J) = ref3
       save_ref4(I,J) = ref4
       save_ref5(I,J) = ref5

       ! Environmental lapse rates at reference levels
       ! Le = (T_upper - T_surface) / (Z_upper - Z_surface)
       Le1 = 0.0_sp; Le2 = 0.0_sp; Le3 = 0.0_sp; Le4 = 0.0_sp

       IF ( (Zh(ref2) - Zh(ref1)) /= 0.0_sp ) THEN
          Le1 = ( REAL(ExtState%TK%Arr%Val(I,J,ref2), sp)           &
                - REAL(ExtState%TK%Arr%Val(I,J,ref1), sp) )          &
              / ( Zh(ref2) - Zh(ref1) )   ! surface to 850 hPa
          Le2 = ( REAL(ExtState%TK%Arr%Val(I,J,ref3), sp)           &
                - REAL(ExtState%TK%Arr%Val(I,J,ref1), sp) )          &
              / ( Zh(ref3) - Zh(ref1) )   ! surface to 700 hPa
          Le3 = ( REAL(ExtState%TK%Arr%Val(I,J,ref4), sp)           &
                - REAL(ExtState%TK%Arr%Val(I,J,ref1), sp) )          &
              / ( Zh(ref4) - Zh(ref1) )   ! surface to 500 hPa
          Le4 = ( REAL(ExtState%TK%Arr%Val(I,J,ref5), sp)           &
                - REAL(ExtState%TK%Arr%Val(I,J,ref1), sp) )          &
              / ( Zh(ref5) - Zh(ref1) )   ! surface to 200 hPa
       ENDIF

!       save_Le1(I,J) = Le1
!       save_Le2(I,J) = Le2
!       save_Le3(I,J) = Le3
!       save_Le4(I,J) = Le4

       ! Preamble calculations done, now start the Anderson algorithm

       ! Initalize initial variables
       tsurf  = REAL( ExtState%TK%Arr%Val(I,J,1), sp )     ! surface-layer T [K]
       ps     = ExtState%PSC2_WET%Arr%Val(I,J) * 100.0_sp  ! hPa -> Pa
       area   = Area_ij * 1.0e6_sp                          ! km2 -> m2

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
          IF ( Qatm > Qo_ij ) THEN
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

       !-----------------------------------------------------------------
       ! Distribute Emissions
       !-----------------------------------------------------------------

       ! Find layer corresponding to plume top.
       Weight    = 0.0_sp
       zplm_topL = 1     ! default: surface only

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
          Emis3D(I,J,L) = CO_ij * Weight(L)
       ENDDO

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
    CALL Diagn_Update( HcoState, cName='ZPlume_Calc', Array2D=ZPlume2D, RC=RC )
 !   CALL Diagn_Update( HcoState, cName='Le1', Array2D=save_Le1, RC=RC )
 !   CALL Diagn_Update( HcoState, cName='Le2', Array2D=save_Le2, RC=RC )
!    CALL Diagn_Update( HcoState, cName='Le3', Array2D=save_Le3, RC=RC )
!    CALL Diagn_Update( HcoState, cName='Le4', Array2D=save_Le3, RC=RC )
    CALL Diagn_Update( HcoState, cName='ref1', Array2D=save_ref1, RC=RC )
    CALL Diagn_Update( HcoState, cName='ref2', Array2D=save_ref2, RC=RC )
    CALL Diagn_Update( HcoState, cName='ref3', Array2D=save_ref3, RC=RC )
    CALL Diagn_Update( HcoState, cName='ref4', Array2D=save_ref4, RC=RC )
    CALL Diagn_Update( HcoState, cName='ref5', Array2D=save_ref5, RC=RC )
    CALL Diagn_Update( HcoState, cName='PressureF', Array3D=save_Pres, RC=RC )
    IF (RC /= HCO_SUCCESS ) THEN
       MSG = 'Diagn_Update error: plume rise'
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF

    ! Cleanup
    DEALLOCATE( Ph, Zh, Rho, Weight, Emis3D, ZPlume2D )
!    DEALLOCATE( save_Le1, save_Le2, save_Le3, save_Le4 )
    DEALLOCATE( save_ref1, save_ref2, save_ref3, save_ref4, save_ref5 )
    DEALLOCATE( save_Pres )
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
    CHARACTER(LEN=255)    :: MSG, LOC

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
       CALL HCO_ERROR( 'ERROR 3', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Per-species scale field names (optional)
    CALL GetExtSpcVal( HcoState%Config, ExtNr, Inst%nSpc, &
                       Inst%SpcNames, 'ScaleField', HCOX_NOSCALE, Inst%SpcScalFldNme, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'ERROR 4', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Activate required meteorological fields in ExtState
    ExtState%TK%DoUse       = .TRUE.    ! 3D temperature [K]
    ExtState%AIRDEN%DoUse   = .TRUE.    ! 3D dry air density [kg/m3]
    ExtState%PSC2_WET%DoUse = .TRUE.    ! 2D surface pressure [hPa]

    ! Register the plume top height diagnostic. Created here in code (not in
    ! HEMCO_Diagn.rc) so we can set AutoFill=0, preventing HEMCO from
    ! auto-populating it with CO emission data. OutOper='Mean' disables
    ! automatic unit conversion (data in metres, not kg/m2/s).
    CALL Diagn_Create( HcoState,                                            &
                       cName     = 'ZPlume_Calc',                           &
                       long_name = 'thermodynamic_plume_injection_height',   &
                       ExtNr     = Inst%ExtNr,                              &
                       SpaceDim  = 2,                                        &
                       OutUnit   = 'm',                                      &
                       OutOper   = 'Mean',                                   &
                       AutoFill  = 0,                                        &
                       COL       = HcoState%Diagn%HcoDiagnIDDefault,         &
                       RC        = RC                                        )

!    CALL Diagn_Create( HcoState, &
!                       cName = 'Le1', &
!                       ExtNr = Inst%ExtNr, &
!                       SpaceDim = 2, &
!                       OutUnit = 'K/m', &
!                       OutOper = 'Mean', &
!                       AutoFill = 0, &
!                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
!                       RC = RC )
 !   CALL Diagn_Create( HcoState, &
 !                      cName = 'Le2', &
 !                      ExtNr = Inst%ExtNr, &
 !                      SpaceDim = 2, &
 !                      OutUnit = 'm', &
 !                      OutOper = 'Mean', &
 !                      AutoFill = 0, &
!                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
!                       RC = RC )
!    CALL Diagn_Create( HcoState, &
!                       cName = 'Le3', &
!                       ExtNr = Inst%ExtNr, &
!                       SpaceDim = 2, &
!                       OutUnit = 'm', &
!                       OutOper = 'Mean', &
!                       AutoFill = 0, &
!                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
!                       RC = RC )
!    CALL Diagn_Create( HcoState, &
!                       cName = 'Le4', &
!                       ExtNr = Inst%ExtNr, &
!                       SpaceDim = 2, &
!                       OutUnit = 'm', &
!                       OutOper = 'Mean', &
!                       AutoFill = 0, &
!                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
!                       RC = RC )

    CALL Diagn_Create( HcoState, &
                       cName = 'ref1', &
                       ExtNr = Inst%ExtNr, &
                       SpaceDim = 2, &
                       OutUnit = 'm', &
                       OutOper = 'Mean', &
                       AutoFill = 0, &
                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
                       RC = RC )
    CALL Diagn_Create( HcoState, &
                       cName = 'ref2', &
                       ExtNr = Inst%ExtNr, &
                       SpaceDim = 2, &
                       OutUnit = 'm', &
                       OutOper = 'Mean', &
                       AutoFill = 0, &
                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
                       RC = RC )
    CALL Diagn_Create( HcoState, &
                       cName = 'ref3', &
                       ExtNr = Inst%ExtNr, &
                       SpaceDim = 2, &
                       OutUnit = 'm', &
                       OutOper = 'Mean', &
                       AutoFill = 0, &
                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
                       RC = RC )
    CALL Diagn_Create( HcoState, &
                       cName = 'ref4', &
                       ExtNr = Inst%ExtNr, &
                       SpaceDim = 2, &
                       OutUnit = 'm', &
                       OutOper = 'Mean', &
                       AutoFill = 0, &
                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
                       RC = RC )
    CALL Diagn_Create( HcoState, &
                       cName = 'ref5', &
                       ExtNr = Inst%ExtNr, &
                       SpaceDim = 2, &
                       OutUnit = 'm', &
                       OutOper = 'Mean', &
                       AutoFill = 0, &
                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
                       RC = RC )
    CALL Diagn_Create( HcoState, &
                       cName = 'PressureF', &
                       ExtNr = Inst%ExtNr, &
                       SpaceDim =3, &
                       OutUnit = 'Pa', &
                       OutOper = 'Mean', &
                       AutoFill = 0, &
                       COL = HcoState%Diagn%HcoDiagnIDDefault, &
                       RC = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot create ZPlume_Calc diagnostic', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Log extension activation
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Cannot create ZPlume_Calc diagnostic', RC, THISLOC=LOC )
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
       IF ( ASSOCIATED(Inst%Fire_QO)          ) DEALLOCATE(Inst%Fire_QO)
       IF ( ASSOCIATED(Inst%Fire_BurnAreaTot) ) DEALLOCATE(Inst%Fire_BurnAreaTot)

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
