! Harsh Kamath (NOAA GSL/CIRES), September 2026


module WaterVarInTransferMod

  use Machine, only: kind_noahmp
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpVarType, only: noahmp_type

  implicit none

contains

!=== initialize with input data or table values

  subroutine WaterVarInTransfer(noahmp, NoahmpIO, N)

    implicit none

    type(noahmp_type),        intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO
	integer,                  intent(in)    :: N

    ! local variables 
    integer                            :: IndexSoilLayer
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilSand
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilClay
    real(kind=kind_noahmp), allocatable, dimension(:) :: SoilOrg

! -------------------------------------------------------------------------
    associate(                                                         &
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer    ,&
              VegType         => noahmp%config%domain%VegType         ,&
              SoilType        => noahmp%config%domain%SoilType        ,&
              FlagUrban       => noahmp%config%domain%FlagUrban       ,&
              RunoffSlopeType => noahmp%config%domain%RunoffSlopeType ,&
              NumSnowLayerNeg => noahmp%config%domain%NumSnowLayerNeg  &
             )
! -------------------------------------------------------------------------

    ! water state variables
    noahmp%water%state%CanopyLiqWater                     = NoahmpIO%CANLIQXY   (N)
    noahmp%water%state%CanopyIce                          = NoahmpIO%CANICEXY   (N)
    noahmp%water%state%CanopyWetFrac                      = NoahmpIO%FWETXY     (N)
    noahmp%water%state%SnowWaterEquiv                     = NoahmpIO%SNOW       (N)
    noahmp%water%state%SnowWaterEquivPrev                 = NoahmpIO%SNEQVOXY   (N) 
    noahmp%water%state%SnowDepth                          = NoahmpIO%SNOWH      (N)
    noahmp%water%state%IrrigationFracFlood                = NoahmpIO%FIFRACT    (N)
    noahmp%water%state%IrrigationAmtFlood                 = NoahmpIO%IRWATFI    (N)
    noahmp%water%state%IrrigationFracMicro                = NoahmpIO%MIFRACT    (N)
    noahmp%water%state%IrrigationAmtMicro                 = NoahmpIO%IRWATMI    (N) 
    noahmp%water%state%IrrigationFracSprinkler            = NoahmpIO%SIFRACT    (N)
    noahmp%water%state%IrrigationAmtSprinkler             = NoahmpIO%IRWATSI    (N)  
    noahmp%water%state%WaterTableDepth                    = NoahmpIO%ZWTXY      (N) 
    noahmp%water%state%SoilMoistureToWT                   = NoahmpIO%SMCWTDXY   (N)
    noahmp%water%state%TileDrainFrac                      = NoahmpIO%TD_FRACTION(N)
    noahmp%water%state%WaterStorageAquifer                = NoahmpIO%WAXY       (N)
    noahmp%water%state%WaterStorageSoilAqf                = NoahmpIO%WTXY       (N)
    noahmp%water%state%WaterStorageLake                   = NoahmpIO%WSLAKEXY   (N)
    noahmp%water%state%IrrigationFracGrid                 = NoahmpIO%IRFRACT    (N)
    noahmp%water%state%IrrigationCntSprinkler             = NoahmpIO%IRNUMSI    (N)     
    noahmp%water%state%IrrigationCntMicro                 = NoahmpIO%IRNUMMI    (N)
    noahmp%water%state%IrrigationCntFlood                 = NoahmpIO%IRNUMFI    (N)
    noahmp%water%state%SnowIce     (-NumSnowLayerMax+1:0) = NoahmpIO%SNICEXY    (-NumSnowLayerMax+1:0,N)
    noahmp%water%state%SnowLiqWater(-NumSnowLayerMax+1:0) = NoahmpIO%SNLIQXY    (-NumSnowLayerMax+1:0,N)
    noahmp%water%state%SoilLiqWater      (1:NumSoilLayer) = NoahmpIO%SH2O       (1:NumSoilLayer,N)
    noahmp%water%state%SoilMoisture      (1:NumSoilLayer) = NoahmpIO%SMOIS      (1:NumSoilLayer,N)    
    noahmp%water%state%SoilMoistureEqui  (1:NumSoilLayer) = NoahmpIO%SMOISEQ    (1:NumSoilLayer,N)
    noahmp%water%state%RechargeGwDeepWT                   = 0.0
    noahmp%water%state%RechargeGwShallowWT                = 0.0
    if ( noahmp%config%nmlist%OptWetlandModel > 0 ) then
       noahmp%water%state%SoilSaturateFrac                = NoahmpIO%FSATXY     (N)
       noahmp%water%state%WaterStorageWetland             = NoahmpIO%WSURFXY    (N)
    else
       ! The modular balance checker includes wetland storage for every land
       ! point.  Do not leave its initialization sentinel in that sum when
       ! the wetland option is disabled in CCPP-SCM.
       noahmp%water%state%SoilSaturateFrac                = 0.0_kind_noahmp
       noahmp%water%state%WaterStorageWetland             = 0.0_kind_noahmp
    endif
    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       noahmp%water%state%SnowRadius  (-NumSnowLayerMax+1:0)       = NoahmpIO%SNRDSXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassBChydrophi(-NumSnowLayerMax+1:0)     = NoahmpIO%BCPHIXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassBChydropho(-NumSnowLayerMax+1:0)     = NoahmpIO%BCPHOXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassOChydrophi(-NumSnowLayerMax+1:0)     = NoahmpIO%OCPHIXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassOChydropho(-NumSnowLayerMax+1:0)     = NoahmpIO%OCPHOXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassDust1(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST1XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassDust2(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST2XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassDust3(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST3XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassDust4(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST4XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassDust5(-NumSnowLayerMax+1:0)          = NoahmpIO%DUST5XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcBChydrophi(-NumSnowLayerMax+1:0) = NoahmpIO%MassConcBCPHIXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcBChydropho(-NumSnowLayerMax+1:0) = NoahmpIO%MassConcBCPHOXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcOChydrophi(-NumSnowLayerMax+1:0) = NoahmpIO%MassConcOCPHIXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcOChydropho(-NumSnowLayerMax+1:0) = NoahmpIO%MassConcOCPHOXY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcDust1(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST1XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcDust2(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST2XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcDust3(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST3XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcDust4(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST4XY (-NumSnowLayerMax+1:0,N)
       noahmp%water%state%MassConcDust5(-NumSnowLayerMax+1:0)      = NoahmpIO%MassConcDUST5XY (-NumSnowLayerMax+1:0,N)
    endif


    ! water flux variables
    noahmp%water%flux%EvapSoilSfcLiqAcc                   = NoahmpIO%ACC_QSEVAXY (N)
    noahmp%water%flux%SoilSfcInflowAcc                    = NoahmpIO%ACC_QINSURXY(N)
    noahmp%water%flux%SfcWaterTotChgAcc                   = NoahmpIO%ACC_DWATERXY(N)
    noahmp%water%flux%PrecipTotAcc                        = NoahmpIO%ACC_PRCPXY  (N)
    noahmp%water%flux%EvapCanopyNetAcc                    = NoahmpIO%ACC_ECANXY  (N)
    noahmp%water%flux%TranspirationAcc                    = NoahmpIO%ACC_ETRANXY (N)
    noahmp%water%flux%EvapGroundNetAcc                    = NoahmpIO%ACC_EDIRXY  (N)
    noahmp%water%flux%TranspWatLossSoilAcc(1:NumSoilLayer)= NoahmpIO%ACC_ETRANIXY(1:NumSoilLayer,N)
    noahmp%water%flux%GlacierExcessFlowAcc                = NoahmpIO%ACC_GLAFLWXY(N)
    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       noahmp%water%flux%SnowFreezeRate(-NumSnowLayerMax+1:0) = NoahmpIO%SNFRXY(-NumSnowLayerMax+1:0,N)
    endif


    ! water parameter variables
    noahmp%water%param%DrainSoilLayerInd                  = NoahmpIO%DRAIN_LAYER_OPT_TABLE
    noahmp%water%param%CanopyLiqHoldCap                   = NoahmpIO%CH2OP_TABLE(VegType)
    noahmp%water%param%SnowCompactBurdenFac               = NoahmpIO%C2_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactAgingFac1               = NoahmpIO%C3_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactAgingFac2               = NoahmpIO%C4_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactAgingFac3               = NoahmpIO%C5_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactAgingMax                = NoahmpIO%DM_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowViscosityCoeff                 = NoahmpIO%ETA0_SNOWCOMPACT_TABLE
    noahmp%water%param%SnowCompactmAR24                   = NoahmpIO%SNOWCOMPACTm_AR24_TABLE
    noahmp%water%param%SnowCompactbAR24                   = NoahmpIO%SNOWCOMPACTb_AR24_TABLE
    noahmp%water%param%SnowCompactP1AR24                  = NoahmpIO%SNOWCOMPACT_P1_AR24_TABLE
    noahmp%water%param%SnowCompactP2AR24                  = NoahmpIO%SNOWCOMPACT_P2_AR24_TABLE
    noahmp%water%param%SnowCompactP3AR24                  = NoahmpIO%SNOWCOMPACT_P3_AR24_TABLE
    noahmp%water%param%SnowCoverM1AR25                    = NoahmpIO%SCFm1_AR25_TABLE
    noahmp%water%param%SnowCoverM2AR25                    = NoahmpIO%SCFm2_AR25_TABLE
    noahmp%water%param%SnowCoverFac1AR25                  = NoahmpIO%SCfac1_AR25_TABLE
    noahmp%water%param%SnowCoverFac2AR25                  = NoahmpIO%SCfac2_AR25_TABLE
    noahmp%water%param%BurdenFacUpAR24                    = NoahmpIO%SNOWCOMPACT_Up_AR24_TABLE
    noahmp%water%param%SnowLiqFracMax                     = NoahmpIO%SNLIQMAXFRAC_TABLE
    noahmp%water%param%SnowLiqHoldCap                     = NoahmpIO%SSI_TABLE
    noahmp%water%param%SnowLiqReleaseFac                  = NoahmpIO%SNOW_RET_FAC_TABLE
    noahmp%water%param%IrriFloodRateFac                   = NoahmpIO%FIRTFAC_TABLE
    noahmp%water%param%IrriMicroRate                      = NoahmpIO%MICIR_RATE_TABLE
    noahmp%water%param%SoilConductivityRef                = NoahmpIO%REFDK_TABLE
    noahmp%water%param%SoilInfilFacRef                    = NoahmpIO%REFKDT_TABLE
    noahmp%water%param%GroundFrzCoeff                     = NoahmpIO%FRZK_TABLE
    noahmp%water%param%GridTopoIndex                      = NoahmpIO%TIMEAN_TABLE
    noahmp%water%param%SoilSfcSatFracMax                  = NoahmpIO%FSATMX_TABLE
    noahmp%water%param%SpecYieldGw                        = NoahmpIO%ROUS_TABLE
    noahmp%water%param%MicroPoreContent                   = NoahmpIO%CMIC_TABLE
    noahmp%water%param%WaterStorageLakeMax                = NoahmpIO%WSLMAX_TABLE
    noahmp%water%param%SnoWatEqvMaxGlacier                = NoahmpIO%SWEMAXGLA_TABLE
    noahmp%water%param%IrriStopDayBfHarvest               = NoahmpIO%IRR_HAR_TABLE
    noahmp%water%param%IrriTriggerLaiMin                  = NoahmpIO%IRR_LAI_TABLE
    noahmp%water%param%SoilWatDeficitAllow                = NoahmpIO%IRR_MAD_TABLE
    noahmp%water%param%IrriFloodLossFrac                  = NoahmpIO%FILOSS_TABLE
    noahmp%water%param%IrriSprinklerRate                  = NoahmpIO%SPRIR_RATE_TABLE
    noahmp%water%param%IrriFracThreshold                  = NoahmpIO%IRR_FRAC_TABLE
    noahmp%water%param%IrriStopPrecipThr                  = NoahmpIO%IR_RAIN_TABLE
    noahmp%water%param%SnowfallDensityMax                 = NoahmpIO%SNOWDEN_MAX_TABLE
    noahmp%water%param%SnowMassFullCoverOld               = NoahmpIO%SWEMX_TABLE
    noahmp%water%param%SoilMatPotentialWilt               = NoahmpIO%PSIWLT_TABLE
    noahmp%water%param%SnowMeltFac                        = NoahmpIO%MFSNO_TABLE(VegType)
    noahmp%water%param%SnowCoverFac                       = NoahmpIO%SCFFAC_TABLE(VegType)
    noahmp%water%param%InfilFacVic                        = NoahmpIO%BVIC_TABLE(SoilType(1))
    noahmp%water%param%TensionWatDistrInfl                = NoahmpIO%AXAJ_TABLE(SoilType(1))
    noahmp%water%param%TensionWatDistrShp                 = NoahmpIO%BXAJ_TABLE(SoilType(1))
    noahmp%water%param%FreeWatDistrShp                    = NoahmpIO%XXAJ_TABLE(SoilType(1))
    noahmp%water%param%InfilHeteroDynVic                  = NoahmpIO%BBVIC_TABLE(SoilType(1))
    noahmp%water%param%InfilCapillaryDynVic               = NoahmpIO%GDVIC_TABLE(SoilType(1))
    noahmp%water%param%InfilFacDynVic                     = NoahmpIO%BDVIC_TABLE(SoilType(1))
    noahmp%water%param%TileDrainCoeffSp                   = NoahmpIO%TD_DC_TABLE(SoilType(1))
    noahmp%water%param%TileDrainTubeDepth                 = NoahmpIO%TD_DEPTH_TABLE(SoilType(1))
    noahmp%water%param%DrainFacSoilWat                    = NoahmpIO%TDSMC_FAC_TABLE(SoilType(1))
    noahmp%water%param%TileDrainCoeff                     = NoahmpIO%TD_DCOEF_TABLE(SoilType(1))
    noahmp%water%param%DrainDepthToImperv                 = NoahmpIO%TD_ADEPTH_TABLE(SoilType(1))
    noahmp%water%param%LateralWatCondFac                  = NoahmpIO%KLAT_FAC_TABLE(SoilType(1))
    noahmp%water%param%TileDrainDepth                     = NoahmpIO%TD_DDRAIN_TABLE(SoilType(1))
    noahmp%water%param%DrainTubeDist                      = NoahmpIO%TD_SPAC_TABLE(SoilType(1))
    noahmp%water%param%DrainTubeRadius                    = NoahmpIO%TD_RADI_TABLE(SoilType(1))
    noahmp%water%param%DrainWatDepToImperv                = NoahmpIO%TD_D_TABLE(SoilType(1))
    noahmp%water%param%NumSoilLayerRoot                   = NoahmpIO%NROOT_TABLE(VegType)
    noahmp%water%param%SoilDrainSlope                     = NoahmpIO%SLOPE_TABLE(RunoffSlopeType)
    noahmp%water%param%WetlandCapMax                      = NoahmpIO%WCAP_TABLE

    ! SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 )then
       noahmp%water%param%snowage_tau                     = NoahmpIO%snowage_tau
       noahmp%water%param%snowage_kappa                   = NoahmpIO%snowage_kappa
       noahmp%water%param%snowage_drdt0                   = NoahmpIO%snowage_drdt0
       noahmp%water%param%SnowRadiusMin                   = NoahmpIO%SnowRadiusMin_TABLE
       noahmp%water%param%FreshSnowRadiusMax              = NoahmpIO%FreshSnowRadiusMax_TABLE
       noahmp%water%param%SnowRadiusRefrz                 = NoahmpIO%SnowRadiusRefrz_TABLE
       noahmp%water%param%ScavEffMeltScale                = NoahmpIO%ScavEffMeltScale_TABLE
       noahmp%water%param%ScavEffMeltBCphi                = NoahmpIO%ScavEffMeltBCphi_TABLE
       noahmp%water%param%ScavEffMeltBCpho                = NoahmpIO%ScavEffMeltBCpho_TABLE
       noahmp%water%param%ScavEffMeltOCphi                = NoahmpIO%ScavEffMeltOCphi_TABLE
       noahmp%water%param%ScavEffMeltOCpho                = NoahmpIO%ScavEffMeltOCpho_TABLE
       noahmp%water%param%ScavEffMeltDust1                = NoahmpIO%ScavEffMeltDust1_TABLE
       noahmp%water%param%ScavEffMeltDust2                = NoahmpIO%ScavEffMeltDust2_TABLE
       noahmp%water%param%ScavEffMeltDust3                = NoahmpIO%ScavEffMeltDust3_TABLE
       noahmp%water%param%ScavEffMeltDust4                = NoahmpIO%ScavEffMeltDust4_TABLE
       noahmp%water%param%ScavEffMeltDust5                = NoahmpIO%ScavEffMeltDust5_TABLE
       noahmp%water%param%SnowRadiusMax                   = NoahmpIO%SnowRadiusMax_TABLE
       noahmp%water%param%SnowWetAgeC1Brun89              = NoahmpIO%SnowWetAgeC1Brun89_TABLE
       noahmp%water%param%SnowWetAgeC2Brun89              = NoahmpIO%SnowWetAgeC2Brun89_TABLE
       noahmp%water%param%SnowAgeScaleFac                 = NoahmpIO%SnowAgeScaleFac_TABLE
    endif

    ! soil properties
    do IndexSoilLayer = 1, size(SoilType)
       noahmp%water%param%SoilMoistureSat       (IndexSoilLayer) = NoahmpIO%SMCMAX_TABLE(SoilType(IndexSoilLayer))
       noahmp%water%param%SoilMoistureWilt      (IndexSoilLayer) = NoahmpIO%SMCWLT_TABLE(SoilType(IndexSoilLayer))
       noahmp%water%param%SoilMoistureFieldCap  (IndexSoilLayer) = NoahmpIO%SMCREF_TABLE(SoilType(IndexSoilLayer))
       noahmp%water%param%SoilMoistureDry       (IndexSoilLayer) = NoahmpIO%SMCDRY_TABLE(SoilType(IndexSoilLayer))
       noahmp%water%param%SoilWatDiffusivitySat (IndexSoilLayer) = NoahmpIO%DWSAT_TABLE (SoilType(IndexSoilLayer))
       noahmp%water%param%SoilWatConductivitySat(IndexSoilLayer) = NoahmpIO%DKSAT_TABLE (SoilType(IndexSoilLayer))
       noahmp%water%param%SoilExpCoeffB         (IndexSoilLayer) = NoahmpIO%BEXP_TABLE  (SoilType(IndexSoilLayer))
       noahmp%water%param%SoilMatPotentialSat   (IndexSoilLayer) = NoahmpIO%PSISAT_TABLE(SoilType(IndexSoilLayer))
    enddo
   
    ! spatial varying soil texture and properties directly from input
    if ( noahmp%config%nmlist%OptSoilProperty == 4 ) then
       ! 3D soil properties
       noahmp%water%param%SoilExpCoeffB          = NoahmpIO%BEXP_3D  (1:NumSoilLayer,N) ! C-H B exponent
       noahmp%water%param%SoilMoistureDry        = NoahmpIO%SMCDRY_3D(1:NumSoilLayer,N) ! Soil Moisture Limit: Dry
       noahmp%water%param%SoilMoistureWilt       = NoahmpIO%SMCWLT_3D(1:NumSoilLayer,N) ! Soil Moisture Limit: Wilt
       noahmp%water%param%SoilMoistureFieldCap   = NoahmpIO%SMCREF_3D(1:NumSoilLayer,N) ! Soil Moisture Limit: Reference
       noahmp%water%param%SoilMoistureSat        = NoahmpIO%SMCMAX_3D(1:NumSoilLayer,N) ! Soil Moisture Limit: Max
       noahmp%water%param%SoilWatConductivitySat = NoahmpIO%DKSAT_3D (1:NumSoilLayer,N) ! Saturated Soil Conductivity
       noahmp%water%param%SoilWatDiffusivitySat  = NoahmpIO%DWSAT_3D (1:NumSoilLayer,N) ! Saturated Soil Diffusivity
       noahmp%water%param%SoilMatPotentialSat    = NoahmpIO%PSISAT_3D(1:NumSoilLayer,N) ! Saturated Matric Potential
       noahmp%water%param%SoilConductivityRef    = NoahmpIO%REFDK_2D (N)                ! Reference Soil Conductivity
       noahmp%water%param%SoilInfilFacRef        = NoahmpIO%REFKDT_2D(N)                ! Soil Infiltration Parameter
       ! 2D additional runoff6~8 parameters
       noahmp%water%param%InfilFacVic            = NoahmpIO%BVIC_2D (N)                 ! VIC model infiltration parameter
       noahmp%water%param%TensionWatDistrInfl    = NoahmpIO%AXAJ_2D (N)                 ! Xinanjiang: Tension water distribution inflection parameter
       noahmp%water%param%TensionWatDistrShp     = NoahmpIO%BXAJ_2D (N)                 ! Xinanjiang: Tension water distribution shape parameter
       noahmp%water%param%FreeWatDistrShp        = NoahmpIO%XXAJ_2D (N)                 ! Xinanjiang: Free water distribution shape parameter
       noahmp%water%param%InfilFacDynVic         = NoahmpIO%BDVIC_2D(N)                 ! VIC model infiltration parameter
       noahmp%water%param%InfilCapillaryDynVic   = NoahmpIO%GDVIC_2D(N)                 ! Mean Capillary Drive for infiltration models
       noahmp%water%param%InfilHeteroDynVic      = NoahmpIO%BBVIC_2D(N)                 ! DVIC heterogeniety parameter for infiltraton
       ! 2D irrigation params
       noahmp%water%param%IrriFracThreshold      = NoahmpIO%IRR_FRAC_2D  (N)            ! irrigation Fraction
       noahmp%water%param%IrriStopDayBfHarvest   = NoahmpIO%IRR_HAR_2D   (N)            ! number of days before harvest date to stop irrigation 
       noahmp%water%param%IrriTriggerLaiMin      = NoahmpIO%IRR_LAI_2D   (N)            ! Minimum lai to trigger irrigation
       noahmp%water%param%SoilWatDeficitAllow    = NoahmpIO%IRR_MAD_2D   (N)            ! management allowable deficit (0-1)
       noahmp%water%param%IrriFloodLossFrac      = NoahmpIO%FILOSS_2D    (N)            ! fraction of flood irrigation loss (0-1) 
       noahmp%water%param%IrriSprinklerRate      = NoahmpIO%SPRIR_RATE_2D(N)            ! mm/h, sprinkler irrigation rate
       noahmp%water%param%IrriMicroRate          = NoahmpIO%MICIR_RATE_2D(N)            ! mm/h, micro irrigation rate
       noahmp%water%param%IrriFloodRateFac       = NoahmpIO%FIRTFAC_2D   (N)            ! flood application rate factor
       noahmp%water%param%IrriStopPrecipThr      = NoahmpIO%IR_RAIN_2D   (N)            ! maximum precipitation to stop irrigation trigger
       ! 2D tile drainage parameters
       noahmp%water%param%LateralWatCondFac      = NoahmpIO%KLAT_FAC (N)                ! factor multiplier to hydraulic conductivity
       noahmp%water%param%DrainFacSoilWat        = NoahmpIO%TDSMC_FAC(N)                ! factor multiplier to field capacity
       noahmp%water%param%TileDrainCoeffSp       = NoahmpIO%TD_DC    (N)                ! drainage coefficient for simple
       noahmp%water%param%TileDrainCoeff         = NoahmpIO%TD_DCOEF (N)                ! drainge coefficient for Hooghoudt 
       noahmp%water%param%TileDrainDepth         = NoahmpIO%TD_DDRAIN(N)                ! depth of drain
       noahmp%water%param%DrainTubeRadius        = NoahmpIO%TD_RADI  (N)                ! tile tube radius
       noahmp%water%param%DrainTubeDist          = NoahmpIO%TD_SPAC  (N)                ! tile spacing
    endif

    ! spatial varying wetland parameters from input
    if ( noahmp%config%nmlist%OptWetlandModel == 2 ) then
       noahmp%water%param%SoilSfcSatFracMax      = NoahmpIO%FSATMX(N)
       noahmp%water%param%WetlandCapMax          = NoahmpIO%WCAP(N)
    endif 

    ! derived water parameters
    noahmp%water%param%SoilInfilMaxCoeff  = noahmp%water%param%SoilInfilFacRef *           &
                                            noahmp%water%param%SoilWatConductivitySat(1) / &
                                            noahmp%water%param%SoilConductivityRef
    if ( FlagUrban .eqv. .true. ) then
       noahmp%water%param%SoilMoistureSat      = 0.45
       noahmp%water%param%SoilMoistureFieldCap = 0.42
       noahmp%water%param%SoilMoistureWilt     = 0.40
       noahmp%water%param%SoilMoistureDry      = 0.40
    endif

    if ( SoilType(1) /= 14 ) then
       noahmp%water%param%SoilImpervFracCoeff = noahmp%water%param%GroundFrzCoeff *       &
                                                ((noahmp%water%param%SoilMoistureSat(1) / &
                                                 noahmp%water%param%SoilMoistureFieldCap(1)) * (0.412/0.468))
    endif

    noahmp%water%state%SnowIceFracPrev = 0.0
    noahmp%water%state%SnowIceFracPrev(NumSnowLayerNeg+1:0) = NoahmpIO%SNICEXY(NumSnowLayerNeg+1:0,N) /  & 
                                                              (NoahmpIO%SNICEXY(NumSnowLayerNeg+1:0,N) + &
                                                               NoahmpIO%SNLIQXY(NumSnowLayerNeg+1:0,N))

    if ( (noahmp%config%nmlist%OptSoilProperty == 3) .and. (.not. noahmp%config%domain%FlagUrban) ) then
		stop 'WaterVarInTransfer'
       ! if (.not. allocated(SoilSand)) allocate( SoilSand(1:NumSoilLayer) )
       ! if (.not. allocated(SoilClay)) allocate( SoilClay(1:NumSoilLayer) )
       ! if (.not. allocated(SoilOrg) ) allocate( SoilOrg (1:NumSoilLayer) )
       ! SoilSand = 0.01 * NoahmpIO%soilcomp(1:NumSoilLayer,N)
       ! SoilClay = 0.01 * NoahmpIO%soilcomp((NumSoilLayer+1):(NumSoilLayer*2),N)
       ! SoilOrg  = 0.0
       ! if (noahmp%config%nmlist%OptPedotransfer == 1) &
          ! call PedoTransferSR2006(block,noahmp,SoilSand,SoilClay,SoilOrg)
       ! deallocate(SoilSand)
       ! deallocate(SoilClay)
       ! deallocate(SoilOrg )
    endif

    end associate

  end subroutine WaterVarInTransfer

end module WaterVarInTransferMod
