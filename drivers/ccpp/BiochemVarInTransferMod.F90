! Harsh Kamath (NOAA GSL/CIRES), September 2026

module BiochemVarInTransferMod

  use Machine, only: kind_noahmp
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpVarType, only: noahmp_type

  implicit none

contains

!=== initialize with input data or table values

  subroutine BiochemVarInTransfer(noahmp, NoahmpIO, N)

    implicit none

    type(noahmp_type),        intent(inout) :: noahmp
    type(NoahmpIO_type),      intent(inout) :: NoahmpIO
	integer,                  intent(in)    :: N

! -------------------------------------------------------------------------
    associate(                                                   &
              VegType      => noahmp%config%domain%VegType      ,&
              CropType     => noahmp%config%domain%CropType     ,&
              OptCropModel => noahmp%config%nmlist%OptCropModel  &
             )
! -------------------------------------------------------------------------

    ! biochem state variables
    noahmp%biochem%state%PlantGrowStage             = NoahmpIO%PGSXY   (N)   
    noahmp%biochem%state%LeafMass                   = NoahmpIO%LFMASSXY(N)
    noahmp%biochem%state%RootMass                   = NoahmpIO%RTMASSXY(N)
    noahmp%biochem%state%StemMass                   = NoahmpIO%STMASSXY(N) 
    noahmp%biochem%state%WoodMass                   = NoahmpIO%WOODXY  (N) 
    noahmp%biochem%state%CarbonMassDeepSoil         = NoahmpIO%STBLCPXY(N) 
    noahmp%biochem%state%CarbonMassShallowSoil      = NoahmpIO%FASTCPXY(N)
    noahmp%biochem%state%GrainMass                  = NoahmpIO%GRAINXY (N)  
    noahmp%biochem%state%GrowDegreeDay              = NoahmpIO%GDDXY   (N)  
    noahmp%biochem%state%NitrogenConcFoliage        = 1.0  ! for now, set to nitrogen saturation

    ! biochem parameter variables
    noahmp%biochem%param%NitrogenConcFoliageMax     = NoahmpIO%FOLNMX_TABLE (VegType)
    noahmp%biochem%param%QuantumEfficiency25C       = NoahmpIO%QE25_TABLE   (VegType)
    noahmp%biochem%param%CarboxylRateMax25C         = NoahmpIO%VCMX25_TABLE (VegType)
    noahmp%biochem%param%CarboxylRateMaxQ10         = NoahmpIO%AVCMX_TABLE  (VegType)
    noahmp%biochem%param%PhotosynPathC3             = NoahmpIO%C3PSN_TABLE  (VegType)
    noahmp%biochem%param%SlopeConductToPhotosyn     = NoahmpIO%MP_TABLE     (VegType)
    noahmp%biochem%param%RespMaintQ10               = NoahmpIO%ARM_TABLE    (VegType)
    noahmp%biochem%param%RespMaintLeaf25C           = NoahmpIO%RMF25_TABLE  (VegType)
    noahmp%biochem%param%RespMaintStem25C           = NoahmpIO%RMS25_TABLE  (VegType)
    noahmp%biochem%param%RespMaintRoot25C           = NoahmpIO%RMR25_TABLE  (VegType)
    noahmp%biochem%param%WoodToRootRatio            = NoahmpIO%WRRAT_TABLE  (VegType)
    noahmp%biochem%param%WoodPoolIndex              = NoahmpIO%WDPOOL_TABLE (VegType)
    noahmp%biochem%param%TurnoverCoeffLeafVeg       = NoahmpIO%LTOVRC_TABLE (VegType)
    noahmp%biochem%param%TemperaureLeafFreeze       = NoahmpIO%TDLEF_TABLE  (VegType)
    noahmp%biochem%param%LeafDeathWaterCoeffVeg     = NoahmpIO%DILEFW_TABLE (VegType)
    noahmp%biochem%param%LeafDeathTempCoeffVeg      = NoahmpIO%DILEFC_TABLE (VegType)
    noahmp%biochem%param%GrowthRespFrac             = NoahmpIO%FRAGR_TABLE  (VegType)
    noahmp%biochem%param%MicroRespCoeff             = NoahmpIO%MRP_TABLE    (VegType)
    noahmp%biochem%param%TemperatureMinPhotosyn     = NoahmpIO%TMIN_TABLE   (VegType)
    noahmp%biochem%param%LeafAreaPerMass1side       = NoahmpIO%SLA_TABLE    (VegType)
    noahmp%biochem%param%StemAreaIndexMin           = NoahmpIO%XSAMIN_TABLE (VegType)
    noahmp%biochem%param%WoodAllocFac               = NoahmpIO%BF_TABLE     (VegType)
    noahmp%biochem%param%WaterStressCoeff           = NoahmpIO%WSTRC_TABLE  (VegType)
    noahmp%biochem%param%LeafAreaIndexMin           = NoahmpIO%LAIMIN_TABLE (VegType)
    noahmp%biochem%param%TurnoverCoeffRootVeg       = NoahmpIO%RTOVRC_TABLE (VegType)
    noahmp%biochem%param%WoodRespCoeff              = NoahmpIO%RSWOODC_TABLE(VegType)
    ! crop model specific parameters
    if ( (OptCropModel > 0) .and. (CropType > 0) ) then
       noahmp%biochem%param%DatePlanting            = NoahmpIO%PLTDAY_TABLE   (CropType)
       noahmp%biochem%param%DateHarvest             = NoahmpIO%HSDAY_TABLE    (CropType)
       noahmp%biochem%param%NitrogenConcFoliageMax  = NoahmpIO%FOLNMXI_TABLE  (CropType)
       noahmp%biochem%param%QuantumEfficiency25C    = NoahmpIO%QE25I_TABLE    (CropType)
       noahmp%biochem%param%CarboxylRateMax25C      = NoahmpIO%VCMX25I_TABLE  (CropType)
       noahmp%biochem%param%CarboxylRateMaxQ10      = NoahmpIO%AVCMXI_TABLE   (CropType)
       noahmp%biochem%param%PhotosynPathC3          = NoahmpIO%C3PSNI_TABLE   (CropType)
       noahmp%biochem%param%SlopeConductToPhotosyn  = NoahmpIO%MPI_TABLE      (CropType)
       noahmp%biochem%param%RespMaintQ10            = NoahmpIO%Q10MR_TABLE    (CropType)
       noahmp%biochem%param%RespMaintLeaf25C        = NoahmpIO%LFMR25_TABLE   (CropType)
       noahmp%biochem%param%RespMaintStem25C        = NoahmpIO%STMR25_TABLE   (CropType)
       noahmp%biochem%param%RespMaintRoot25C        = NoahmpIO%RTMR25_TABLE   (CropType)
       noahmp%biochem%param%GrowthRespFrac          = NoahmpIO%FRA_GR_TABLE   (CropType)
       noahmp%biochem%param%TemperaureLeafFreeze    = NoahmpIO%LEFREEZ_TABLE  (CropType)
       noahmp%biochem%param%LeafAreaPerBiomass      = NoahmpIO%BIO2LAI_TABLE  (CropType)
       noahmp%biochem%param%TempBaseGrowDegDay      = NoahmpIO%GDDTBASE_TABLE (CropType)
       noahmp%biochem%param%TempMaxGrowDegDay       = NoahmpIO%GDDTCUT_TABLE  (CropType)
       noahmp%biochem%param%GrowDegDayEmerg         = NoahmpIO%GDDS1_TABLE    (CropType)
       noahmp%biochem%param%GrowDegDayInitVeg       = NoahmpIO%GDDS2_TABLE    (CropType)
       noahmp%biochem%param%GrowDegDayPostVeg       = NoahmpIO%GDDS3_TABLE    (CropType)
       noahmp%biochem%param%GrowDegDayInitReprod    = NoahmpIO%GDDS4_TABLE    (CropType)
       noahmp%biochem%param%GrowDegDayMature        = NoahmpIO%GDDS5_TABLE    (CropType)
       noahmp%biochem%param%PhotosynRadFrac         = NoahmpIO%I2PAR_TABLE    (CropType)
       noahmp%biochem%param%TempMinCarbonAssim      = NoahmpIO%TASSIM0_TABLE  (CropType)
       noahmp%biochem%param%TempMaxCarbonAssim      = NoahmpIO%TASSIM1_TABLE  (CropType)
       noahmp%biochem%param%TempMaxCarbonAssimMax   = NoahmpIO%TASSIM2_TABLE  (CropType)
       noahmp%biochem%param%CarbonAssimRefMax       = NoahmpIO%AREF_TABLE     (CropType)
       noahmp%biochem%param%LightExtCoeff           = NoahmpIO%K_TABLE        (CropType)
       noahmp%biochem%param%LightUseEfficiency      = NoahmpIO%EPSI_TABLE     (CropType)
       noahmp%biochem%param%CarbonAssimReducFac     = NoahmpIO%PSNRF_TABLE    (CropType)
       noahmp%biochem%param%RespMaintGrain25C       = NoahmpIO%GRAINMR25_TABLE(CropType)
       noahmp%biochem%param%LeafDeathTempCoeffCrop  = NoahmpIO%DILE_FC_TABLE  (CropType,:)
       noahmp%biochem%param%LeafDeathWaterCoeffCrop = NoahmpIO%DILE_FW_TABLE  (CropType,:)
       noahmp%biochem%param%CarbohydrLeafToGrain    = NoahmpIO%LFCT_TABLE     (CropType,:)
       noahmp%biochem%param%CarbohydrStemToGrain    = NoahmpIO%STCT_TABLE     (CropType,:)
       noahmp%biochem%param%CarbohydrRootToGrain    = NoahmpIO%RTCT_TABLE     (CropType,:)
       noahmp%biochem%param%CarbohydrFracToLeaf     = NoahmpIO%LFPT_TABLE     (CropType,:)
       noahmp%biochem%param%CarbohydrFracToStem     = NoahmpIO%STPT_TABLE     (CropType,:)
       noahmp%biochem%param%CarbohydrFracToRoot     = NoahmpIO%RTPT_TABLE     (CropType,:)
       noahmp%biochem%param%CarbohydrFracToGrain    = NoahmpIO%GRAINPT_TABLE  (CropType,:)
       noahmp%biochem%param%TurnoverCoeffLeafCrop   = NoahmpIO%LF_OVRC_TABLE  (CropType,:)
       noahmp%biochem%param%TurnoverCoeffStemCrop   = NoahmpIO%ST_OVRC_TABLE  (CropType,:)
       noahmp%biochem%param%TurnoverCoeffRootCrop   = NoahmpIO%RT_OVRC_TABLE  (CropType,:)

       if ( OptCropModel == 1 ) then
          if ( (NoahmpIO%PLANTING(N)>0) .and. (NoahmpIO%PLANTING(N)<367) ) then
             noahmp%biochem%param%DatePlanting      = NoahmpIO%PLANTING(N)
          endif ! 2D input map exist
          if ( (NoahmpIO%HARVEST(N)>0) .and. (NoahmpIO%HARVEST(N)<367) ) then
             noahmp%biochem%param%DateHarvest       = NoahmpIO%HARVEST(N)
          endif ! 2D input map exist
          if ( (NoahmpIO%SEASON_GDD(N)>0.0) .and. (NoahmpIO%SEASON_GDD(N)<10000.0) ) then
             noahmp%biochem%param%GrowDegDayEmerg   = NoahmpIO%SEASON_GDD(N) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayEmerg
             noahmp%biochem%param%GrowDegDayInitVeg = NoahmpIO%SEASON_GDD(N) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayInitVeg
             noahmp%biochem%param%GrowDegDayPostVeg = NoahmpIO%SEASON_GDD(N) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayPostVeg
             noahmp%biochem%param%GrowDegDayInitReprod = NoahmpIO%SEASON_GDD(N) / 1770.0 * &
                                                         noahmp%biochem%param%GrowDegDayInitReprod
             noahmp%biochem%param%GrowDegDayMature  = NoahmpIO%SEASON_GDD(N) / 1770.0 * &
                                                      noahmp%biochem%param%GrowDegDayMature
          endif ! 2D input map exist
       endif ! OptCropModel == 1
    endif ! activate crop parameters

    if ( noahmp%config%nmlist%OptIrrigation == 2 ) then
       if ( (NoahmpIO%PLANTING(N)>0) .and. (NoahmpIO%PLANTING(N)<367) ) then
          noahmp%biochem%param%DatePlanting = NoahmpIO%PLANTING(N)
       endif ! 2D input map exist
       if ( (NoahmpIO%HARVEST(N)>0) .and. (NoahmpIO%HARVEST(N)<367) ) then
          noahmp%biochem%param%DateHarvest  = NoahmpIO%HARVEST(N)
       endif ! 2D input map exist
    endif
    
    end associate

  end subroutine BiochemVarInTransfer

end module BiochemVarInTransferMod
