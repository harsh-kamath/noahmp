! Harsh Kamath (NOAA GSL/CIRES), September 2026


module EnergyVarOutTransferMod

  use Machine
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpVarType, only: noahmp_type

  implicit none

contains

!=== Transfer model states to output =====

  subroutine EnergyVarOutTransfer(noahmp, NoahmpIO, N)

    implicit none

    type(noahmp_type),        intent(inout) :: noahmp
    type(NoahmpIO_type),      intent(inout) :: NoahmpIO
	integer,                  intent(in)    :: N

    ! local variables
    integer                          :: LoopInd                   ! snow/soil layer loop index
    real(kind=kind_noahmp)           :: LeafAreaIndSunlit         ! sunlit leaf area index [m2/m2]
    real(kind=kind_noahmp)           :: LeafAreaIndShade          ! shaded leaf area index [m2/m2]
    real(kind=kind_noahmp)           :: ResistanceLeafBoundary    ! leaf boundary layer resistance [s/m]
    real(kind=kind_noahmp)           :: ThicknessSnowSoilLayer    ! temporary snow/soil layer thickness [m]

!-----------------------------------------------------------------------
    associate(                                                         &
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer    ,&
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSnowLayerNeg => noahmp%config%domain%NumSnowLayerNeg ,&
              NumSwRadBand    => noahmp%config%domain%NumSwRadBand    ,&
              IndicatorIceSfc => noahmp%config%domain%IndicatorIceSfc  &
             )
!-----------------------------------------------------------------------

    ! special treatment for glacier point output
    if ( IndicatorIceSfc == -1 ) then ! land ice point
       noahmp%energy%state%VegFrac             = 0.0
       noahmp%energy%state%RoughLenMomSfcToAtm = 0.002
       noahmp%energy%flux%RadSwAbsVeg          = 0.0
       noahmp%energy%flux%RadLwNetCanopy       = 0.0
       noahmp%energy%flux%RadLwNetVegGrd       = 0.0
       noahmp%energy%flux%HeatSensibleCanopy   = 0.0
       noahmp%energy%flux%HeatSensibleVegGrd   = 0.0
       noahmp%energy%flux%HeatLatentVegGrd     = 0.0
       noahmp%energy%flux%HeatGroundVegGrd     = 0.0
       noahmp%energy%flux%HeatCanStorageChg    = 0.0
       noahmp%energy%flux%HeatLatentCanTransp  = 0.0
       noahmp%energy%flux%HeatLatentCanEvap    = 0.0
       noahmp%energy%flux%HeatPrecipAdvCanopy  = 0.0
       noahmp%energy%flux%HeatPrecipAdvVegGrd  = 0.0
       noahmp%energy%flux%HeatLatentCanopy     = 0.0
       noahmp%energy%flux%HeatLatentTransp     = 0.0
       noahmp%energy%flux%RadLwNetBareGrd      = noahmp%energy%flux%RadLwNetSfc
       noahmp%energy%flux%HeatSensibleBareGrd  = noahmp%energy%flux%HeatSensibleSfc
       noahmp%energy%flux%HeatLatentBareGrd    = noahmp%energy%flux%HeatLatentGrd
       noahmp%energy%flux%HeatGroundBareGrd    = noahmp%energy%flux%HeatGroundTot
       noahmp%energy%state%TemperatureGrdBare  = noahmp%energy%state%TemperatureGrd
       noahmp%energy%state%ExchCoeffShBare     = noahmp%energy%state%ExchCoeffShSfc
       NoahmpIO%LH(N)                          = noahmp%energy%flux%HeatLatentGrd
    endif

    if ( IndicatorIceSfc == 0 ) then ! land soil point
       NoahmpIO%LH(N) = noahmp%energy%flux%HeatLatentGrd + noahmp%energy%flux%HeatLatentCanopy + &
                          noahmp%energy%flux%HeatLatentTransp + noahmp%energy%flux%HeatLatentIrriEvap 
    endif

! energy flux variables
    NoahmpIO%HFX        (N) = noahmp%energy%flux%HeatSensibleSfc
    NoahmpIO%GRDFLX     (N) = noahmp%energy%flux%HeatGroundTot
    NoahmpIO%FSAXY      (N) = noahmp%energy%flux%RadSwAbsSfc
    NoahmpIO%FIRAXY     (N) = noahmp%energy%flux%RadLwNetSfc
    NoahmpIO%APARXY     (N) = noahmp%energy%flux%RadPhotoActAbsCan
    NoahmpIO%SAVXY      (N) = noahmp%energy%flux%RadSwAbsVeg
    NoahmpIO%SAGXY      (N) = noahmp%energy%flux%RadSwAbsGrd
    NoahmpIO%IRCXY      (N) = noahmp%energy%flux%RadLwNetCanopy
    NoahmpIO%IRGXY      (N) = noahmp%energy%flux%RadLwNetVegGrd
    NoahmpIO%SHCXY      (N) = noahmp%energy%flux%HeatSensibleCanopy
    NoahmpIO%SHGXY      (N) = noahmp%energy%flux%HeatSensibleVegGrd
    NoahmpIO%EVGXY      (N) = noahmp%energy%flux%HeatLatentVegGrd
    NoahmpIO%GHVXY      (N) = noahmp%energy%flux%HeatGroundVegGrd
    NoahmpIO%IRBXY      (N) = noahmp%energy%flux%RadLwNetBareGrd
    NoahmpIO%SHBXY      (N) = noahmp%energy%flux%HeatSensibleBareGrd
    NoahmpIO%EVBXY      (N) = noahmp%energy%flux%HeatLatentBareGrd
    NoahmpIO%GHBXY      (N) = noahmp%energy%flux%HeatGroundBareGrd
    NoahmpIO%TRXY       (N) = noahmp%energy%flux%HeatLatentCanTransp
    NoahmpIO%EVCXY      (N) = noahmp%energy%flux%HeatLatentCanEvap
    NoahmpIO%CANHSXY    (N) = noahmp%energy%flux%HeatCanStorageChg
    NoahmpIO%PAHXY      (N) = noahmp%energy%flux%HeatPrecipAdvSfc
    NoahmpIO%PAHGXY     (N) = noahmp%energy%flux%HeatPrecipAdvVegGrd
    NoahmpIO%PAHVXY     (N) = noahmp%energy%flux%HeatPrecipAdvCanopy
    NoahmpIO%PAHBXY     (N) = noahmp%energy%flux%HeatPrecipAdvBareGrd
    NoahmpIO%ACC_SSOILXY(N) = noahmp%energy%flux%HeatGroundTotAcc
    NoahmpIO%EFLXBXY    (N) = noahmp%energy%flux%HeatFromSoilBot

    ! energy state variables
    NoahmpIO%TSK     (N) = noahmp%energy%state%TemperatureRadSfc
    NoahmpIO%EMISS   (N) = noahmp%energy%state%EmissivitySfc
    NoahmpIO%QSFC    (N) = noahmp%energy%state%SpecHumiditySfcMean
    NoahmpIO%TVXY    (N) = noahmp%energy%state%TemperatureCanopy
    NoahmpIO%TGXY    (N) = noahmp%energy%state%TemperatureGrd
    NoahmpIO%EAHXY   (N) = noahmp%energy%state%PressureVaporCanAir
    NoahmpIO%TAHXY   (N) = noahmp%energy%state%TemperatureCanopyAir
    NoahmpIO%CMXY    (N) = noahmp%energy%state%ExchCoeffMomSfc
    NoahmpIO%CHXY    (N) = noahmp%energy%state%ExchCoeffShSfc
    NoahmpIO%ALBOLDXY(N) = noahmp%energy%state%AlbedoSnowPrev
    NoahmpIO%LAI     (N) = noahmp%energy%state%LeafAreaIndex
    NoahmpIO%XSAIXY  (N) = noahmp%energy%state%StemAreaIndex
    NoahmpIO%TAUSSXY (N) = noahmp%energy%state%SnowAgeNondim
    NoahmpIO%Z0      (N) = noahmp%energy%state%RoughLenMomSfcToAtm
    NoahmpIO%ZNT     (N) = noahmp%energy%state%RoughLenMomSfcToAtm
    NoahmpIO%T2MVXY  (N) = noahmp%energy%state%TemperatureAir2mVeg
    NoahmpIO%T2MBXY  (N) = noahmp%energy%state%TemperatureAir2mBare
    NoahmpIO%TRADXY  (N) = noahmp%energy%state%TemperatureRadSfc
    NoahmpIO%FVEGXY  (N) = noahmp%energy%state%VegFrac
    NoahmpIO%RSSUNXY (N) = noahmp%energy%state%ResistanceStomataSunlit
    NoahmpIO%RSSHAXY (N) = noahmp%energy%state%ResistanceStomataShade
    NoahmpIO%BGAPXY  (N) = noahmp%energy%state%GapBtwCanopy
    NoahmpIO%WGAPXY  (N) = noahmp%energy%state%GapInCanopy
    NoahmpIO%TGVXY   (N) = noahmp%energy%state%TemperatureGrdVeg
    NoahmpIO%TGBXY   (N) = noahmp%energy%state%TemperatureGrdBare
    NoahmpIO%CHVXY   (N) = noahmp%energy%state%ExchCoeffShAbvCan
    NoahmpIO%CHBXY   (N) = noahmp%energy%state%ExchCoeffShBare
    NoahmpIO%CHLEAFXY(N) = noahmp%energy%state%ExchCoeffShLeaf
    NoahmpIO%CHUCXY  (N) = noahmp%energy%state%ExchCoeffShUndCan
    NoahmpIO%CHV2XY  (N) = noahmp%energy%state%ExchCoeffSh2mVeg
    NoahmpIO%CHB2XY  (N) = noahmp%energy%state%ExchCoeffSh2mBare
    NoahmpIO%Q2MVXY  (N) = noahmp%energy%state%SpecHumidity2mVeg /(1.0-noahmp%energy%state%SpecHumidity2mVeg)  ! spec humidity to mixing ratio
    NoahmpIO%Q2MBXY  (N) = noahmp%energy%state%SpecHumidity2mBare/(1.0-noahmp%energy%state%SpecHumidity2mBare)
    NoahmpIO%ALBEDO  (N) = noahmp%energy%state%AlbedoSfc
    NoahmpIO%IRRSPLH (N) = NoahmpIO%IRRSPLH(N) + &
                          (noahmp%energy%flux%HeatLatentIrriEvap * noahmp%config%domain%MainTimeStep)
    NoahmpIO%TSLB    (1:NumSoilLayer,N)        = noahmp%energy%state%TemperatureSoilSnow(1:NumSoilLayer)
    NoahmpIO%TSNOXY  (-NumSnowLayerMax+1:0,N) = noahmp%energy%state%TemperatureSoilSnow(-NumSnowLayerMax+1:0)

    NoahmpIO%ALBSOILDIRXY(1:NumSwRadBand,N) = noahmp%energy%state%AlbedoSoilDir(1:NumSwRadBand)
    NoahmpIO%ALBSOILDIFXY(1:NumSwRadBand,N) = noahmp%energy%state%AlbedoSoilDif(1:NumSwRadBand)
    NoahmpIO%ALBSFCDIRXY (1:NumSwRadBand,N) = noahmp%energy%state%AlbedoSfcDir (1:NumSwRadBand)
    NoahmpIO%ALBSFCDIFXY (1:NumSwRadBand,N) = noahmp%energy%state%AlbedoSfcDif (1:NumSwRadBand)
    NoahmpIO%ALBSNOWDIRXY(1:NumSwRadBand,N) = noahmp%energy%state%AlbedoSnowDir(1:NumSwRadBand)
    NoahmpIO%ALBSNOWDIFXY(1:NumSwRadBand,N) = noahmp%energy%state%AlbedoSnowDif(1:NumSwRadBand)

    ! New Calculation of total Canopy/Stomatal Conductance Based on Bonan et al. (2011), Inverse of Canopy Resistance (below)
    LeafAreaIndSunlit      = max(noahmp%energy%state%LeafAreaIndSunlit, 0.0)
    LeafAreaIndShade       = max(noahmp%energy%state%LeafAreaIndShade, 0.0)
    ResistanceLeafBoundary = max(noahmp%energy%state%ResistanceLeafBoundary, 0.0)
    if ( (noahmp%energy%state%ResistanceStomataSunlit <= 0.0) .or. (noahmp%energy%state%ResistanceStomataShade <= 0.0) .or. &
         (LeafAreaIndSunlit == 0.0) .or. (LeafAreaIndShade == 0.0)       .or. &
         (noahmp%energy%state%ResistanceStomataSunlit == undefined_real) .or. &
         (noahmp%energy%state%ResistanceStomataShade == undefined_real) ) then
       NoahmpIO%RS   (N) = 0.0
    else
       NoahmpIO%RS   (N) = ((1.0 / (noahmp%energy%state%ResistanceStomataSunlit + ResistanceLeafBoundary) * &
                            noahmp%energy%state%LeafAreaIndSunlit) + &
                           ((1.0 / (noahmp%energy%state%ResistanceStomataShade + ResistanceLeafBoundary)) * &
                            noahmp%energy%state%LeafAreaIndShade))
       NoahmpIO%RS   (N) = 1.0 / NoahmpIO%RS (N) ! Resistance
    endif

    ! calculation of snow and soil energy storage
    NoahmpIO%SNOWENERGY(N) = 0.0
    NoahmpIO%SOILENERGY(N) = 0.0
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       if ( LoopInd == NumSnowLayerNeg+1 ) then
          ThicknessSnowSoilLayer = -noahmp%config%domain%DepthSnowSoilLayer(LoopInd)
       else
          ThicknessSnowSoilLayer = noahmp%config%domain%DepthSnowSoilLayer(LoopInd-1) - &
                                   noahmp%config%domain%DepthSnowSoilLayer(LoopInd)
       endif
       if ( LoopInd >= 1 ) then
          NoahmpIO%SOILENERGY(N) = NoahmpIO%SOILENERGY(N) + ThicknessSnowSoilLayer * &
                                     noahmp%energy%state%HeatCapacSoilSnow(LoopInd) * &
                                     (noahmp%energy%state%TemperatureSoilSnow(LoopInd) - 273.16) * 0.001
       else
          NoahmpIO%SNOWENERGY(N) = NoahmpIO%SNOWENERGY(N) + ThicknessSnowSoilLayer * &
                                     noahmp%energy%state%HeatCapacSoilSnow(LoopInd) * &
                                     (noahmp%energy%state%TemperatureSoilSnow(LoopInd) - 273.16) * 0.001
       endif
    enddo

    end associate

  end subroutine EnergyVarOutTransfer

end module EnergyVarOutTransferMod
