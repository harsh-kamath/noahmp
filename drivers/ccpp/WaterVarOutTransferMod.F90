! Harsh Kamath (NOAA GSL/CIRES), September 2026


module WaterVarOutTransferMod

  use Machine
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpVarType, only: noahmp_type

  implicit none

contains

!=== Transfer model states to output =====

  subroutine WaterVarOutTransfer(noahmp, NoahmpIO, N)

    implicit none

    type(noahmp_type),        intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO
	integer,                  intent(in)    :: N

! -------------------------------------------------------------------------
    associate(                                                         &
              ! I               => noahmp%config%domain%GridIndexI      ,&
              ! J               => noahmp%config%domain%GridIndexJ      ,&
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer    ,&
              IndicatorIceSfc => noahmp%config%domain%IndicatorIceSfc  &
             )
! -------------------------------------------------------------------------

    ! special treatment for glacier point output
    if ( IndicatorIceSfc == -1 ) then ! land ice point
       noahmp%water%state%SnowCoverFrac      = 1.0
       noahmp%water%flux%EvapCanopyNet       = 0.0
       noahmp%water%flux%Transpiration       = 0.0
       noahmp%water%flux%InterceptCanopySnow = 0.0
       noahmp%water%flux%InterceptCanopyRain = 0.0
       noahmp%water%flux%DripCanopySnow      = 0.0
       noahmp%water%flux%DripCanopyRain      = 0.0
       noahmp%water%flux%ThroughfallSnow     = noahmp%water%flux%SnowfallRefHeight
       noahmp%water%flux%ThroughfallRain     = noahmp%water%flux%RainfallRefHeight
       noahmp%water%flux%SublimCanopyIce     = 0.0
       noahmp%water%flux%FrostCanopyIce      = 0.0
       noahmp%water%flux%FreezeCanopyLiq     = 0.0
       noahmp%water%flux%MeltCanopyIce       = 0.0
       noahmp%water%flux%EvapCanopyLiq       = 0.0
       noahmp%water%flux%DewCanopyLiq        = 0.0
       noahmp%water%state%CanopyIce          = 0.0
       noahmp%water%state%CanopyLiqWater     = 0.0
       noahmp%water%flux%TileDrain           = 0.0
       noahmp%water%flux%RunoffSurface       = noahmp%water%flux%RunoffSurface * noahmp%config%domain%MainTimeStep
       noahmp%water%flux%RunoffSubsurface    = noahmp%water%flux%RunoffSubsurface * noahmp%config%domain%MainTimeStep
       NoahmpIO%QFX(N)                          = noahmp%water%flux%EvapGroundNet
    endif

    if ( IndicatorIceSfc == 0 ) then ! land soil point
       NoahmpIO%QFX(N) = noahmp%water%flux%EvapCanopyNet + noahmp%water%flux%EvapGroundNet + &
                           noahmp%water%flux%Transpiration + noahmp%water%flux%EvapIrriSprinkler
    endif

    NoahmpIO%SMSTAV      (N) = 0.0  ! [maintained as Noah consistency] water
    NoahmpIO%SMSTOT      (N) = 0.0  ! [maintained as Noah consistency] water
    NoahmpIO%SFCRUNOFF   (N) = NoahmpIO%SFCRUNOFF(N) + noahmp%water%flux%RunoffSurface
    NoahmpIO%UDRUNOFF    (N) = NoahmpIO%UDRUNOFF (N) + noahmp%water%flux%RunoffSubsurface
    NoahmpIO%QTDRAIN     (N) = NoahmpIO%QTDRAIN  (N) + noahmp%water%flux%TileDrain
    NoahmpIO%SNOWC       (N) = noahmp%water%state%SnowCoverFrac
    NoahmpIO%SNOW        (N) = noahmp%water%state%SnowWaterEquiv
    NoahmpIO%SNOWH       (N) = noahmp%water%state%SnowDepth
    NoahmpIO%CANWAT      (N) = noahmp%water%state%CanopyLiqWater + noahmp%water%state%CanopyIce
    NoahmpIO%ACSNOW      (N) = NoahmpIO%ACSNOW(N) + NoahmpIO%RAINBL (N) * noahmp%water%state%FrozenPrecipFrac
    NoahmpIO%ACSNOM      (N) = NoahmpIO%ACSNOM(N) + noahmp%water%flux%MeltGroundSnow * NoahmpIO%DTBL
    NoahmpIO%CANLIQXY    (N) = noahmp%water%state%CanopyLiqWater
    NoahmpIO%CANICEXY    (N) = noahmp%water%state%CanopyIce
    NoahmpIO%FWETXY      (N) = noahmp%water%state%CanopyWetFrac
    NoahmpIO%SNEQVOXY    (N) = noahmp%water%state%SnowWaterEquivPrev
    NoahmpIO%QSNOWXY     (N) = noahmp%water%flux%SnowfallGround
    NoahmpIO%QRAINXY     (N) = noahmp%water%flux%RainfallGround
    NoahmpIO%WSLAKEXY    (N) = noahmp%water%state%WaterStorageLake
    NoahmpIO%ZWTXY       (N) = noahmp%water%state%WaterTableDepth
    NoahmpIO%WAXY        (N) = noahmp%water%state%WaterStorageAquifer
    NoahmpIO%WTXY        (N) = noahmp%water%state%WaterStorageSoilAqf
    NoahmpIO%RUNSFXY     (N) = noahmp%water%flux%RunoffSurface
    NoahmpIO%RUNSBXY     (N) = noahmp%water%flux%RunoffSubsurface
    NoahmpIO%ECANXY      (N) = noahmp%water%flux%EvapCanopyNet
    NoahmpIO%EDIRXY      (N) = noahmp%water%flux%EvapGroundNet
    NoahmpIO%ETRANXY     (N) = noahmp%water%flux%Transpiration
    NoahmpIO%QINTSXY     (N) = noahmp%water%flux%InterceptCanopySnow
    NoahmpIO%QINTRXY     (N) = noahmp%water%flux%InterceptCanopyRain
    NoahmpIO%QDRIPSXY    (N) = noahmp%water%flux%DripCanopySnow
    NoahmpIO%QDRIPRXY    (N) = noahmp%water%flux%DripCanopyRain
    NoahmpIO%QTHROSXY    (N) = noahmp%water%flux%ThroughfallSnow
    NoahmpIO%QTHRORXY    (N) = noahmp%water%flux%ThroughfallRain
    NoahmpIO%QSNSUBXY    (N) = noahmp%water%flux%SublimSnowSfcIce
    NoahmpIO%QSNFROXY    (N) = noahmp%water%flux%FrostSnowSfcIce
    NoahmpIO%QSUBCXY     (N) = noahmp%water%flux%SublimCanopyIce
    NoahmpIO%QFROCXY     (N) = noahmp%water%flux%FrostCanopyIce
    NoahmpIO%QEVACXY     (N) = noahmp%water%flux%EvapCanopyLiq
    NoahmpIO%QDEWCXY     (N) = noahmp%water%flux%DewCanopyLiq
    NoahmpIO%QFRZCXY     (N) = noahmp%water%flux%FreezeCanopyLiq
    NoahmpIO%QMELTCXY    (N) = noahmp%water%flux%MeltCanopyIce
    NoahmpIO%QSNBOTXY    (N) = noahmp%water%flux%SnowBotOutflow
    NoahmpIO%QMELTXY     (N) = noahmp%water%flux%MeltGroundSnow
    NoahmpIO%PONDINGXY   (N) = noahmp%water%state%PondSfcThinSnwTrans + &
                            noahmp%water%state%PondSfcThinSnwComb + noahmp%water%state%PondSfcThinSnwMelt
    NoahmpIO%FPICEXY     (N) = noahmp%water%state%FrozenPrecipFrac
    NoahmpIO%RAINLSM     (N) = noahmp%water%flux%RainfallRefHeight
    NoahmpIO%SNOWLSM     (N) = noahmp%water%flux%SnowfallRefHeight
    NoahmpIO%ACC_QINSURXY(N) = noahmp%water%flux%SoilSfcInflowAcc
    NoahmpIO%ACC_QSEVAXY (N) = noahmp%water%flux%EvapSoilSfcLiqAcc
    NoahmpIO%ACC_DWATERXY(N) = noahmp%water%flux%SfcWaterTotChgAcc
    NoahmpIO%ACC_PRCPXY  (N) = noahmp%water%flux%PrecipTotAcc
    NoahmpIO%ACC_ECANXY  (N) = noahmp%water%flux%EvapCanopyNetAcc
    NoahmpIO%ACC_ETRANXY (N) = noahmp%water%flux%TranspirationAcc
    NoahmpIO%ACC_EDIRXY  (N) = noahmp%water%flux%EvapGroundNetAcc
    NoahmpIO%ACC_GLAFLWXY(N) = noahmp%water%flux%GlacierExcessFlowAcc
    NoahmpIO%RECHXY      (N) = NoahmpIO%RECHXY(N) + (noahmp%water%state%RechargeGwShallowWT*1.0e3)
    NoahmpIO%DEEPRECHXY  (N) = NoahmpIO%DEEPRECHXY(N) + noahmp%water%state%RechargeGwDeepWT
    NoahmpIO%SMCWTDXY    (N) = noahmp%water%state%SoilMoistureToWT
    NoahmpIO%SMOIS       (1:NumSoilLayer,N)       = noahmp%water%state%SoilMoisture(1:NumSoilLayer)
    NoahmpIO%SH2O        (1:NumSoilLayer,N)       = noahmp%water%state%SoilLiqWater(1:NumSoilLayer)
    NoahmpIO%ACC_ETRANIXY(1:NumSoilLayer,N)       = noahmp%water%flux%TranspWatLossSoilAcc(1:NumSoilLayer)
    NoahmpIO%SNICEXY     (-NumSnowLayerMax+1:0,N) = noahmp%water%state%SnowIce(-NumSnowLayerMax+1:0)
    NoahmpIO%SNLIQXY     (-NumSnowLayerMax+1:0,N) = noahmp%water%state%SnowLiqWater(-NumSnowLayerMax+1:0)

    !SNICAR
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       NoahmpIO%SNRDSXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%SnowRadius(-NumSnowLayerMax+1:0)
       NoahmpIO%SNFRXY (-NumSnowLayerMax+1:0,N) = noahmp%water%flux%SnowFreezeRate(-NumSnowLayerMax+1:0)
       NoahmpIO%BCPHIXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassBChydrophi(-NumSnowLayerMax+1:0)
       NoahmpIO%BCPHOXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassBChydropho(-NumSnowLayerMax+1:0)
       NoahmpIO%OCPHIXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassOChydrophi(-NumSnowLayerMax+1:0)
       NoahmpIO%OCPHOXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassOChydropho(-NumSnowLayerMax+1:0)
       NoahmpIO%DUST1XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassDust1(-NumSnowLayerMax+1:0)
       NoahmpIO%DUST2XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassDust2(-NumSnowLayerMax+1:0)
       NoahmpIO%DUST3XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassDust3(-NumSnowLayerMax+1:0)
       NoahmpIO%DUST4XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassDust4(-NumSnowLayerMax+1:0)
       NoahmpIO%DUST5XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassDust5(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcBCPHIXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcBChydrophi(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcBCPHOXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcBChydropho(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcOCPHIXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcOChydrophi(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcOCPHOXY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcOChydropho(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcDUST1XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcDust1(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcDUST2XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcDust2(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcDUST3XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcDust3(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcDUST4XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcDust4(-NumSnowLayerMax+1:0)
       NoahmpIO%MassConcDUST5XY(-NumSnowLayerMax+1:0,N) = noahmp%water%state%MassConcDust5(-NumSnowLayerMax+1:0)
    endif

    ! irrigation
    NoahmpIO%IRNUMSI(N) = noahmp%water%state%IrrigationCntSprinkler
    NoahmpIO%IRNUMMI(N) = noahmp%water%state%IrrigationCntMicro
    NoahmpIO%IRNUMFI(N) = noahmp%water%state%IrrigationCntFlood
    NoahmpIO%IRWATSI(N) = noahmp%water%state%IrrigationAmtSprinkler
    NoahmpIO%IRWATMI(N) = noahmp%water%state%IrrigationAmtMicro
    NoahmpIO%IRWATFI(N) = noahmp%water%state%IrrigationAmtFlood
    NoahmpIO%IRSIVOL(N) = NoahmpIO%IRSIVOL(N) + (noahmp%water%flux%IrrigationRateSprinkler*1000.0)
    NoahmpIO%IRMIVOL(N) = NoahmpIO%IRMIVOL(N) + (noahmp%water%flux%IrrigationRateMicro*1000.0)
    NoahmpIO%IRFIVOL(N) = NoahmpIO%IRFIVOL(N) + (noahmp%water%flux%IrrigationRateFlood*1000.0)
    NoahmpIO%IRELOSS(N) = NoahmpIO%IRELOSS(N) + (noahmp%water%flux%EvapIrriSprinkler*NoahmpIO%DTBL)

    ! wetland (Zhang2022)
    if ( noahmp%config%nmlist%OptWetlandModel > 0 ) then
       NoahmpIO%WSURFXY(N) = noahmp%water%state%WaterStorageWetland
       NoahmpIO%FSATXY (N) = noahmp%water%state%SoilSaturateFrac
    endif


    end associate

  end subroutine WaterVarOutTransfer

end module WaterVarOutTransferMod
