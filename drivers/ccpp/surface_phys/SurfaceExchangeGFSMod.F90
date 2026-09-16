module SurfaceExchangeGFSMod

! GFS surface-layer exchange coefficients for the modular Noah-MP CCPP
! driver.  Conductance outputs use Noah-MP units [m/s], rather than the
! dimensionless coefficients returned by the legacy sfcdif3 interface.
! Harsh Kamath (NOAA GSL/CIRES), September 2026


  use Machine, only : kind_noahmp
  use NoahmpVarType, only : noahmp_type
  use ConstantDefineMod, only : ConstGravityAcc, ConstGasDryAir, &
       ConstHeatCapacAir

  implicit none
  private

  public :: SurfaceExchangeGFS
  public :: UpdateNoahmpSurfaceExchangeGFS

contains

  subroutine UpdateNoahmpSurfaceExchangeGFS(noahmp, SurfaceTemperature, IsVegetated)

    type(noahmp_type),      intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: SurfaceTemperature
    logical,                intent(in)    :: IsVegetated

    real(kind=kind_noahmp) :: MomentumRoughnessLength
    real(kind=kind_noahmp) :: HeatRoughnessLength
    real(kind=kind_noahmp) :: ReferenceHeight
    real(kind=kind_noahmp) :: VirtualTemperatureFactor
    real(kind=kind_noahmp) :: VirtualAirTemperature
    real(kind=kind_noahmp) :: VirtualAirPotentialTemperature
    real(kind=kind_noahmp) :: VirtualSurfaceTemperature
    real(kind=kind_noahmp) :: PotentialTemperatureFactor
    real(kind=kind_noahmp) :: WindSpeed
    real(kind=kind_noahmp) :: BulkRichardsonNumber
    real(kind=kind_noahmp) :: MomentumDenominator
    real(kind=kind_noahmp) :: HeatDenominator
    real(kind=kind_noahmp) :: MomentumDenominator10m
    real(kind=kind_noahmp) :: HeatDenominator2m
    real(kind=kind_noahmp) :: MomentumConductance
    real(kind=kind_noahmp) :: HeatConductance
    real(kind=kind_noahmp) :: SurfaceStress
    real(kind=kind_noahmp) :: FrictionVelocity

    ReferenceHeight = max(0.01_kind_noahmp, noahmp%energy%state%RefHeightAboveGrd - &
         noahmp%energy%state%ZeroPlaneDispSfc)
    WindSpeed = max(0.1_kind_noahmp, noahmp%energy%state%WindSpdRefHeight)
    VirtualTemperatureFactor = 1.0_kind_noahmp + &
         0.61_kind_noahmp*max(noahmp%forcing%SpecHumidityRefHeight, 1.0e-8_kind_noahmp)
    VirtualAirTemperature = noahmp%forcing%TemperatureAirRefHeight * &
         VirtualTemperatureFactor

    if ( noahmp%config%domain%FlagUseLocalPotTemp ) then
       PotentialTemperatureFactor = (noahmp%forcing%PressureAirSurface / &
            noahmp%forcing%PressureAirRefHeight)** &
            (ConstGasDryAir/ConstHeatCapacAir)
       VirtualSurfaceTemperature = SurfaceTemperature*VirtualTemperatureFactor
    else
       PotentialTemperatureFactor = (100000.0_kind_noahmp / &
            noahmp%forcing%PressureAirRefHeight)** &
            (ConstGasDryAir/ConstHeatCapacAir)
       VirtualSurfaceTemperature = SurfaceTemperature * &
            (100000.0_kind_noahmp/noahmp%forcing%PressureAirSurface)** &
            (ConstGasDryAir/ConstHeatCapacAir) * VirtualTemperatureFactor
    endif
    VirtualAirPotentialTemperature = noahmp%forcing%TemperatureAirRefHeight * &
         PotentialTemperatureFactor*VirtualTemperatureFactor

    if ( IsVegetated ) then
       MomentumRoughnessLength = noahmp%energy%state%RoughLenMomSfc
       HeatRoughnessLength = noahmp%energy%state%RoughLenShCanopy
    else
       MomentumRoughnessLength = noahmp%energy%state%RoughLenMomGrd
       HeatRoughnessLength = noahmp%energy%state%RoughLenShBareGrd
    endif

    call SurfaceExchangeGFS(ReferenceHeight, 1.0_kind_noahmp, &
         3000.0_kind_noahmp, VirtualAirTemperature, &
         VirtualAirPotentialTemperature, WindSpeed, MomentumRoughnessLength, &
         HeatRoughnessLength, VirtualSurfaceTemperature, ConstGravityAcc, &
         noahmp%config%domain%FlagUseLocalPotTemp, BulkRichardsonNumber, &
         MomentumDenominator, HeatDenominator, MomentumDenominator10m, &
         HeatDenominator2m, MomentumConductance, HeatConductance, &
         SurfaceStress, FrictionVelocity)

    if ( IsVegetated ) then
       noahmp%energy%state%ExchCoeffMomAbvCan = MomentumConductance/WindSpeed
       noahmp%energy%state%ExchCoeffShAbvCan = HeatConductance/WindSpeed
       noahmp%energy%state%ResistanceMomAbvCan = max(1.0_kind_noahmp, &
            1.0_kind_noahmp/MomentumConductance)
       noahmp%energy%state%ResistanceShAbvCan = max(1.0_kind_noahmp, &
            1.0_kind_noahmp/HeatConductance)
       noahmp%energy%state%ResistanceLhAbvCan = &
            noahmp%energy%state%ResistanceShAbvCan
       noahmp%energy%state%FrictionVelVeg = FrictionVelocity
       noahmp%energy%state%ExchCoeffSh2mVeg = &
            FrictionVelocity*0.4_kind_noahmp/HeatDenominator2m
    else
       noahmp%energy%state%ExchCoeffMomBare = MomentumConductance/WindSpeed
       noahmp%energy%state%ExchCoeffShBare = HeatConductance/WindSpeed
       noahmp%energy%state%ResistanceMomBareGrd = max(1.0_kind_noahmp, &
            1.0_kind_noahmp/MomentumConductance)
       noahmp%energy%state%ResistanceShBareGrd = max(1.0_kind_noahmp, &
            1.0_kind_noahmp/HeatConductance)
       noahmp%energy%state%ResistanceLhBareGrd = &
            noahmp%energy%state%ResistanceShBareGrd
       noahmp%energy%state%FrictionVelBare = FrictionVelocity
       noahmp%energy%state%ExchCoeffSh2mBare = &
            FrictionVelocity*0.4_kind_noahmp/HeatDenominator2m
    endif

  end subroutine UpdateNoahmpSurfaceExchangeGFS

  subroutine SurfaceExchangeGFS(ReferenceHeight, VegetationAdjustment, &
       HorizontalGridLength, VirtualAirTemperature, &
       VirtualAirPotentialTemperature, WindSpeed, MomentumRoughnessLength, &
       HeatRoughnessLength, VirtualSurfaceTemperature, GravityAcceleration, &
       UseLocalPotentialTemperature, BulkRichardsonNumber, &
       MomentumSimilarityDenominator, HeatSimilarityDenominator, &
       MomentumSimilarityDenominator10m, HeatSimilarityDenominator2m, &
       MomentumConductance, HeatConductance, SurfaceStress, FrictionVelocity)

    real(kind=kind_noahmp), intent(in)  :: ReferenceHeight
    real(kind=kind_noahmp), intent(in)  :: VegetationAdjustment
    real(kind=kind_noahmp), intent(in)  :: HorizontalGridLength
    real(kind=kind_noahmp), intent(in)  :: VirtualAirTemperature
    real(kind=kind_noahmp), intent(in)  :: VirtualAirPotentialTemperature
    real(kind=kind_noahmp), intent(in)  :: WindSpeed
    real(kind=kind_noahmp), intent(in)  :: MomentumRoughnessLength
    real(kind=kind_noahmp), intent(in)  :: HeatRoughnessLength
    real(kind=kind_noahmp), intent(in)  :: VirtualSurfaceTemperature
    real(kind=kind_noahmp), intent(in)  :: GravityAcceleration
    logical,                intent(in)  :: UseLocalPotentialTemperature
    real(kind=kind_noahmp), intent(out) :: BulkRichardsonNumber
    real(kind=kind_noahmp), intent(out) :: MomentumSimilarityDenominator
    real(kind=kind_noahmp), intent(out) :: HeatSimilarityDenominator
    real(kind=kind_noahmp), intent(out) :: MomentumSimilarityDenominator10m
    real(kind=kind_noahmp), intent(out) :: HeatSimilarityDenominator2m
    real(kind=kind_noahmp), intent(out) :: MomentumConductance
    real(kind=kind_noahmp), intent(out) :: HeatConductance
    real(kind=kind_noahmp), intent(out) :: SurfaceStress
    real(kind=kind_noahmp), intent(out) :: FrictionVelocity

    real(kind=kind_noahmp), parameter :: VonKarmanConstant = 0.4_kind_noahmp
    real(kind=kind_noahmp), parameter :: StableAlpha = 5.0_kind_noahmp
    real(kind=kind_noahmp), parameter :: ReferenceMaximumStability = 0.3_kind_noahmp
    real(kind=kind_noahmp), parameter :: MinimumDiffusivityFactor = 0.05_kind_noahmp
    real(kind=kind_noahmp), parameter :: CriticalGridLength = 3000.0_kind_noahmp
    real(kind=kind_noahmp), parameter :: MinimumStability = -10.0_kind_noahmp
    real(kind=kind_noahmp), parameter :: A0Momentum = -3.975_kind_noahmp
    real(kind=kind_noahmp), parameter :: A1Momentum = 12.32_kind_noahmp
    real(kind=kind_noahmp), parameter :: B1Momentum = -7.755_kind_noahmp
    real(kind=kind_noahmp), parameter :: B2Momentum = 6.041_kind_noahmp
    real(kind=kind_noahmp), parameter :: A0Heat = -7.941_kind_noahmp
    real(kind=kind_noahmp), parameter :: A1Heat = 24.75_kind_noahmp
    real(kind=kind_noahmp), parameter :: B1Heat = -8.705_kind_noahmp
    real(kind=kind_noahmp), parameter :: B2Heat = 7.899_kind_noahmp

    real(kind=kind_noahmp) :: AirSurfaceTemperatureDifference
    real(kind=kind_noahmp) :: AbsoluteTemperatureDifference
    real(kind=kind_noahmp) :: DiffusivityFactor
    real(kind=kind_noahmp) :: MaximumStability
    real(kind=kind_noahmp) :: MoninObukhovParameter
    real(kind=kind_noahmp) :: MomentumCorrection
    real(kind=kind_noahmp) :: HeatCorrection
    real(kind=kind_noahmp) :: MomentumCorrection10m
    real(kind=kind_noahmp) :: HeatCorrection2m
    real(kind=kind_noahmp) :: ParameterAtMomentumRoughness
    real(kind=kind_noahmp) :: ParameterAtHeatRoughness
    real(kind=kind_noahmp) :: InverseReferenceHeight
    real(kind=kind_noahmp) :: RootAtReferenceHeight
    real(kind=kind_noahmp) :: RootAtMomentumRoughness
    real(kind=kind_noahmp) :: RootAtHeatRoughness
    real(kind=kind_noahmp) :: DimensionlessMomentumCoefficient
    real(kind=kind_noahmp) :: DimensionlessHeatCoefficient
    real(kind=kind_noahmp) :: Temporary

    InverseReferenceHeight = 1.0_kind_noahmp / ReferenceHeight
    DiffusivityFactor = min(1.0_kind_noahmp, &
         HorizontalGridLength/CriticalGridLength)
    if ( VirtualAirTemperature > VirtualSurfaceTemperature ) then
       DiffusivityFactor = min(max(DiffusivityFactor*VegetationAdjustment, &
            MinimumDiffusivityFactor), DiffusivityFactor)
    endif
    MaximumStability = ReferenceMaximumStability / sqrt(DiffusivityFactor)

    AirSurfaceTemperatureDifference = VirtualAirPotentialTemperature - &
         VirtualSurfaceTemperature
    AbsoluteTemperatureDifference = max(abs(AirSurfaceTemperatureDifference), &
         0.001_kind_noahmp)
    AirSurfaceTemperatureDifference = sign(AbsoluteTemperatureDifference, &
         AirSurfaceTemperatureDifference)

    if ( UseLocalPotentialTemperature ) then
       BulkRichardsonNumber = max(-5000.0_kind_noahmp, &
            2.0_kind_noahmp*GravityAcceleration*AirSurfaceTemperatureDifference* &
            ReferenceHeight / ((VirtualAirPotentialTemperature + &
            VirtualSurfaceTemperature)*WindSpeed**2))
    else
       BulkRichardsonNumber = max(-5000.0_kind_noahmp, &
            GravityAcceleration*AirSurfaceTemperatureDifference*ReferenceHeight / &
            (VirtualAirTemperature*WindSpeed**2))
    endif

    MomentumSimilarityDenominator = log((MomentumRoughnessLength+ReferenceHeight) / &
         MomentumRoughnessLength)
    HeatSimilarityDenominator = log((HeatRoughnessLength+ReferenceHeight) / &
         HeatRoughnessLength)
    MomentumSimilarityDenominator10m = log((MomentumRoughnessLength+10.0_kind_noahmp) / &
         MomentumRoughnessLength)
    HeatSimilarityDenominator2m = log((HeatRoughnessLength+2.0_kind_noahmp) / &
         HeatRoughnessLength)

    MoninObukhovParameter = BulkRichardsonNumber * &
         MomentumSimilarityDenominator**2 / HeatSimilarityDenominator
    MoninObukhovParameter = min(max(MoninObukhovParameter, MinimumStability), &
         MaximumStability)

    if ( AirSurfaceTemperatureDifference >= 0.0_kind_noahmp ) then
       if ( MoninObukhovParameter > 0.25_kind_noahmp ) then
          ParameterAtMomentumRoughness = MomentumRoughnessLength * &
               MoninObukhovParameter*InverseReferenceHeight
          ParameterAtHeatRoughness = HeatRoughnessLength * &
               MoninObukhovParameter*InverseReferenceHeight
          RootAtReferenceHeight = sqrt(1.0_kind_noahmp + &
               4.0_kind_noahmp*StableAlpha*MoninObukhovParameter)
          RootAtMomentumRoughness = sqrt(1.0_kind_noahmp + &
               4.0_kind_noahmp*StableAlpha*ParameterAtMomentumRoughness)
          RootAtHeatRoughness = sqrt(1.0_kind_noahmp + &
               4.0_kind_noahmp*StableAlpha*ParameterAtHeatRoughness)
          MomentumCorrection = RootAtMomentumRoughness-RootAtReferenceHeight + &
               log((RootAtReferenceHeight+1.0_kind_noahmp) / &
                   (RootAtMomentumRoughness+1.0_kind_noahmp))
          HeatCorrection = RootAtHeatRoughness-RootAtReferenceHeight + &
               log((RootAtReferenceHeight+1.0_kind_noahmp) / &
                   (RootAtHeatRoughness+1.0_kind_noahmp))
          MoninObukhovParameter = min(MaximumStability, &
               BulkRichardsonNumber * &
               (MomentumSimilarityDenominator-MomentumCorrection)**2 / &
               (HeatSimilarityDenominator-HeatCorrection))
       endif

       ParameterAtMomentumRoughness = MomentumRoughnessLength * &
            MoninObukhovParameter*InverseReferenceHeight
       ParameterAtHeatRoughness = HeatRoughnessLength * &
            MoninObukhovParameter*InverseReferenceHeight
       RootAtReferenceHeight = sqrt(1.0_kind_noahmp + &
            4.0_kind_noahmp*StableAlpha*MoninObukhovParameter)
       RootAtMomentumRoughness = sqrt(1.0_kind_noahmp + &
            4.0_kind_noahmp*StableAlpha*ParameterAtMomentumRoughness)
       RootAtHeatRoughness = sqrt(1.0_kind_noahmp + &
            4.0_kind_noahmp*StableAlpha*ParameterAtHeatRoughness)
       MomentumCorrection = RootAtMomentumRoughness-RootAtReferenceHeight + &
            log((RootAtReferenceHeight+1.0_kind_noahmp) / &
                (RootAtMomentumRoughness+1.0_kind_noahmp))
       HeatCorrection = RootAtHeatRoughness-RootAtReferenceHeight + &
            log((RootAtReferenceHeight+1.0_kind_noahmp) / &
                (RootAtHeatRoughness+1.0_kind_noahmp))
       Temporary = MoninObukhovParameter*10.0_kind_noahmp*InverseReferenceHeight
       RootAtReferenceHeight = sqrt(1.0_kind_noahmp + &
            4.0_kind_noahmp*StableAlpha*Temporary)
       MomentumCorrection10m = RootAtMomentumRoughness-RootAtReferenceHeight + &
            log((RootAtReferenceHeight+1.0_kind_noahmp) / &
                (RootAtMomentumRoughness+1.0_kind_noahmp))
       Temporary = MoninObukhovParameter*2.0_kind_noahmp*InverseReferenceHeight
       RootAtReferenceHeight = sqrt(1.0_kind_noahmp + &
            4.0_kind_noahmp*StableAlpha*Temporary)
       HeatCorrection2m = RootAtHeatRoughness-RootAtReferenceHeight + &
            log((RootAtReferenceHeight+1.0_kind_noahmp) / &
                (RootAtHeatRoughness+1.0_kind_noahmp))
    else
       Temporary = ReferenceHeight/MoninObukhovParameter
       if ( abs(Temporary) <= 50.0_kind_noahmp*MomentumRoughnessLength ) then
          MoninObukhovParameter = max(MinimumStability, &
               -ReferenceHeight/(50.0_kind_noahmp*MomentumRoughnessLength))
       endif
       if ( MoninObukhovParameter >= -0.5_kind_noahmp ) then
          MomentumCorrection = (A0Momentum+A1Momentum*MoninObukhovParameter) * &
               MoninObukhovParameter / (1.0_kind_noahmp + &
               (B1Momentum+B2Momentum*MoninObukhovParameter)*MoninObukhovParameter)
          HeatCorrection = (A0Heat+A1Heat*MoninObukhovParameter) * &
               MoninObukhovParameter / (1.0_kind_noahmp + &
               (B1Heat+B2Heat*MoninObukhovParameter)*MoninObukhovParameter)
          Temporary = MoninObukhovParameter*10.0_kind_noahmp*InverseReferenceHeight
          MomentumCorrection10m = (A0Momentum+A1Momentum*Temporary)*Temporary / &
               (1.0_kind_noahmp+(B1Momentum+B2Momentum*Temporary)*Temporary)
          Temporary = MoninObukhovParameter*2.0_kind_noahmp*InverseReferenceHeight
          HeatCorrection2m = (A0Heat+A1Heat*Temporary)*Temporary / &
               (1.0_kind_noahmp+(B1Heat+B2Heat*Temporary)*Temporary)
       else
          Temporary = -MoninObukhovParameter
          MomentumCorrection = log(Temporary) + &
               2.0_kind_noahmp/sqrt(sqrt(Temporary)) - 0.8776_kind_noahmp
          HeatCorrection = log(Temporary) + &
               0.5_kind_noahmp/sqrt(Temporary) + 1.386_kind_noahmp
          Temporary = -MoninObukhovParameter*10.0_kind_noahmp*InverseReferenceHeight
          MomentumCorrection10m = log(Temporary) + &
               2.0_kind_noahmp/sqrt(sqrt(Temporary)) - 0.8776_kind_noahmp
          Temporary = -MoninObukhovParameter*2.0_kind_noahmp*InverseReferenceHeight
          HeatCorrection2m = log(Temporary) + &
               0.5_kind_noahmp/sqrt(Temporary) + 1.386_kind_noahmp
       endif
    endif

    MomentumSimilarityDenominator = MomentumSimilarityDenominator - MomentumCorrection
    HeatSimilarityDenominator = HeatSimilarityDenominator - HeatCorrection
    MomentumSimilarityDenominator10m = MomentumSimilarityDenominator10m - &
         MomentumCorrection10m
    HeatSimilarityDenominator2m = HeatSimilarityDenominator2m - HeatCorrection2m

    DimensionlessMomentumCoefficient = max(1.0e-5_kind_noahmp/ReferenceHeight, &
         VonKarmanConstant**2/MomentumSimilarityDenominator**2)
    DimensionlessHeatCoefficient = max(1.0e-5_kind_noahmp/ReferenceHeight, &
         VonKarmanConstant**2/(MomentumSimilarityDenominator*HeatSimilarityDenominator))
    MomentumConductance = DimensionlessMomentumCoefficient*WindSpeed
    HeatConductance = DimensionlessHeatCoefficient*WindSpeed
    SurfaceStress = DimensionlessMomentumCoefficient*WindSpeed**2
    FrictionVelocity = sqrt(SurfaceStress)

  end subroutine SurfaceExchangeGFS

end module SurfaceExchangeGFSMod
