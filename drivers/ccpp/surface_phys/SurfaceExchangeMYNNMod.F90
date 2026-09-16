module SurfaceExchangeMYNNMod

! MYNN surface-layer exchange coefficients for the modular Noah-MP CCPP
! driver.  This is the scalar, column form of the legacy sfcdif4 routine.
! Harsh Kamath (NOAA GSL/CIRES), September 2026


  use Machine, only : kind_noahmp
  use NoahmpVarType, only : noahmp_type
  use RoughnessLengthsMod, only : CalculateScalarRoughnessLengths
  use StabilityFunctionsMod, only : MomentumStabilityFunction, &
       HeatStabilityFunction, MoninObukhovParameterFromBulkRichardson, &
       ApproximateMoninObukhovParameter

  implicit none
  private

  public :: SurfaceExchangeMYNN
  public :: UpdateNoahmpSurfaceExchangeMYNN

contains

  subroutine UpdateNoahmpSurfaceExchangeMYNN(noahmp, SurfaceTemperature, &
       SensibleHeatFlux, MoistureFlux, IsVegetated)

    type(noahmp_type),      intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: SurfaceTemperature
    real(kind=kind_noahmp), intent(in)    :: SensibleHeatFlux
    real(kind=kind_noahmp), intent(in)    :: MoistureFlux
    logical,                intent(in)    :: IsVegetated

    real(kind=kind_noahmp) :: MomentumRoughnessLength
    real(kind=kind_noahmp) :: HeatRoughnessLength
    real(kind=kind_noahmp) :: MoistureRoughnessLength
    real(kind=kind_noahmp) :: SurfaceSpecificHumidity
    real(kind=kind_noahmp) :: FrictionVelocity
    real(kind=kind_noahmp) :: BulkRichardsonNumber
    real(kind=kind_noahmp) :: ReciprocalMoninObukhovLength
    real(kind=kind_noahmp) :: MomentumDenominator
    real(kind=kind_noahmp) :: HeatDenominator
    real(kind=kind_noahmp) :: MomentumDenominator10m
    real(kind=kind_noahmp) :: HeatDenominator2m
    real(kind=kind_noahmp) :: MomentumConductance
    real(kind=kind_noahmp) :: HeatConductance
    real(kind=kind_noahmp) :: HeatConductance2m
    real(kind=kind_noahmp) :: MoistureConductance2m
    real(kind=kind_noahmp) :: SurfaceStress
    real(kind=kind_noahmp) :: EffectiveWindSpeed
    real(kind=kind_noahmp) :: HeatFluxCoefficient
    real(kind=kind_noahmp) :: MoistureFluxCoefficient
    logical                :: IsIceSurface

    if ( IsVegetated ) then
       MomentumRoughnessLength = noahmp%energy%state%RoughLenMomSfc
       FrictionVelocity = noahmp%energy%state%FrictionVelVeg
    else
       MomentumRoughnessLength = noahmp%energy%state%RoughLenMomGrd
       FrictionVelocity = noahmp%energy%state%FrictionVelBare
    endif
    SurfaceSpecificHumidity = noahmp%energy%state%SpecHumiditySfc
    IsIceSurface = noahmp%config%domain%IndicatorIceSfc /= 0

    call SurfaceExchangeMYNN(noahmp%config%domain%TimeStepIndex, &
         noahmp%forcing%WindEastwardRefHeight, &
         noahmp%forcing%WindNorthwardRefHeight, &
         noahmp%forcing%TemperatureAirRefHeight, &
         noahmp%forcing%PressureAirRefHeight, &
         noahmp%forcing%PressureAirSurface, &
         noahmp%forcing%BoundaryLayerHeight, noahmp%config%domain%GridSize, &
         MomentumRoughnessLength, noahmp%water%state%SnowDepth, IsIceSurface, &
         noahmp%config%nmlist%OptSurfaceStabilityFunction, SurfaceTemperature, &
         noahmp%forcing%SpecHumidityRefHeight, &
         noahmp%config%domain%RefHeightAboveSfc, &
         noahmp%config%nmlist%OptSurfaceThermalRoughness, &
         SurfaceSpecificHumidity, SensibleHeatFlux, MoistureFlux, &
         FrictionVelocity, BulkRichardsonNumber, ReciprocalMoninObukhovLength, &
         MomentumDenominator, HeatDenominator, MomentumDenominator10m, &
         HeatDenominator2m, MomentumConductance, HeatConductance, &
         HeatConductance2m, MoistureConductance2m, SurfaceStress, &
         EffectiveWindSpeed, HeatFluxCoefficient, MoistureFluxCoefficient, &
         HeatRoughnessLength, MoistureRoughnessLength)

    noahmp%energy%state%SpecHumiditySfc = SurfaceSpecificHumidity
    if ( IsVegetated ) then
       noahmp%energy%state%RoughLenShCanopy = HeatRoughnessLength
       noahmp%energy%state%ExchCoeffMomAbvCan = MomentumConductance/EffectiveWindSpeed
       noahmp%energy%state%ExchCoeffShAbvCan = HeatConductance/EffectiveWindSpeed
       noahmp%energy%state%ResistanceMomAbvCan = max(1.0_kind_noahmp, &
            1.0_kind_noahmp/MomentumConductance)
       noahmp%energy%state%ResistanceShAbvCan = max(1.0_kind_noahmp, &
            1.0_kind_noahmp/HeatConductance)
       noahmp%energy%state%ResistanceLhAbvCan = &
            noahmp%energy%state%ResistanceShAbvCan
       noahmp%energy%state%FrictionVelVeg = FrictionVelocity
       noahmp%energy%state%MoStabParaAbvCan = &
            ReciprocalMoninObukhovLength*noahmp%config%domain%RefHeightAboveSfc
       noahmp%energy%state%ExchCoeffSh2mVeg = HeatConductance2m
    else
       if ( noahmp%water%state%SnowDepth > 0.0_kind_noahmp ) then
          MomentumConductance = min(0.01_kind_noahmp*EffectiveWindSpeed, &
               MomentumConductance)
          HeatConductance = min(0.01_kind_noahmp*EffectiveWindSpeed, &
               HeatConductance)
          HeatConductance2m = min(0.01_kind_noahmp*EffectiveWindSpeed, &
               HeatConductance2m)
          MoistureConductance2m = min(0.01_kind_noahmp*EffectiveWindSpeed, &
               MoistureConductance2m)
       endif
       noahmp%energy%state%RoughLenShBareGrd = HeatRoughnessLength
       noahmp%energy%state%ExchCoeffMomBare = MomentumConductance/EffectiveWindSpeed
       noahmp%energy%state%ExchCoeffShBare = HeatConductance/EffectiveWindSpeed
       noahmp%energy%state%ResistanceMomBareGrd = max(1.0_kind_noahmp, &
            1.0_kind_noahmp/MomentumConductance)
       noahmp%energy%state%ResistanceShBareGrd = max(1.0_kind_noahmp, &
            1.0_kind_noahmp/HeatConductance)
       noahmp%energy%state%ResistanceLhBareGrd = &
            noahmp%energy%state%ResistanceShBareGrd
       noahmp%energy%state%FrictionVelBare = FrictionVelocity
       noahmp%energy%state%MoStabParaBare = &
            ReciprocalMoninObukhovLength*noahmp%config%domain%RefHeightAboveSfc
       noahmp%energy%state%ExchCoeffSh2mBare = HeatConductance2m
    endif

  end subroutine UpdateNoahmpSurfaceExchangeMYNN

  subroutine SurfaceExchangeMYNN(TimeStepIndex, WindEastward, WindNorthward, &
       AirTemperature, AirPressure, SurfaceAirPressure, BoundaryLayerHeight, &
       HorizontalGridLength, MomentumRoughnessLength, SnowDepth, IsIceSurface, &
       StabilityFunctionOption, SurfaceTemperature, AirSpecificHumidity, &
       ReferenceHeight, ScalarRoughnessOption, SurfaceSpecificHumidity, &
       SensibleHeatFlux, MoistureFlux, FrictionVelocity, BulkRichardsonNumber, &
       ReciprocalMoninObukhovLength, MomentumSimilarityDenominator, &
       HeatSimilarityDenominator, MomentumSimilarityDenominator10m, &
       HeatSimilarityDenominator2m, MomentumConductance, HeatConductance, &
       HeatConductance2m, MoistureConductance2m, SurfaceStress, &
       EffectiveWindSpeed, HeatFluxCoefficient, MoistureFluxCoefficient, &
       HeatRoughnessLength, MoistureRoughnessLength)

    integer,                intent(in)    :: TimeStepIndex
    real(kind=kind_noahmp), intent(in)    :: WindEastward
    real(kind=kind_noahmp), intent(in)    :: WindNorthward
    real(kind=kind_noahmp), intent(in)    :: AirTemperature
    real(kind=kind_noahmp), intent(in)    :: AirPressure
    real(kind=kind_noahmp), intent(in)    :: SurfaceAirPressure
    real(kind=kind_noahmp), intent(in)    :: BoundaryLayerHeight
    real(kind=kind_noahmp), intent(in)    :: HorizontalGridLength
    real(kind=kind_noahmp), intent(in)    :: MomentumRoughnessLength
    real(kind=kind_noahmp), intent(in)    :: SnowDepth
    logical,                intent(in)    :: IsIceSurface
    integer,                intent(in)    :: StabilityFunctionOption
    real(kind=kind_noahmp), intent(in)    :: SurfaceTemperature
    real(kind=kind_noahmp), intent(in)    :: AirSpecificHumidity
    real(kind=kind_noahmp), intent(in)    :: ReferenceHeight
    integer,                intent(in)    :: ScalarRoughnessOption
    real(kind=kind_noahmp), intent(inout) :: SurfaceSpecificHumidity
    real(kind=kind_noahmp), intent(in)    :: SensibleHeatFlux
    real(kind=kind_noahmp), intent(in)    :: MoistureFlux
    real(kind=kind_noahmp), intent(inout) :: FrictionVelocity
    real(kind=kind_noahmp), intent(out)   :: BulkRichardsonNumber
    real(kind=kind_noahmp), intent(out)   :: ReciprocalMoninObukhovLength
    real(kind=kind_noahmp), intent(out)   :: MomentumSimilarityDenominator
    real(kind=kind_noahmp), intent(out)   :: HeatSimilarityDenominator
    real(kind=kind_noahmp), intent(out)   :: MomentumSimilarityDenominator10m
    real(kind=kind_noahmp), intent(out)   :: HeatSimilarityDenominator2m
    real(kind=kind_noahmp), intent(out)   :: MomentumConductance
    real(kind=kind_noahmp), intent(out)   :: HeatConductance
    real(kind=kind_noahmp), intent(out)   :: HeatConductance2m
    real(kind=kind_noahmp), intent(out)   :: MoistureConductance2m
    real(kind=kind_noahmp), intent(out)   :: SurfaceStress
    real(kind=kind_noahmp), intent(out)   :: EffectiveWindSpeed
    real(kind=kind_noahmp), intent(out)   :: HeatFluxCoefficient
    real(kind=kind_noahmp), intent(out)   :: MoistureFluxCoefficient
    real(kind=kind_noahmp), intent(out)   :: HeatRoughnessLength
    real(kind=kind_noahmp), intent(out)   :: MoistureRoughnessLength

    real(kind=kind_noahmp), parameter :: VonKarmanConstant = 0.4_kind_noahmp
    real(kind=kind_noahmp), parameter :: GasConstantDryAir = 287.04_kind_noahmp
    real(kind=kind_noahmp), parameter :: HeatCapacityDryAir = 1004.64_kind_noahmp
    real(kind=kind_noahmp), parameter :: GravityAcceleration = 9.80616_kind_noahmp
    real(kind=kind_noahmp), parameter :: ReferencePressure = 100000.0_kind_noahmp
    real(kind=kind_noahmp), parameter :: Epsilon = 0.622_kind_noahmp
    real(kind=kind_noahmp), parameter :: VirtualTemperatureFactor = &
         1.0_kind_noahmp/Epsilon - 1.0_kind_noahmp
    real(kind=kind_noahmp), parameter :: MinimumWindSpeed = 0.1_kind_noahmp
    real(kind=kind_noahmp), parameter :: ConvectiveVelocityCoefficient = 1.25_kind_noahmp

    real(kind=kind_noahmp) :: AirPotentialTemperature
    real(kind=kind_noahmp) :: VirtualAirPotentialTemperature
    real(kind=kind_noahmp) :: VirtualAirTemperature
    real(kind=kind_noahmp) :: VirtualSurfacePotentialTemperature
    real(kind=kind_noahmp) :: SurfacePotentialTemperature
    real(kind=kind_noahmp) :: AirDensity
    real(kind=kind_noahmp) :: MoistAirHeatCapacity
    real(kind=kind_noahmp) :: ConvectiveVelocity
    real(kind=kind_noahmp) :: SubgridVelocity
    real(kind=kind_noahmp) :: BuoyancyFlux
    real(kind=kind_noahmp) :: KinematicViscosity
    real(kind=kind_noahmp) :: TemperatureScale
    real(kind=kind_noahmp) :: HumidityScale
    real(kind=kind_noahmp) :: InitialMoninObukhovParameter
    real(kind=kind_noahmp) :: MoninObukhovParameter
    real(kind=kind_noahmp) :: ParameterAtMomentumRoughness
    real(kind=kind_noahmp) :: ParameterAtHeatRoughness
    real(kind=kind_noahmp) :: ParameterAtReferenceHeight
    real(kind=kind_noahmp) :: ParameterAt2m
    real(kind=kind_noahmp) :: ParameterAt10m
    real(kind=kind_noahmp) :: MomentumCorrection
    real(kind=kind_noahmp) :: HeatCorrection
    real(kind=kind_noahmp) :: MomentumCorrection10m
    real(kind=kind_noahmp) :: HeatCorrection2m
    real(kind=kind_noahmp) :: MomentumNeutralDenominator
    real(kind=kind_noahmp) :: HeatNeutralDenominator
    real(kind=kind_noahmp) :: MomentumNeutralDenominator10m
    real(kind=kind_noahmp) :: HeatNeutralDenominator2m
    real(kind=kind_noahmp) :: MoistureDenominator
    real(kind=kind_noahmp) :: MoistureDenominator2m
    real(kind=kind_noahmp) :: SaturationVaporPressure
    real(kind=kind_noahmp) :: SaturationTemperature
    real(kind=kind_noahmp) :: SurfacePressureKPa

    SurfacePressureKPa = SurfaceAirPressure/1000.0_kind_noahmp
    SaturationTemperature = SurfaceTemperature
    if ( .not. IsIceSurface ) &
       SaturationTemperature = 0.5_kind_noahmp*(SurfaceTemperature+AirTemperature)
    call SaturationSpecificHumidity(SaturationTemperature, SurfacePressureKPa, &
         SurfaceSpecificHumidity)

    SurfacePotentialTemperature = SurfaceTemperature * &
         (ReferencePressure/SurfaceAirPressure)**(GasConstantDryAir/HeatCapacityDryAir)
    AirPotentialTemperature = AirTemperature * &
         (ReferencePressure/AirPressure)**(GasConstantDryAir/HeatCapacityDryAir)
    VirtualAirPotentialTemperature = AirPotentialTemperature * &
         (1.0_kind_noahmp+VirtualTemperatureFactor*AirSpecificHumidity)
    VirtualAirTemperature = AirTemperature * &
         (1.0_kind_noahmp+VirtualTemperatureFactor*AirSpecificHumidity)
    VirtualSurfacePotentialTemperature = SurfacePotentialTemperature * &
         (1.0_kind_noahmp+VirtualTemperatureFactor*SurfaceSpecificHumidity)

    AirDensity = SurfaceAirPressure/(GasConstantDryAir*VirtualAirTemperature)
    MoistAirHeatCapacity = HeatCapacityDryAir * &
         (1.0_kind_noahmp+0.84_kind_noahmp*AirSpecificHumidity / &
          max(1.0e-8_kind_noahmp, 1.0_kind_noahmp-AirSpecificHumidity))

    BuoyancyFlux = max(SensibleHeatFlux/(AirDensity*HeatCapacityDryAir) + &
         VirtualTemperatureFactor*VirtualSurfacePotentialTemperature* &
         MoistureFlux/AirDensity, 0.0_kind_noahmp)
    ConvectiveVelocity = ConvectiveVelocityCoefficient * &
         (GravityAcceleration/SurfaceTemperature * &
          min(1.5_kind_noahmp*BoundaryLayerHeight, 4000.0_kind_noahmp) * &
          BuoyancyFlux)**(1.0_kind_noahmp/3.0_kind_noahmp)
    SubgridVelocity = min(0.32_kind_noahmp * &
         max(HorizontalGridLength/5000.0_kind_noahmp-1.0_kind_noahmp, &
             0.0_kind_noahmp)**(1.0_kind_noahmp/3.0_kind_noahmp), &
         0.5_kind_noahmp)
    EffectiveWindSpeed = max(MinimumWindSpeed, sqrt(WindEastward**2 + &
         WindNorthward**2 + ConvectiveVelocity**2 + SubgridVelocity**2))

    BulkRichardsonNumber = GravityAcceleration/AirPotentialTemperature * &
         ReferenceHeight*(VirtualAirPotentialTemperature - &
         VirtualSurfacePotentialTemperature)/EffectiveWindSpeed**2
    if ( TimeStepIndex == 1 ) then
       BulkRichardsonNumber = min(2.0_kind_noahmp, &
            max(-2.0_kind_noahmp, BulkRichardsonNumber))
    else
       BulkRichardsonNumber = min(4.0_kind_noahmp, &
            max(-4.0_kind_noahmp, BulkRichardsonNumber))
    endif

    KinematicViscosity = 1.326e-5_kind_noahmp * &
         (1.0_kind_noahmp + 6.542e-3_kind_noahmp*(AirTemperature-273.15_kind_noahmp) + &
          8.301e-6_kind_noahmp*(AirTemperature-273.15_kind_noahmp)**2 - &
          4.84e-9_kind_noahmp*(AirTemperature-273.15_kind_noahmp)**3)

    TemperatureScale = 0.0_kind_noahmp
    HumidityScale = 0.0_kind_noahmp
    call CalculateScalarRoughnessLengths(MomentumRoughnessLength, &
         KinematicViscosity, FrictionVelocity, TemperatureScale, HumidityScale, &
         SnowDepth, IsIceSurface, ScalarRoughnessOption, HeatRoughnessLength, &
         MoistureRoughnessLength)

    MomentumNeutralDenominator = log((ReferenceHeight+MomentumRoughnessLength) / &
         MomentumRoughnessLength)
    HeatNeutralDenominator = log((ReferenceHeight+MomentumRoughnessLength) / &
         HeatRoughnessLength)
    MomentumNeutralDenominator10m = log((10.0_kind_noahmp+MomentumRoughnessLength) / &
         MomentumRoughnessLength)
    HeatNeutralDenominator2m = log((2.0_kind_noahmp+MomentumRoughnessLength) / &
         HeatRoughnessLength)

    call ApproximateMoninObukhovParameter(InitialMoninObukhovParameter, &
         BulkRichardsonNumber, &
         ReferenceHeight/MomentumRoughnessLength, &
         MomentumRoughnessLength/HeatRoughnessLength)
    if ( BulkRichardsonNumber > 0.0_kind_noahmp ) then
       InitialMoninObukhovParameter = min(20.0_kind_noahmp, &
            max(0.0_kind_noahmp, InitialMoninObukhovParameter))
    else
       InitialMoninObukhovParameter = min(0.0_kind_noahmp, &
            max(-20.0_kind_noahmp, InitialMoninObukhovParameter))
    endif
    MoninObukhovParameter = MoninObukhovParameterFromBulkRichardson( &
         BulkRichardsonNumber, ReferenceHeight, MomentumRoughnessLength, &
         HeatRoughnessLength, InitialMoninObukhovParameter, &
         StabilityFunctionOption)
    MoninObukhovParameter = min(20.0_kind_noahmp, &
         max(-20.0_kind_noahmp, MoninObukhovParameter))

    ParameterAtMomentumRoughness = MoninObukhovParameter * &
         MomentumRoughnessLength/ReferenceHeight
    ParameterAtHeatRoughness = MoninObukhovParameter * &
         HeatRoughnessLength/ReferenceHeight
    ParameterAtReferenceHeight = MoninObukhovParameter * &
         (ReferenceHeight+MomentumRoughnessLength)/ReferenceHeight
    ParameterAt10m = MoninObukhovParameter * &
         (10.0_kind_noahmp+MomentumRoughnessLength)/ReferenceHeight
    ParameterAt2m = MoninObukhovParameter * &
         (2.0_kind_noahmp+MomentumRoughnessLength)/ReferenceHeight

    MomentumCorrection = MomentumStabilityFunction(ParameterAtReferenceHeight, &
         StabilityFunctionOption) - MomentumStabilityFunction( &
         ParameterAtMomentumRoughness, StabilityFunctionOption)
    HeatCorrection = HeatStabilityFunction(ParameterAtReferenceHeight, &
         StabilityFunctionOption) - HeatStabilityFunction( &
         ParameterAtHeatRoughness, StabilityFunctionOption)
    MomentumCorrection10m = MomentumStabilityFunction(ParameterAt10m, &
         StabilityFunctionOption) - MomentumStabilityFunction( &
         ParameterAtMomentumRoughness, StabilityFunctionOption)
    HeatCorrection2m = HeatStabilityFunction(ParameterAt2m, &
         StabilityFunctionOption) - HeatStabilityFunction( &
         ParameterAtHeatRoughness, StabilityFunctionOption)

    if ( BulkRichardsonNumber < 0.0_kind_noahmp ) then
       MomentumCorrection = min(MomentumCorrection, &
            0.9_kind_noahmp*MomentumNeutralDenominator)
       HeatCorrection = min(HeatCorrection, &
            0.9_kind_noahmp*HeatNeutralDenominator)
       MomentumCorrection10m = min(MomentumCorrection10m, &
            0.9_kind_noahmp*MomentumNeutralDenominator10m)
       HeatCorrection2m = min(HeatCorrection2m, &
            0.9_kind_noahmp*HeatNeutralDenominator2m)
    endif

    MomentumSimilarityDenominator = max(1.0_kind_noahmp, &
         MomentumNeutralDenominator-MomentumCorrection)
    HeatSimilarityDenominator = max(1.0_kind_noahmp, &
         HeatNeutralDenominator-HeatCorrection)
    MomentumSimilarityDenominator10m = max(1.0_kind_noahmp, &
         MomentumNeutralDenominator10m-MomentumCorrection10m)
    HeatSimilarityDenominator2m = max(1.0_kind_noahmp, &
         HeatNeutralDenominator2m-HeatCorrection2m)
    MoistureDenominator = max(1.0_kind_noahmp, &
         log((ReferenceHeight+MoistureRoughnessLength)/MoistureRoughnessLength) - &
         HeatCorrection)
    MoistureDenominator2m = max(1.0_kind_noahmp, &
         log((2.0_kind_noahmp+MoistureRoughnessLength)/MoistureRoughnessLength) - &
         HeatCorrection2m)

    FrictionVelocity = max(0.005_kind_noahmp, &
         0.5_kind_noahmp*FrictionVelocity + 0.5_kind_noahmp* &
         VonKarmanConstant*EffectiveWindSpeed/MomentumSimilarityDenominator)
    TemperatureScale = VonKarmanConstant * &
         (VirtualAirPotentialTemperature-VirtualSurfacePotentialTemperature) / &
         HeatSimilarityDenominator
    HumidityScale = VonKarmanConstant * &
         (AirSpecificHumidity-SurfaceSpecificHumidity)*1000.0_kind_noahmp / &
         MoistureDenominator

    MomentumConductance = (VonKarmanConstant/MomentumSimilarityDenominator)**2 * &
         EffectiveWindSpeed
    HeatConductance = FrictionVelocity*VonKarmanConstant/HeatSimilarityDenominator
    HeatConductance2m = FrictionVelocity*VonKarmanConstant/HeatSimilarityDenominator2m
    MoistureConductance2m = FrictionVelocity*VonKarmanConstant/MoistureDenominator2m
    SurfaceStress = FrictionVelocity**2
    ReciprocalMoninObukhovLength = MoninObukhovParameter/ReferenceHeight
    HeatFluxCoefficient = AirDensity*MoistAirHeatCapacity*HeatConductance
    MoistureFluxCoefficient = AirDensity*FrictionVelocity*VonKarmanConstant / &
         MoistureDenominator

  end subroutine SurfaceExchangeMYNN


  subroutine SaturationSpecificHumidity(Temperature, PressureKPa, SpecificHumidity)

    real(kind=kind_noahmp), intent(in)  :: Temperature
    real(kind=kind_noahmp), intent(in)  :: PressureKPa
    real(kind=kind_noahmp), intent(out) :: SpecificHumidity

    real(kind=kind_noahmp), parameter :: Epsilon = 0.622_kind_noahmp
    real(kind=kind_noahmp) :: SaturationVaporPressure

    if ( Temperature < 273.15_kind_noahmp ) then
       SaturationVaporPressure = 0.6112_kind_noahmp * exp( &
            4648.0_kind_noahmp*(1.0_kind_noahmp/273.15_kind_noahmp - &
            1.0_kind_noahmp/Temperature) - &
            11.64_kind_noahmp*log(273.15_kind_noahmp/Temperature) + &
            0.02265_kind_noahmp*(273.15_kind_noahmp-Temperature))
    else
       SaturationVaporPressure = 0.6112_kind_noahmp * exp( &
            17.67_kind_noahmp*(Temperature-273.15_kind_noahmp) / &
            (Temperature-29.65_kind_noahmp))
    endif
    SpecificHumidity = Epsilon*SaturationVaporPressure / &
         (PressureKPa-(1.0_kind_noahmp-Epsilon)*SaturationVaporPressure)

  end subroutine SaturationSpecificHumidity

end module SurfaceExchangeMYNNMod
