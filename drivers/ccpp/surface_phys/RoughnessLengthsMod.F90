module RoughnessLengthsMod

! Driver-side thermal and moisture roughness-length formulations used by
! the CCPP MYNN surface-exchange option.
! Harsh Kamath (NOAA GSL/CIRES), September 2026


  use Machine, only : kind_noahmp

  implicit none
  private

  integer, parameter, public :: ScalarRoughnessZilitinkevich = 0
  integer, parameter, public :: ScalarRoughnessChenZhang     = 1
  integer, parameter, public :: ScalarRoughnessYang          = 2
  integer, parameter, public :: ScalarRoughnessGarratt       = 3

  public :: CalculateScalarRoughnessLengths

contains

  subroutine CalculateScalarRoughnessLengths( &
       MomentumRoughnessLength, KinematicViscosity, FrictionVelocity, &
       TemperatureScale, HumidityScale, SnowDepth, IsIceSurface, &
       ScalarRoughnessOption, HeatRoughnessLength, MoistureRoughnessLength)

    real(kind=kind_noahmp), intent(in)  :: MomentumRoughnessLength
    real(kind=kind_noahmp), intent(in)  :: KinematicViscosity
    real(kind=kind_noahmp), intent(in)  :: FrictionVelocity
    real(kind=kind_noahmp), intent(in)  :: TemperatureScale
    real(kind=kind_noahmp), intent(in)  :: HumidityScale
    real(kind=kind_noahmp), intent(in)  :: SnowDepth
    logical,                intent(in)  :: IsIceSurface
    integer,                intent(in)  :: ScalarRoughnessOption
    real(kind=kind_noahmp), intent(out) :: HeatRoughnessLength
    real(kind=kind_noahmp), intent(out) :: MoistureRoughnessLength

    real(kind=kind_noahmp) :: RoughnessReynoldsNumber

    RoughnessReynoldsNumber = max(0.1_kind_noahmp, &
         FrictionVelocity*MomentumRoughnessLength/KinematicViscosity)

    if ( SnowDepth > 0.05_kind_noahmp .or. IsIceSurface ) then
       call Andreas2002(MomentumRoughnessLength, KinematicViscosity, &
            FrictionVelocity, HeatRoughnessLength, MoistureRoughnessLength)
    elseif ( ScalarRoughnessOption == ScalarRoughnessYang ) then
       call Yang2008(MomentumRoughnessLength, KinematicViscosity, &
            FrictionVelocity, TemperatureScale, HumidityScale, &
            HeatRoughnessLength, MoistureRoughnessLength)
    elseif ( ScalarRoughnessOption == ScalarRoughnessGarratt ) then
       call Garratt1992(MomentumRoughnessLength, HeatRoughnessLength, &
            MoistureRoughnessLength)
    else
       call Zilitinkevich1995(MomentumRoughnessLength, &
            RoughnessReynoldsNumber, ScalarRoughnessOption, &
            HeatRoughnessLength, MoistureRoughnessLength)
    endif

    HeatRoughnessLength = max(2.0e-9_kind_noahmp, HeatRoughnessLength)
    MoistureRoughnessLength = max(2.0e-9_kind_noahmp, MoistureRoughnessLength)

  end subroutine CalculateScalarRoughnessLengths


  subroutine Zilitinkevich1995(MomentumRoughnessLength, &
       RoughnessReynoldsNumber, ScalarRoughnessOption, &
       HeatRoughnessLength, MoistureRoughnessLength)

    real(kind=kind_noahmp), intent(in)  :: MomentumRoughnessLength
    real(kind=kind_noahmp), intent(in)  :: RoughnessReynoldsNumber
    integer,                intent(in)  :: ScalarRoughnessOption
    real(kind=kind_noahmp), intent(out) :: HeatRoughnessLength
    real(kind=kind_noahmp), intent(out) :: MoistureRoughnessLength

    real(kind=kind_noahmp) :: ZilitinkevichCoefficient
    real(kind=kind_noahmp), parameter :: VonKarmanConstant = 0.4_kind_noahmp

    if ( ScalarRoughnessOption == ScalarRoughnessChenZhang ) then
       ZilitinkevichCoefficient = 10.0_kind_noahmp** &
            (-0.4_kind_noahmp*(MomentumRoughnessLength/0.07_kind_noahmp))
    else
       ZilitinkevichCoefficient = 0.085_kind_noahmp
    endif

    HeatRoughnessLength = MomentumRoughnessLength * exp( &
         -VonKarmanConstant*ZilitinkevichCoefficient*sqrt(RoughnessReynoldsNumber))
    HeatRoughnessLength = min(HeatRoughnessLength, &
         0.75_kind_noahmp*MomentumRoughnessLength)
    MoistureRoughnessLength = HeatRoughnessLength

  end subroutine Zilitinkevich1995


  subroutine Garratt1992(MomentumRoughnessLength, HeatRoughnessLength, &
       MoistureRoughnessLength)

    real(kind=kind_noahmp), intent(in)  :: MomentumRoughnessLength
    real(kind=kind_noahmp), intent(out) :: HeatRoughnessLength
    real(kind=kind_noahmp), intent(out) :: MoistureRoughnessLength

    HeatRoughnessLength = MomentumRoughnessLength / exp(2.0_kind_noahmp)
    MoistureRoughnessLength = HeatRoughnessLength

  end subroutine Garratt1992


  subroutine Yang2008(MomentumRoughnessLength, KinematicViscosity, &
       FrictionVelocity, TemperatureScale, HumidityScale, &
       HeatRoughnessLength, MoistureRoughnessLength)

    real(kind=kind_noahmp), intent(in)  :: MomentumRoughnessLength
    real(kind=kind_noahmp), intent(in)  :: KinematicViscosity
    real(kind=kind_noahmp), intent(in)  :: FrictionVelocity
    real(kind=kind_noahmp), intent(in)  :: TemperatureScale
    real(kind=kind_noahmp), intent(in)  :: HumidityScale
    real(kind=kind_noahmp), intent(out) :: HeatRoughnessLength
    real(kind=kind_noahmp), intent(out) :: MoistureRoughnessLength

    real(kind=kind_noahmp) :: BoundedMomentumRoughnessLength
    real(kind=kind_noahmp) :: CriticalRoughnessReynoldsNumber
    real(kind=kind_noahmp) :: CriticalRoughnessHeight
    real(kind=kind_noahmp), parameter :: Beta = 1.5_kind_noahmp

    BoundedMomentumRoughnessLength = min(0.5_kind_noahmp, max(0.04_kind_noahmp, MomentumRoughnessLength))
    CriticalRoughnessReynoldsNumber = 691.0_kind_noahmp + 170.0_kind_noahmp*log(BoundedMomentumRoughnessLength)
    CriticalRoughnessHeight = CriticalRoughnessReynoldsNumber * KinematicViscosity/max(FrictionVelocity, 0.01_kind_noahmp)

    HeatRoughnessLength = CriticalRoughnessHeight * exp(-Beta*sqrt(FrictionVelocity)*abs(min(TemperatureScale, 0.0_kind_noahmp)))
    MoistureRoughnessLength = CriticalRoughnessHeight * exp(-Beta*sqrt(FrictionVelocity)*abs(min(HumidityScale, 0.0_kind_noahmp)))
    HeatRoughnessLength = min(HeatRoughnessLength, 0.5_kind_noahmp*MomentumRoughnessLength)
    MoistureRoughnessLength = min(MoistureRoughnessLength, 0.5_kind_noahmp*MomentumRoughnessLength)

  end subroutine Yang2008


  subroutine Andreas2002(MomentumRoughnessLength, KinematicViscosity, &
       FrictionVelocity, HeatRoughnessLength, MoistureRoughnessLength)

    real(kind=kind_noahmp), intent(in)  :: MomentumRoughnessLength
    real(kind=kind_noahmp), intent(in)  :: KinematicViscosity
    real(kind=kind_noahmp), intent(in)  :: FrictionVelocity
    real(kind=kind_noahmp), intent(out) :: HeatRoughnessLength
    real(kind=kind_noahmp), intent(out) :: MoistureRoughnessLength

    real(kind=kind_noahmp) :: SnowMomentumRoughnessLength
    real(kind=kind_noahmp) :: SnowRoughnessReynoldsNumber
    real(kind=kind_noahmp) :: LogReynoldsNumber
    real(kind=kind_noahmp) :: B0Heat
    real(kind=kind_noahmp) :: B1Heat
    real(kind=kind_noahmp) :: B2Heat
    real(kind=kind_noahmp) :: B0Moisture
    real(kind=kind_noahmp) :: B1Moisture
    real(kind=kind_noahmp) :: B2Moisture

    SnowMomentumRoughnessLength = 0.135_kind_noahmp*KinematicViscosity / max(FrictionVelocity, 0.005_kind_noahmp) + &
         0.035_kind_noahmp*FrictionVelocity**2/9.8_kind_noahmp * &
         (5.0_kind_noahmp*exp(-((FrictionVelocity-0.18_kind_noahmp)/ &
         0.1_kind_noahmp)**2) + 1.0_kind_noahmp)
    SnowRoughnessReynoldsNumber = min(1000.0_kind_noahmp, max(1.0e-12_kind_noahmp, &
         FrictionVelocity*SnowMomentumRoughnessLength/KinematicViscosity))

    if ( SnowRoughnessReynoldsNumber <= 0.135_kind_noahmp ) then
       B0Heat = 1.25_kind_noahmp
       B1Heat = 0.0_kind_noahmp
       B2Heat = 0.0_kind_noahmp
       B0Moisture = 1.61_kind_noahmp
       B1Moisture = 0.0_kind_noahmp
       B2Moisture = 0.0_kind_noahmp
    elseif ( SnowRoughnessReynoldsNumber < 2.5_kind_noahmp ) then
       B0Heat = 0.149_kind_noahmp
       B1Heat = -0.55_kind_noahmp
       B2Heat = 0.0_kind_noahmp
       B0Moisture = 0.351_kind_noahmp
       B1Moisture = -0.628_kind_noahmp
       B2Moisture = 0.0_kind_noahmp
    else
       B0Heat = 0.317_kind_noahmp
       B1Heat = -0.565_kind_noahmp
       B2Heat = -0.183_kind_noahmp
       B0Moisture = 0.396_kind_noahmp
       B1Moisture = -0.512_kind_noahmp
       B2Moisture = -0.180_kind_noahmp
    endif

    LogReynoldsNumber = log(SnowRoughnessReynoldsNumber)
    HeatRoughnessLength = SnowMomentumRoughnessLength * exp( &
         B0Heat + B1Heat*LogReynoldsNumber + B2Heat*LogReynoldsNumber**2)
    MoistureRoughnessLength = SnowMomentumRoughnessLength * exp( &
         B0Moisture + B1Moisture*LogReynoldsNumber + B2Moisture*LogReynoldsNumber**2)

    if ( MomentumRoughnessLength <= 0.0_kind_noahmp ) then
       HeatRoughnessLength = 2.0e-9_kind_noahmp
       MoistureRoughnessLength = 2.0e-9_kind_noahmp
    endif

  end subroutine Andreas2002

end module RoughnessLengthsMod
