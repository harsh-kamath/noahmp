module StabilityFunctionsMod

! Integrated surface-layer stability functions used by the CCPP MYNN
! surface-exchange option.  The misspelling in this source filename is
! retained to match the requested driver layout; the module name is spelled
! correctly.

! Harsh Kamath (NOAA GSL/CIRES), September 2026


  use Machine, only : kind_noahmp

  implicit none
  private

  integer, parameter, public :: StabilityFunctionMynn = 0
  integer, parameter, public :: StabilityFunctionGfs  = 1

  ! Initialized during the serialized CCPP init phase, then read-only.
  real(kind=kind_noahmp), save :: MomentumStableTable( &
       0:1000, StabilityFunctionMynn:StabilityFunctionGfs)
  real(kind=kind_noahmp), save :: MomentumUnstableTable( &
       0:1000, StabilityFunctionMynn:StabilityFunctionGfs)
  real(kind=kind_noahmp), save :: HeatStableTable( &
       0:1000, StabilityFunctionMynn:StabilityFunctionGfs)
  real(kind=kind_noahmp), save :: HeatUnstableTable( &
       0:1000, StabilityFunctionMynn:StabilityFunctionGfs)
  logical, save :: StabilityTablesInitialized = .false.

  public :: InitializeStabilityFunctions
  public :: MomentumStabilityFunction
  public :: HeatStabilityFunction
  public :: MoninObukhovParameterFromBulkRichardson
  public :: ApproximateMoninObukhovParameter

contains

  subroutine InitializeStabilityFunctions(StabilityFunctionOption, ErrMsg, ErrFlg)

    integer,          intent(in)  :: StabilityFunctionOption
    character(len=*), intent(out) :: ErrMsg
    integer,          intent(out) :: ErrFlg

    integer :: TableIndex
    real(kind=kind_noahmp) :: MoninObukhovParameter

    ErrMsg = ''
    ErrFlg = 0

    if ( StabilityFunctionOption /= StabilityFunctionMynn .and. &
         StabilityFunctionOption /= StabilityFunctionGfs ) then
       ErrFlg = 1
       write(ErrMsg, '(a,i0)') &
            'InitializeStabilityFunctions: unsupported stability-function option ', &
            StabilityFunctionOption
       return
    endif

    if ( StabilityTablesInitialized ) return

    do TableIndex = 0, 1000
       MoninObukhovParameter = 0.01_kind_noahmp*real(TableIndex, kind_noahmp)
       MomentumStableTable(TableIndex, StabilityFunctionMynn) = &
            MynnMomentumStable(MoninObukhovParameter)
       HeatStableTable(TableIndex, StabilityFunctionMynn) = &
            MynnHeatStable(MoninObukhovParameter)
       MomentumUnstableTable(TableIndex, StabilityFunctionMynn) = &
            MynnMomentumUnstable(-MoninObukhovParameter)
       HeatUnstableTable(TableIndex, StabilityFunctionMynn) = &
            MynnHeatUnstable(-MoninObukhovParameter)
       MomentumStableTable(TableIndex, StabilityFunctionGfs) = &
            GfsMomentumStable(MoninObukhovParameter)
       HeatStableTable(TableIndex, StabilityFunctionGfs) = &
            GfsHeatStable(MoninObukhovParameter)
       MomentumUnstableTable(TableIndex, StabilityFunctionGfs) = &
            GfsMomentumUnstable(-MoninObukhovParameter)
       HeatUnstableTable(TableIndex, StabilityFunctionGfs) = &
            GfsHeatUnstable(-MoninObukhovParameter)
    enddo
    StabilityTablesInitialized = .true.

  end subroutine InitializeStabilityFunctions


  real(kind=kind_noahmp) function MomentumStabilityFunction( &
       MoninObukhovParameter, StabilityFunctionOption) result(StabilityCorrection)

    real(kind=kind_noahmp), intent(in) :: MoninObukhovParameter
    integer,                intent(in) :: StabilityFunctionOption

    integer :: TableIndex
    real(kind=kind_noahmp) :: TableFraction

    TableIndex = int(abs(MoninObukhovParameter)*100.0_kind_noahmp)
    TableFraction = abs(MoninObukhovParameter)*100.0_kind_noahmp - TableIndex
    if ( StabilityTablesInitialized .and. TableIndex+1 < 1000 ) then
       if ( MoninObukhovParameter >= 0.0_kind_noahmp ) then
          StabilityCorrection = MomentumStableTable(TableIndex, StabilityFunctionOption) + &
               TableFraction * (MomentumStableTable(TableIndex+1, StabilityFunctionOption) - &
               MomentumStableTable(TableIndex, StabilityFunctionOption))
       else
          StabilityCorrection = MomentumUnstableTable(TableIndex, StabilityFunctionOption) + &
               TableFraction * (MomentumUnstableTable(TableIndex+1, StabilityFunctionOption) - &
               MomentumUnstableTable(TableIndex, StabilityFunctionOption))
       endif
       return
    endif

    if ( MoninObukhovParameter >= 0.0_kind_noahmp ) then
       if ( StabilityFunctionOption == StabilityFunctionMynn ) then
          StabilityCorrection = MynnMomentumStable(MoninObukhovParameter)
       else
          StabilityCorrection = GfsMomentumStable(MoninObukhovParameter)
       endif
    else
       if ( StabilityFunctionOption == StabilityFunctionMynn ) then
          StabilityCorrection = MynnMomentumUnstable(MoninObukhovParameter)
       else
          StabilityCorrection = GfsMomentumUnstable(MoninObukhovParameter)
       endif
    endif

  end function MomentumStabilityFunction


  real(kind=kind_noahmp) function HeatStabilityFunction( &
       MoninObukhovParameter, StabilityFunctionOption) result(StabilityCorrection)

    real(kind=kind_noahmp), intent(in) :: MoninObukhovParameter
    integer,                intent(in) :: StabilityFunctionOption

    integer :: TableIndex
    real(kind=kind_noahmp) :: TableFraction

    TableIndex = int(abs(MoninObukhovParameter)*100.0_kind_noahmp)
    TableFraction = abs(MoninObukhovParameter)*100.0_kind_noahmp - TableIndex
    if ( StabilityTablesInitialized .and. TableIndex+1 < 1000 ) then
       if ( MoninObukhovParameter >= 0.0_kind_noahmp ) then
          StabilityCorrection = HeatStableTable(TableIndex, StabilityFunctionOption) + &
               TableFraction * (HeatStableTable(TableIndex+1, StabilityFunctionOption) - &
               HeatStableTable(TableIndex, StabilityFunctionOption))
       else
          StabilityCorrection = HeatUnstableTable(TableIndex, StabilityFunctionOption) + &
               TableFraction * (HeatUnstableTable(TableIndex+1, StabilityFunctionOption) - &
               HeatUnstableTable(TableIndex, StabilityFunctionOption))
       endif
       return
    endif

    if ( MoninObukhovParameter >= 0.0_kind_noahmp ) then
       if ( StabilityFunctionOption == StabilityFunctionMynn ) then
          StabilityCorrection = MynnHeatStable(MoninObukhovParameter)
       else
          StabilityCorrection = GfsHeatStable(MoninObukhovParameter)
       endif
    else
       if ( StabilityFunctionOption == StabilityFunctionMynn ) then
          StabilityCorrection = MynnHeatUnstable(MoninObukhovParameter)
       else
          StabilityCorrection = GfsHeatUnstable(MoninObukhovParameter)
       endif
    endif

  end function HeatStabilityFunction


  real(kind=kind_noahmp) function MoninObukhovParameterFromBulkRichardson( &
       BulkRichardsonNumber, ReferenceHeight, RoughnessLengthMomentum, &
       RoughnessLengthHeat, InitialMoninObukhovParameter, &
       StabilityFunctionOption) result(MoninObukhovParameter)

    real(kind=kind_noahmp), intent(in) :: BulkRichardsonNumber
    real(kind=kind_noahmp), intent(in) :: ReferenceHeight
    real(kind=kind_noahmp), intent(in) :: RoughnessLengthMomentum
    real(kind=kind_noahmp), intent(in) :: RoughnessLengthHeat
    real(kind=kind_noahmp), intent(in) :: InitialMoninObukhovParameter
    integer,                intent(in) :: StabilityFunctionOption

    ! The legacy loop advances n=1 through n=19 before invoking its fitted
    ! fallback at n=20.
    integer, parameter :: MaximumIterations = 19
    integer            :: Iteration
    real(kind=kind_noahmp) :: PreviousParameter
    real(kind=kind_noahmp) :: MomentumDenominator
    real(kind=kind_noahmp) :: HeatDenominator
    real(kind=kind_noahmp) :: ParameterAtMomentumRoughness
    real(kind=kind_noahmp) :: ParameterAtHeatRoughness
    real(kind=kind_noahmp) :: ParameterAtReferenceHeight
    real(kind=kind_noahmp) :: NeutralMomentumDenominator
    real(kind=kind_noahmp) :: NeutralHeatDenominator

    NeutralMomentumDenominator = log((ReferenceHeight + RoughnessLengthMomentum) / &
                                      RoughnessLengthMomentum)
    NeutralHeatDenominator = log((ReferenceHeight + RoughnessLengthMomentum) / &
                                  RoughnessLengthHeat)

    MoninObukhovParameter = InitialMoninObukhovParameter
    if ( MoninObukhovParameter * BulkRichardsonNumber < 0.0_kind_noahmp ) &
       MoninObukhovParameter = 0.0_kind_noahmp

    do Iteration = 1, MaximumIterations
       PreviousParameter = MoninObukhovParameter
       ParameterAtMomentumRoughness = PreviousParameter * RoughnessLengthMomentum / ReferenceHeight
       ParameterAtHeatRoughness = PreviousParameter * RoughnessLengthHeat / ReferenceHeight
       ParameterAtReferenceHeight = PreviousParameter * &
            (ReferenceHeight + RoughnessLengthMomentum) / ReferenceHeight

       MomentumDenominator = max(1.0_kind_noahmp, NeutralMomentumDenominator - &
            (MomentumStabilityFunction(ParameterAtReferenceHeight, StabilityFunctionOption) - &
             MomentumStabilityFunction(ParameterAtMomentumRoughness, StabilityFunctionOption)))
       HeatDenominator = max(1.0_kind_noahmp, NeutralHeatDenominator - &
            (HeatStabilityFunction(ParameterAtReferenceHeight, StabilityFunctionOption) - &
             HeatStabilityFunction(ParameterAtHeatRoughness, StabilityFunctionOption)))

       MoninObukhovParameter = BulkRichardsonNumber * &
            MomentumDenominator**2 / HeatDenominator
       if ( abs(MoninObukhovParameter - PreviousParameter) <= 0.01_kind_noahmp ) return
    enddo

    ! The fitted Li et al. (2010) relation is the legacy MYNN fallback when
    ! the fixed-point iteration does not converge.
    call ApproximateMoninObukhovParameter(MoninObukhovParameter, BulkRichardsonNumber, &
         ReferenceHeight / RoughnessLengthMomentum, &
         RoughnessLengthMomentum / RoughnessLengthHeat)

  end function MoninObukhovParameterFromBulkRichardson


  subroutine ApproximateMoninObukhovParameter(MoninObukhovParameter, BulkRichardsonNumber, &
       HeightRoughnessRatio, RoughnessRatio)

    real(kind=kind_noahmp), intent(out) :: MoninObukhovParameter
    real(kind=kind_noahmp), intent(in)  :: BulkRichardsonNumber
    real(kind=kind_noahmp), intent(in)  :: HeightRoughnessRatio
    real(kind=kind_noahmp), intent(in)  :: RoughnessRatio

    real(kind=kind_noahmp) :: LogHeightRoughnessRatio
    real(kind=kind_noahmp) :: LogRoughnessRatio
    real(kind=kind_noahmp) :: BoundedHeightRoughnessRatio
    real(kind=kind_noahmp) :: BoundedRoughnessRatio

    real(kind=kind_noahmp), parameter :: Au11 = 0.045_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bu11 = 0.003_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bu12 = 0.0059_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bu21 = -0.0828_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bu22 = 0.8845_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bu31 = 0.1739_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bu32 = -0.9213_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bu33 = -0.1057_kind_noahmp
    real(kind=kind_noahmp), parameter :: Aw11 = 0.5738_kind_noahmp
    real(kind=kind_noahmp), parameter :: Aw12 = -0.4399_kind_noahmp
    real(kind=kind_noahmp), parameter :: Aw21 = -4.901_kind_noahmp
    real(kind=kind_noahmp), parameter :: Aw22 = 52.50_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bw11 = -0.0539_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bw12 = 1.540_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bw21 = -0.669_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bw22 = -3.282_kind_noahmp
    real(kind=kind_noahmp), parameter :: As11 = 0.7529_kind_noahmp
    real(kind=kind_noahmp), parameter :: As21 = 14.94_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bs11 = 0.1569_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bs21 = -0.3091_kind_noahmp
    real(kind=kind_noahmp), parameter :: Bs22 = -1.303_kind_noahmp

    BoundedHeightRoughnessRatio = min(100000.0_kind_noahmp, &
         max(100.0_kind_noahmp, HeightRoughnessRatio))
    BoundedRoughnessRatio = min(100.0_kind_noahmp, &
         max(0.5_kind_noahmp, RoughnessRatio))
    LogHeightRoughnessRatio = log(BoundedHeightRoughnessRatio)
    LogRoughnessRatio = log(BoundedRoughnessRatio)

    if ( BulkRichardsonNumber <= 0.0_kind_noahmp ) then
       MoninObukhovParameter = Au11 * LogHeightRoughnessRatio * BulkRichardsonNumber**2 + &
            ((Bu11*LogRoughnessRatio + Bu12)*LogHeightRoughnessRatio**2 + &
             (Bu21*LogRoughnessRatio + Bu22)*LogHeightRoughnessRatio + &
             (Bu31*LogRoughnessRatio**2 + Bu32*LogRoughnessRatio + Bu33)) * &
            BulkRichardsonNumber
       MoninObukhovParameter = min(0.0_kind_noahmp, &
            max(-15.0_kind_noahmp, MoninObukhovParameter))
    elseif ( BulkRichardsonNumber <= 0.2_kind_noahmp ) then
       MoninObukhovParameter = &
            ((Aw11*LogRoughnessRatio + Aw12)*LogHeightRoughnessRatio + &
             (Aw21*LogRoughnessRatio + Aw22))*BulkRichardsonNumber**2 + &
            ((Bw11*LogRoughnessRatio + Bw12)*LogHeightRoughnessRatio + &
             (Bw21*LogRoughnessRatio + Bw22))*BulkRichardsonNumber
       MoninObukhovParameter = min(4.0_kind_noahmp, &
            max(0.0_kind_noahmp, MoninObukhovParameter))
    else
       MoninObukhovParameter = (As11*LogHeightRoughnessRatio + As21) * &
            BulkRichardsonNumber + Bs11*LogHeightRoughnessRatio + &
            Bs21*LogRoughnessRatio + Bs22
       MoninObukhovParameter = min(20.0_kind_noahmp, &
            max(1.0_kind_noahmp, MoninObukhovParameter))
    endif

  end subroutine ApproximateMoninObukhovParameter


  pure real(kind=kind_noahmp) function MynnMomentumStable(Parameter) result(Correction)
    real(kind=kind_noahmp), intent(in) :: Parameter
    Correction = -6.1_kind_noahmp * log(Parameter + &
         (1.0_kind_noahmp + Parameter**2.5_kind_noahmp)**0.4_kind_noahmp)
  end function MynnMomentumStable


  pure real(kind=kind_noahmp) function MynnHeatStable(Parameter) result(Correction)
    real(kind=kind_noahmp), intent(in) :: Parameter
    Correction = -5.3_kind_noahmp * log(Parameter + &
         (1.0_kind_noahmp + Parameter**1.1_kind_noahmp)**(1.0_kind_noahmp/1.1_kind_noahmp))
  end function MynnHeatStable


  pure real(kind=kind_noahmp) function MynnMomentumUnstable(Parameter) result(Correction)
    real(kind=kind_noahmp), intent(in) :: Parameter
    real(kind=kind_noahmp) :: BusingerArgument
    real(kind=kind_noahmp) :: ConvectiveArgument
    real(kind=kind_noahmp) :: BusingerCorrection
    real(kind=kind_noahmp) :: ConvectiveCorrection
    real(kind=kind_noahmp), parameter :: SqrtThree = 1.7320508075688773_kind_noahmp
    real(kind=kind_noahmp), parameter :: ArcTangentOne = 0.7853981633974483_kind_noahmp

    BusingerArgument = (1.0_kind_noahmp - 16.0_kind_noahmp*Parameter)**0.25_kind_noahmp
    BusingerCorrection = 2.0_kind_noahmp*log(0.5_kind_noahmp*(1.0_kind_noahmp+BusingerArgument)) + &
         log(0.5_kind_noahmp*(1.0_kind_noahmp+BusingerArgument**2)) - &
         2.0_kind_noahmp*atan(BusingerArgument) + 2.0_kind_noahmp*ArcTangentOne
    ConvectiveArgument = (1.0_kind_noahmp - 10.0_kind_noahmp*Parameter)**(1.0_kind_noahmp/3.0_kind_noahmp)
    ConvectiveCorrection = 1.5_kind_noahmp*log((ConvectiveArgument**2+ConvectiveArgument+1.0_kind_noahmp)/3.0_kind_noahmp) - &
         SqrtThree*atan((2.0_kind_noahmp*ConvectiveArgument+1.0_kind_noahmp)/SqrtThree) + &
         4.0_kind_noahmp*ArcTangentOne/SqrtThree
    Correction = (BusingerCorrection + Parameter**2*ConvectiveCorrection) / &
         (1.0_kind_noahmp + Parameter**2)
  end function MynnMomentumUnstable


  pure real(kind=kind_noahmp) function MynnHeatUnstable(Parameter) result(Correction)
    real(kind=kind_noahmp), intent(in) :: Parameter
    real(kind=kind_noahmp) :: BusingerArgument
    real(kind=kind_noahmp) :: ConvectiveArgument
    real(kind=kind_noahmp) :: BusingerCorrection
    real(kind=kind_noahmp) :: ConvectiveCorrection
    real(kind=kind_noahmp), parameter :: SqrtThree = 1.7320508075688773_kind_noahmp
    real(kind=kind_noahmp), parameter :: ArcTangentOne = 0.7853981633974483_kind_noahmp

    BusingerArgument = sqrt(1.0_kind_noahmp - 16.0_kind_noahmp*Parameter)
    BusingerCorrection = 2.0_kind_noahmp*log(0.5_kind_noahmp*(1.0_kind_noahmp+BusingerArgument))
    ConvectiveArgument = (1.0_kind_noahmp - 34.0_kind_noahmp*Parameter)**(1.0_kind_noahmp/3.0_kind_noahmp)
    ConvectiveCorrection = 1.5_kind_noahmp*log((ConvectiveArgument**2+ConvectiveArgument+1.0_kind_noahmp)/3.0_kind_noahmp) - &
         SqrtThree*atan((2.0_kind_noahmp*ConvectiveArgument+1.0_kind_noahmp)/SqrtThree) + &
         4.0_kind_noahmp*ArcTangentOne/SqrtThree
    Correction = (BusingerCorrection + Parameter**2*ConvectiveCorrection) / &
         (1.0_kind_noahmp + Parameter**2)
  end function MynnHeatUnstable


  pure real(kind=kind_noahmp) function GfsMomentumStable(Parameter) result(Correction)
    real(kind=kind_noahmp), intent(in) :: Parameter
    real(kind=kind_noahmp) :: RootTerm
    RootTerm = sqrt(1.0_kind_noahmp + 20.0_kind_noahmp*Parameter)
    Correction = -RootTerm + log(RootTerm + 1.0_kind_noahmp)
  end function GfsMomentumStable


  pure real(kind=kind_noahmp) function GfsHeatStable(Parameter) result(Correction)
    real(kind=kind_noahmp), intent(in) :: Parameter
    Correction = GfsMomentumStable(Parameter)
  end function GfsHeatStable


  pure real(kind=kind_noahmp) function GfsMomentumUnstable(Parameter) result(Correction)
    real(kind=kind_noahmp), intent(in) :: Parameter
    real(kind=kind_noahmp) :: PositiveParameter
    real(kind=kind_noahmp) :: RootTerm

    if ( Parameter >= -0.5_kind_noahmp ) then
       Correction = (-3.975_kind_noahmp + 12.32_kind_noahmp*Parameter) * Parameter / &
            (1.0_kind_noahmp + (-7.755_kind_noahmp + 6.041_kind_noahmp*Parameter)*Parameter)
    else
       PositiveParameter = -Parameter
       RootTerm = 1.0_kind_noahmp / sqrt(PositiveParameter)
       Correction = log(PositiveParameter) + 2.0_kind_noahmp*sqrt(RootTerm) - 0.8776_kind_noahmp
    endif
  end function GfsMomentumUnstable


  pure real(kind=kind_noahmp) function GfsHeatUnstable(Parameter) result(Correction)
    real(kind=kind_noahmp), intent(in) :: Parameter
    real(kind=kind_noahmp) :: PositiveParameter
    real(kind=kind_noahmp) :: RootTerm

    if ( Parameter >= -0.5_kind_noahmp ) then
       Correction = (-7.941_kind_noahmp + 24.75_kind_noahmp*Parameter) * Parameter / &
            (1.0_kind_noahmp + (-8.705_kind_noahmp + 7.899_kind_noahmp*Parameter)*Parameter)
    else
       PositiveParameter = -Parameter
       RootTerm = 1.0_kind_noahmp / sqrt(PositiveParameter)
       Correction = log(PositiveParameter) + 0.5_kind_noahmp*RootTerm + 1.386_kind_noahmp
    endif
  end function GfsHeatUnstable

end module StabilityFunctionsMod
