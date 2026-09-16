! Harsh Kamath (NOAA GSL/CIRES), September 2026


module NoahmpColumnPhysicsMod

  use Machine,              only : kind_noahmp
  use NoahmpVarType,        only : noahmp_type
  use NoahmpMainMod,        only : NoahmpMain
  use NoahmpMainGlacierMod, only : NoahmpMainGlacier

  implicit none
  private

  public :: NoahmpRunColumnPhysics

contains

  subroutine NoahmpRunColumnPhysics(noahmp)

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

    if (noahmp%config%domain%VegType == &
        noahmp%config%domain%IndexIcePoint) then

       noahmp%config%domain%IndicatorIceSfc = -1
       noahmp%forcing%TemperatureSoilBottom = &
            min(noahmp%forcing%TemperatureSoilBottom, 263.15_kind_noahmp)

       call NoahmpMainGlacier(noahmp)

    else

       noahmp%config%domain%IndicatorIceSfc = 0
       call NoahmpMain(noahmp)

    end if

  end subroutine NoahmpRunColumnPhysics

end module NoahmpColumnPhysicsMod
