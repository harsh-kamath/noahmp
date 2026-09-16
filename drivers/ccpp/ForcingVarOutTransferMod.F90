! Harsh Kamath (NOAA GSL/CIRES), September 2026


module ForcingVarOutTransferMod

  use Machine
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpVarType, only: noahmp_type

  implicit none

contains

!=== Transfer model states to output =====

  subroutine ForcingVarOutTransfer(noahmp, NoahmpIO, N)

    implicit none

    type(noahmp_type),        intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO
	integer,                  intent(in)    :: N

    NoahmpIO%FORCTLSM  (N) = noahmp%forcing%TemperatureAirRefHeight
    NoahmpIO%FORCQLSM  (N) = noahmp%forcing%SpecHumidityRefHeight
    NoahmpIO%FORCPLSM  (N) = noahmp%forcing%PressureAirRefHeight
    NoahmpIO%FORCWLSM  (N) = sqrt(noahmp%forcing%WindEastwardRefHeight**2 + &
                             noahmp%forcing%WindNorthwardRefHeight**2)
    NoahmpIO%RadSwDirFrac(N) = noahmp%forcing%RadSwDirFrac
    NoahmpIO%RadSwVisFrac(N) = noahmp%forcing%RadSwVisFrac

  end subroutine ForcingVarOutTransfer

end module ForcingVarOutTransferMod
