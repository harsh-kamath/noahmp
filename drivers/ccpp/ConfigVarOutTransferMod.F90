! Harsh Kamath (NOAA GSL/CIRES), September 2026


module ConfigVarOutTransferMod

  use Machine
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpVarType, only: noahmp_type

  implicit none

contains

!=== Transfer model states to output=====

  subroutine ConfigVarOutTransfer(noahmp, NoahmpIO, N)

    implicit none

    type(noahmp_type),        intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO
	integer,                  intent(in)    :: N

! ----------------------------------------------------------------------
    associate(                                                         &
              NumSnowLayerMax => noahmp%config%domain%NumSnowLayerMax ,&
              NumSoilLayer    => noahmp%config%domain%NumSoilLayer     &
             )
! ----------------------------------------------------------------------

    ! config domain variables
    NoahmpIO%ISNOWXY(N)  = noahmp%config%domain%NumSnowLayerNeg
    NoahmpIO%ZSNSOXY(-NumSnowLayerMax+1:NumSoilLayer,N) = &
                            noahmp%config%domain%DepthSnowSoilLayer(-NumSnowLayerMax+1:NumSoilLayer)
    NoahmpIO%FORCZLSM(N) = noahmp%config%domain%RefHeightAboveSfc

    end associate

  end subroutine ConfigVarOutTransfer

end module ConfigVarOutTransferMod
