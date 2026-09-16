! Harsh Kamath (NOAA GSL/CIRES), September 2026


module NoahmpDriverMainMod

  use NoahmpVarType,                  only : noahmp_type
  use NoahmpIOVarType,         only : NoahmpIO_type
  use ConfigVarInitMod,               only : ConfigVarInitDefault
  use ForcingVarInitMod,              only : ForcingVarInitDefault
  use EnergyVarInitMod,               only : EnergyVarInitDefault
  use WaterVarInitMod,                only : WaterVarInitDefault
  use BiochemVarInitMod,              only : BiochemVarInitDefault
  use ConfigVarInTransferMod,    only : ConfigVarInTransfer
  use ForcingVarInTransferMod,   only : ForcingVarInTransfer
  use EnergyVarInTransferMod,    only : EnergyVarInTransfer
  use WaterVarInTransferMod,     only : WaterVarInTransfer
  use BiochemVarInTransferMod,   only : BiochemVarInTransfer
  use ConfigVarOutTransferMod,   only : ConfigVarOutTransfer
  use ForcingVarOutTransferMod,  only : ForcingVarOutTransfer
  use EnergyVarOutTransferMod,   only : EnergyVarOutTransfer
  use WaterVarOutTransferMod,    only : WaterVarOutTransfer
  use BiochemVarOutTransferMod,  only : BiochemVarOutTransfer
  use NoahmpColumnPhysicsMod,         only : NoahmpRunColumnPhysics

  implicit none
  private

  public :: NoahmpDriverMain

contains

  subroutine NoahmpDriverMain(NoahmpIO, apply_urban_irrigation, &
                            errmsg, errflg)

    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    logical, intent(in)                  :: apply_urban_irrigation(:)
    character(len=*), intent(out)        :: errmsg
    integer, intent(out)                 :: errflg

    type(noahmp_type) :: column
    integer           :: n

    errmsg = ''
    errflg = 0

    if (NoahmpIO%ncol < 0 .or. NoahmpIO%ncol > NoahmpIO%capacity) then
       errflg = 1
       errmsg = 'NoahmpDriverMain: ncol is outside [0, capacity]'
       return
    end if

    if (size(apply_urban_irrigation) < NoahmpIO%ncol) then
       errflg = 1
       errmsg = 'NoahmpDriverMain: urban-irrigation mask is smaller than ncol'
       return
    end if

    do n = 1, NoahmpIO%ncol
       call ConfigVarInitDefault(column)
       call ConfigVarInTransfer(column, NoahmpIO, n)

       call ForcingVarInitDefault(column)
       call ForcingVarInTransfer(column, NoahmpIO, n)

       call EnergyVarInitDefault(column)
       call EnergyVarInTransfer(column, NoahmpIO, n)

       call WaterVarInitDefault(column)
       call WaterVarInTransfer(column, NoahmpIO, n)

       call BiochemVarInitDefault(column)
       call BiochemVarInTransfer(column, NoahmpIO, n)

       if (apply_urban_irrigation(n)) then
          if (NoahmpIO%NSOIL < 2) then
             errflg = 1
             errmsg = 'NoahmpDriverMain: urban irrigation requires two soil layers'
             return
          end if

          column%water%state%SoilMoisture(1) = max( &
               column%water%state%SoilMoisture(1), &
               column%water%param%SoilMoistureFieldCap(1))
          column%water%state%SoilMoisture(2) = max( &
               column%water%state%SoilMoisture(2), &
               column%water%param%SoilMoistureFieldCap(2))
       end if

       call NoahmpRunColumnPhysics(column)

       call ConfigVarOutTransfer(column, NoahmpIO, n)
       call ForcingVarOutTransfer(column, NoahmpIO, n)
       call EnergyVarOutTransfer(column, NoahmpIO, n)
       call WaterVarOutTransfer(column, NoahmpIO, n)
       call BiochemVarOutTransfer(column, NoahmpIO, n)
    end do

  end subroutine NoahmpDriverMain

end module NoahmpDriverMainMod
