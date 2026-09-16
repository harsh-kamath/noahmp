! Harsh Kamath (NOAA GSL/CIRES), September 2026


module ForcingVarInTransferMod

  use Machine, only: kind_noahmp
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpVarType, only: noahmp_type
  
  implicit none

contains

!=== initialize with input data or table values

  subroutine ForcingVarInTransfer(noahmp, NoahmpIO, N)

    implicit none

    type(noahmp_type),        intent(inout) :: noahmp
    type(NoahmpIO_type), intent(inout) :: NoahmpIO
    integer,                  intent(in)    :: N
    
    ! local variables
	integer :: K
    real(kind=kind_noahmp)              :: PrecipOtherRefHeight  ! other precipitation, e.g. fog [mm/s] at reference height
    real(kind=kind_noahmp)              :: PrecipTotalRefHeight  ! total precipitation [mm/s] at reference height
	
	K = NoahmpIO%KTS

    noahmp%forcing%TemperatureAirRefHeight = NoahmpIO%T_PHY(K,N)
    noahmp%forcing%WindEastwardRefHeight   = NoahmpIO%U_PHY(K,N)
    noahmp%forcing%WindNorthwardRefHeight  = NoahmpIO%V_PHY(K,N)
    ! The modular interface has already supplied CCPP specific humidity.
    noahmp%forcing%SpecHumidityRefHeight   = NoahmpIO%QV_CURR(K,N)
    noahmp%forcing%PressureAirRefHeight    = NoahmpIO%PRSL1(N)
    noahmp%forcing%PressureAirSurface      = NoahmpIO%PS(N)
    noahmp%forcing%BoundaryLayerHeight     = NoahmpIO%PBLH(N)
    noahmp%forcing%RadLwDownRefHeight      = NoahmpIO%GLW      (N)
    noahmp%forcing%RadSwDownRefHeight      = NoahmpIO%SWDOWN   (N)
    noahmp%forcing%TemperatureSoilBottom   = NoahmpIO%TMN      (N)

    ! treat different precipitation types
    PrecipTotalRefHeight                   = NoahmpIO%RAINBL   (N)
    noahmp%forcing%PrecipConvRefHeight     = NoahmpIO%MP_RAINC (N)
    noahmp%forcing%PrecipNonConvRefHeight  = NoahmpIO%MP_RAINNC(N)
    noahmp%forcing%PrecipShConvRefHeight   = NoahmpIO%MP_SHCV  (N)
    noahmp%forcing%PrecipSnowRefHeight     = NoahmpIO%MP_SNOW  (N)
    noahmp%forcing%PrecipGraupelRefHeight  = NoahmpIO%MP_GRAUP (N)
    noahmp%forcing%PrecipHailRefHeight     = NoahmpIO%MP_HAIL  (N)
    ! treat other precipitation (e.g. fog) contained in total precipitation
    PrecipOtherRefHeight                   = PrecipTotalRefHeight - noahmp%forcing%PrecipConvRefHeight - &
                                             noahmp%forcing%PrecipNonConvRefHeight - noahmp%forcing%PrecipShConvRefHeight
    PrecipOtherRefHeight                   = max(0.0, PrecipOtherRefHeight)
    noahmp%forcing%PrecipNonConvRefHeight  = noahmp%forcing%PrecipNonConvRefHeight + PrecipOtherRefHeight
    noahmp%forcing%PrecipSnowRefHeight     = noahmp%forcing%PrecipSnowRefHeight + PrecipOtherRefHeight * NoahmpIO%SR(N)

    ! downward solar radiation direct/diffuse and visible/NIR partition
    noahmp%forcing%RadSwDirFrac            = NoahmpIO%RadSwDirFrac(N)
    noahmp%forcing%RadSwVisFrac            = NoahmpIO%RadSwVisFrac(N)

    ! SNICAR aerosol deposition flux forcing
    if ( noahmp%config%nmlist%OptSnowAlbedo == 3 ) then
       if ( noahmp%config%nmlist%FlagSnicarAerosolReadTable ) then
          noahmp%forcing%DepBChydropho = NoahmpIO%DepBChydropho_TABLE
          noahmp%forcing%DepBChydrophi = NoahmpIO%DepBChydrophi_TABLE
          noahmp%forcing%DepOChydropho = NoahmpIO%DepOChydropho_TABLE
          noahmp%forcing%DepOChydrophi = NoahmpIO%DepOChydrophi_TABLE
          noahmp%forcing%DepDust1 = NoahmpIO%DepDust1_TABLE
          noahmp%forcing%DepDust2 = NoahmpIO%DepDust2_TABLE
          noahmp%forcing%DepDust3 = NoahmpIO%DepDust3_TABLE
          noahmp%forcing%DepDust4 = NoahmpIO%DepDust4_TABLE
          noahmp%forcing%DepDust5 = NoahmpIO%DepDust5_TABLE
       else
          noahmp%forcing%DepBChydropho = NoahmpIO%DepBChydrophoXY(N)
          noahmp%forcing%DepBChydrophi = NoahmpIO%DepBChydrophiXY(N)
          noahmp%forcing%DepOChydropho = NoahmpIO%DepOChydrophoXY(N)
          noahmp%forcing%DepOChydrophi = NoahmpIO%DepOChydrophiXY(N)
          noahmp%forcing%DepDust1 = NoahmpIO%DepDust1XY(N)
          noahmp%forcing%DepDust2 = NoahmpIO%DepDust2XY(N)
          noahmp%forcing%DepDust3 = NoahmpIO%DepDust3XY(N)
          noahmp%forcing%DepDust4 = NoahmpIO%DepDust4XY(N)
          noahmp%forcing%DepDust5 = NoahmpIO%DepDust5XY(N)
       endif
    endif

 
  end subroutine ForcingVarInTransfer

end module ForcingVarInTransferMod
