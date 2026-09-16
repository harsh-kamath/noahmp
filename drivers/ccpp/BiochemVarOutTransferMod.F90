! Harsh Kamath (NOAA GSL/CIRES), September 2026


module BiochemVarOutTransferMod

  use Machine
  use NoahmpIOVarType, only: NoahmpIO_type
  use NoahmpVarType, only: noahmp_type

  implicit none

contains

!=== Transfer model states to output =====

  subroutine BiochemVarOutTransfer(noahmp, NoahmpIO, N)

    implicit none

    type(noahmp_type),        intent(inout) :: noahmp
    type(NoahmpIO_type),      intent(inout) :: NoahmpIO
	integer,                  intent(in)    :: N

    ! biochem state variables
    NoahmpIO%LFMASSXY(N) = noahmp%biochem%state%LeafMass
    NoahmpIO%RTMASSXY(N) = noahmp%biochem%state%RootMass
    NoahmpIO%STMASSXY(N) = noahmp%biochem%state%StemMass
    NoahmpIO%WOODXY  (N) = noahmp%biochem%state%WoodMass
    NoahmpIO%STBLCPXY(N) = noahmp%biochem%state%CarbonMassDeepSoil
    NoahmpIO%FASTCPXY(N) = noahmp%biochem%state%CarbonMassShallowSoil
    NoahmpIO%GDDXY   (N) = noahmp%biochem%state%GrowDegreeDay
    NoahmpIO%PGSXY   (N) = noahmp%biochem%state%PlantGrowStage
    NoahmpIO%GRAINXY (N) = noahmp%biochem%state%GrainMass

    ! biochem flux variables
    NoahmpIO%NEEXY   (N) = noahmp%biochem%flux%NetEcoExchange
    NoahmpIO%GPPXY   (N) = noahmp%biochem%flux%GrossPriProduction
    NoahmpIO%NPPXY   (N) = noahmp%biochem%flux%NetPriProductionTot
    NoahmpIO%PSNXY   (N) = noahmp%biochem%flux%PhotosynTotal

  end subroutine BiochemVarOutTransfer

end module BiochemVarOutTransferMod
