! Harsh Kamath (NOAA GSL/CIRES), September 2026


module NoahmpIOVarType

  use Machine

  implicit none
  save
  private

  type, public :: NoahmpIO_type

!------------------------------------------------------------------------
! General packed-column Noah-MP variables
!------------------------------------------------------------------------

    ! Number of active columns packed into this CCPP call.
    integer                                                ::  ncol
    ! Allocated capacity for packed columns.
	integer                                                ::  capacity
    ! Index of the lowest atmospheric model level.
	integer 											   ::  kts
    ! Bounds of the allocated atmospheric level dimension.
	integer                                                ::  kms,kme
    ! Original host horizontal I index for each packed column.
	integer, allocatable, dimension(:) 					   ::  source_i
    ! Original host horizontal J index for each packed column.
    integer, allocatable, dimension(:) 					   ::  source_j

   

    ! Configuration and packed-column geometry
    integer                                                ::  NSNOW = 3            ! number of snow layers fixed to 3
    integer                                                ::  NSOIL               ! number of soil layers

    integer                                                ::  IOPT_DVEG           ! dynamic vegetation
    integer                                                ::  IOPT_SNF            ! rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
    integer                                                ::  IOPT_BTR            ! soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
    integer                                                ::  IOPT_RSF            ! surface resistance option (1->Zeng; 2->simple)
    integer                                                ::  IOPT_SFC            ! surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
    integer                                                ::  PSI_OPT = 0         ! MYNN stability functions (0->MYNN; 1->GFS)
    integer                                                ::  IZ0TLND = 0         ! MYNN scalar roughness option (0 through 3)
    integer                                                ::  ITIMESTEP = 1       ! current host-model timestep index
    logical                                                ::  THSFC_LOC = .true.  ! use local pressure reference for potential temperature
    integer                                                ::  IOPT_CRS            ! canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
    integer                                                ::  IOPT_ALB            ! snow surface albedo (1->BATS; 2->CLASS; 3->SNICAR)
    integer                                                ::  IOPT_RAD            ! radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
    integer                                                ::  IOPT_STC            ! snow/soil temperature time scheme
    integer                                                ::  IOPT_TKSNO          ! snow thermal conductivity: 1 -> Stieglitz(yen,1965) scheme (default), 2 -> Anderson, 1976 scheme, 3 -> constant, 4 -> Verseghy (1991) scheme, 5 -> Douvill(Yen, 1981) scheme
    integer                                                ::  IOPT_TBOT           ! lower boundary of soil temperature (1->zero-flux; 2->Noah)
    integer                                                ::  IOPT_FRZ            ! supercooled liquid water (1-> NY06; 2->Koren99)
    integer                                                ::  IOPT_INF            ! frozen soil permeability (1-> NY06; 2->Koren99)
    integer                                                ::  IOPT_INFDV          ! infiltration options for dynamic VIC (1->Philip; 2-> Green-Ampt;3->Smith-Parlange)
    integer                                                ::  IOPT_TDRN           ! drainage option (0->off; 1->simple scheme; 2->Hooghoudt's scheme)
    integer                                                ::  IOPT_IRR            ! irrigation scheme (0->none; >1 irrigation scheme ON)
    integer                                                ::  IOPT_IRRM           ! irrigation method (0->dynamic; 1-> sprinkler; 2-> micro; 3-> flood)
    integer                                                ::  IOPT_CROP           ! crop model option (0->none; 1->Liu et al.)
    integer                                                ::  IOPT_SOIL           ! soil configuration option
    integer                                                ::  IOPT_PEDO           ! soil pedotransfer function option
    integer                                                ::  IOPT_RUNSRF         ! surface runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
    integer                                                ::  IOPT_RUNSUB         ! subsurface runoff option
    integer                                                ::  IOPT_GLA            ! glacier option (1->phase change; 2->simple)
    integer                                                ::  IOPT_COMPACT        ! snowpack compaction (1->Anderson1976; 2->Abolafia-Rosenzweig2024)
    integer                                                ::  IOPT_WETLAND        ! wetland model option (0->off; 1->Zhang2022 fixed parameter; 2->Zhang2022 read in 2D parameter)
    integer                                                ::  IOPT_SCF            ! snow cover fraction (1->NiuYang07; 2->Abolafia-Rosenzweig2025)

    integer                                                ::  SNICAR_SNOWSHAPE_OPT          !option for snow grain shape in SNICAR (He et al. 2017 JC)
    integer                                                ::  SNICAR_RTSOLVER_OPT           !option for two different SNICAR radiative transfer solver
    integer                                                ::  SNICAR_BANDNUMBER_OPT         !number of wavelength bands used in SNICAR snow albedo calculation
    integer                                                ::  SNICAR_SOLARSPEC_OPT          !type of downward solar radiation spectrum for SNICAR snow albedo calculation
    integer                                                ::  SNICAR_SNOWOPTICS_OPT         !snow optics type using different refractive index databases in SNICAR
    integer                                                ::  SNICAR_DUSTOPTICS_OPT         !dust optics type for SNICAR snow albedo calculation
    logical                                                ::  SNICAR_SNOWBC_INTMIX          !option to activate BC-snow internal mixing in SNICAR (He et al. 2017 JC)
    logical                                                ::  SNICAR_SNOWDUST_INTMIX        !option to activate dust-snow internal mixing in SNICAR (He et al. 2017 JC)
    logical                                                ::  SNICAR_USE_AEROSOL            !option to turn on/off aerosol deposition flux effect in snow in SNICAR
    logical                                                ::  SNICAR_USE_OC                 !option to activate OC in snow in SNICAR
    logical                                                ::  SNICAR_AEROSOL_READTABLE      !option to read aerosol deposition fluxes from table (on) or NetCDF forcing file (off)

    integer                                                ::  NUMRAD = 2          ! number of shortwave band
    logical                                                ::  calculate_soil      ! logical index for if do soil calculation
    integer                                                ::  soil_update_steps   ! number of model time steps to update soil process
    integer, allocatable, dimension(:)                   ::  ISNOWXY             ! actual no. of snow layers

    real(kind=kind_noahmp)                                 ::  DTBL                ! timestep [s]
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  DX                  ! horizontal grid spacing [m]
    ! Horizontal grid spacing in the Y direction [m].
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  DY

    CHARACTER(LEN=256)                                     ::  LLANDUSE            ! (=USGS, using USGS landuse classification)

    integer,                allocatable, dimension(:)      ::  IVGTYP              ! vegetation type
    integer, allocatable, dimension(:)                     :: CROPCAT              ! crop category
    integer, allocatable, dimension(:)                     ::  ICE                 ! ice-surface flag (land ice = -1)
    integer, allocatable, dimension(:)                     ::  SOILCOL             ! soil color category

    real(kind=kind_noahmp)                                 ::  JULIAN              ! Julian day
    ! Number of days in the current year.
    integer                                                ::  YEARLEN

    real(kind=kind_noahmp), allocatable, dimension(:)    ::  XLAT                ! latitude [rad]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  DZ8W                ! thickness of atmo layers [m]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  COSZEN              ! cosine zenith angle

    integer                                                :: ISWATER_TABLE             ! water flag
    integer                                                :: ISBARREN_TABLE            ! barren ground flag
    integer                                                :: ISICE_TABLE               ! ice flag
    integer                                                :: ISCROP_TABLE              ! cropland flag
    integer                                                :: EBLFOREST_TABLE           ! evergreen broadleaf forest flag

    ! Runoff slope category for each packed column.
    integer, allocatable, dimension(:)                     ::  SLOPETYP
    real(kind=kind_noahmp)                                 ::  ZBOT_TABLE                ! Depth [m] of lower boundary soil temperature

    integer                                                ::  idx_T_max         = 11        ! maxiumum temperature index used in aging lookup table [idx]
    integer                                                ::  idx_Tgrd_max      = 31        ! maxiumum temperature gradient index used in aging lookup table [idx]
    integer                                                ::  idx_rhos_max      = 8         ! maxiumum snow density index used in aging lookup table [idx]
    integer                                                ::  snicar_numrad_snw             ! wavelength bands used in SNICAR snow albedo calculation
    integer                                                ::  idx_Mie_snw_mx    = 1471      ! number of effective radius indices used in Mie lookup table [idx]

    integer,                allocatable, dimension(:)    ::  ISLTYP              ! soil type
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  soilcl1             ! Soil texture class with depth
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  soilcl2             ! Soil texture class with depth
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  soilcl3             ! Soil texture class with depth
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  soilcl4             ! Soil texture class with depth

    real(kind=kind_noahmp), allocatable, dimension(:)      ::  ZSOIL               ! depth to soil interfaces [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  ZSNSOXY             ! snow layer depth [m]

    integer                                                ::  ISURBAN_TABLE             ! urban flag
    integer                                                ::  URBTYPE_beg         ! urban type start number - 1
    ! Host urban-physics option.
    integer                                                ::  sf_urban_physics
    integer                                                ::  NATURAL_TABLE             ! natural vegetation type

    real(kind=kind_noahmp), allocatable, dimension(:)    ::  GVFMAX              ! annual maximum in vegetation fraction
    integer                                                ::  DEFAULT_CROP_TABLE        ! Default crop index
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  XICE                ! fraction of grid that is seaice
	

    ! forcings in transfers   
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  T_PHY               ! 3D atmospheric temperature valid at mid-levels [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  QV_CURR             ! specific humidity at lowest model layer [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  U_PHY               ! 3D U wind component [m/s]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  V_PHY               ! 3D V wind component [m/s]
	
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SWDOWN              ! solar down at surface [W m-2]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  GLW                 ! longwave down at surface [W m-2]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  PS                  ! surface air pressure [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  PRSL1               ! air pressure at lowest model layer [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  PBLH                ! planetary boundary layer height [m]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  RAINBL              ! total precipitation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SR                  ! frozen precipitation fraction/type [-]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  MP_RAINC            ! convective precipitation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  MP_RAINNC           ! nonconvective precipitation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  MP_SHCV             ! shallow-convective precipitation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  MP_SNOW             ! snow precipitation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  MP_GRAUP            ! graupel precipitation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  MP_HAIL             ! hail precipitation rate [mm/s]
	
    ! deep soil temperature [K]
	real(kind=kind_noahmp), allocatable, dimension(:)    ::  TMN
	real(kind=kind_noahmp), allocatable, dimension(:)    ::  RadSwDirFrac          ! direct shortwave fraction [-]
	real(kind=kind_noahmp), allocatable, dimension(:)    ::  RadSwVisFrac          ! visible shortwave fraction [-]
    

	! forcing out transfers
	real(kind=kind_noahmp), allocatable, dimension(:)    ::  FORCTLSM            ! surface temperature as LSM forcing [K]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FORCQLSM            ! surface specific humidity as LSM forcing [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FORCPLSM            ! surface pressure as LSM forcing [Pa]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FORCZLSM            ! reference height as LSM input [m]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FORCWLSM            ! surface wind speed as LSM forcing [m/s]

	 
	! energy in transfer variables
    ! leaf area index
	real(kind=kind_noahmp), allocatable, dimension(:)   :: LAI
    ! stem area index
    real(kind=kind_noahmp), allocatable, dimension(:)   :: XSAIXY
    ! bulk surface specific humidity
    real(kind=kind_noahmp), allocatable, dimension(:)   :: QSFC
    ! bulk ground surface temperature
    real(kind=kind_noahmp), allocatable, dimension(:)   :: TGXY
    ! vegetation leaf temperature
    real(kind=kind_noahmp), allocatable, dimension(:)   :: TVXY
    ! snow age factor
    real(kind=kind_noahmp), allocatable, dimension(:)   :: TAUSSXY
    ! snow albedo at last time step (-)
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ALBOLDXY
    ! canopy air vapor pressure (pa)
    real(kind=kind_noahmp), allocatable, dimension(:)   :: EAHXY
    ! canopy air temperature (k)
    real(kind=kind_noahmp), allocatable, dimension(:)   :: TAHXY
    ! bulk sensible heat exchange coefficient
    real(kind=kind_noahmp), allocatable, dimension(:)   :: CHXY
    ! bulk momentum drag coefficient
    real(kind=kind_noahmp), allocatable, dimension(:)   :: CMXY
    ! snow temperature [K]
	real(kind=kind_noahmp), allocatable, dimension(:,:) :: TSNOXY
    ! soil temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TSLB
    ! soil albedo (direct)
	real(kind=kind_noahmp), allocatable, dimension(:,:) :: ALBSOILDIRXY
    ! soil albedo (diffuse)
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ALBSOILDIFXY
    ! accumulated ground heat flux [W/m2 * dt_soil/dt_main]
	real(kind=kind_noahmp), allocatable, dimension(:)   :: ACC_SSOILXY
    ! vegetation fraction []
    real(kind=kind_noahmp), allocatable, dimension(:)   :: VEGFRA
	
   ! Scalar table parameters
    ! co2 partial pressure
    real(kind=kind_noahmp) :: CO2_TABLE
    ! o2 partial pressure
    real(kind=kind_noahmp) :: O2_TABLE
    ! Soil heat capacity [J m-3 K-1]
    real(kind=kind_noahmp) :: CSOIL_TABLE
    ! tau0 from Yang97 eqn. 10a
    real(kind=kind_noahmp) :: TAU0_TABLE
    ! growth from vapor diffusion Yang97 eqn. 10b
    real(kind=kind_noahmp) :: GRAIN_GROWTH_TABLE
    ! dirt and soot term Yang97 eqn. 10d
    real(kind=kind_noahmp) :: DIRT_SOOT_TABLE
    ! extra growth near freezing Yang97 eqn. 10c
    real(kind=kind_noahmp) :: EXTRA_GROWTH_TABLE
    ! zenith angle snow albedo adjustment; b in Yang97 eqn. 15
    real(kind=kind_noahmp) :: BATS_COSZ_TABLE
    ! new snow visible albedo
    real(kind=kind_noahmp) :: BATS_VIS_NEW_TABLE
    ! new snow NIR albedo
    real(kind=kind_noahmp) :: BATS_NIR_NEW_TABLE
    ! age factor for diffuse visible snow albedo Yang97 eqn. 17
    real(kind=kind_noahmp) :: BATS_VIS_AGE_TABLE
    ! age factor for diffuse NIR snow albedo Yang97 eqn. 18
    real(kind=kind_noahmp) :: BATS_NIR_AGE_TABLE
    ! cosz factor for direct visible snow albedo Yang97 eqn. 15
    real(kind=kind_noahmp) :: BATS_VIS_DIR_TABLE
    ! cosz factor for direct NIR snow albedo Yang97 eqn. 16
    real(kind=kind_noahmp) :: BATS_NIR_DIR_TABLE
    ! reference snow albedo in CLASS scheme
    real(kind=kind_noahmp) :: CLASS_ALB_REF_TABLE
    ! snow aging e-folding time (s) in CLASS albedo scheme
    real(kind=kind_noahmp) :: CLASS_SNO_AGE_TABLE
    ! fresh snow albedo in CLASS scheme
    real(kind=kind_noahmp) :: CLASS_ALB_NEW_TABLE
    ! two-stream parameter betad for snow
    real(kind=kind_noahmp) :: BETADS_TABLE
    ! two-stream parameter betad for snow
    real(kind=kind_noahmp) :: BETAIS_TABLE
    ! Parameter used in the calculation of the roughness length for heat
    real(kind=kind_noahmp) :: CZIL_TABLE
    ! snow emissivity
    real(kind=kind_noahmp) :: SNOW_EMIS_TABLE
    ! snow surface roughness length (m) (0.002)
    real(kind=kind_noahmp) :: Z0SNO_TABLE
    ! Bare-soil roughness length (m) (i.e., under the canopy)
    real(kind=kind_noahmp) :: Z0SOIL_TABLE
    ! Lake surface roughness length (m)
    real(kind=kind_noahmp) :: Z0LAKE_TABLE
    ! ice surface emissivity
    real(kind=kind_noahmp) :: EICE_TABLE
    ! exponent in the shape parameter for soil resistance option 1
    real(kind=kind_noahmp) :: RSURF_EXP_TABLE
    ! surface resistance for snow(s/m)
    real(kind=kind_noahmp) :: RSURF_SNOW_TABLE

   ! 1D vegetation & surface table parameters
    ! emissivity soil surface
    real(kind=kind_noahmp), allocatable, dimension(:) :: EG_TABLE
    ! albedo land ice: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:) :: ALBICE_TABLE
    ! tree crown radius (m)
    real(kind=kind_noahmp), allocatable, dimension(:) :: RC_TABLE
    ! top of canopy (m)
    real(kind=kind_noahmp), allocatable, dimension(:) :: HVT_TABLE
    ! bottom of canopy (m)
    real(kind=kind_noahmp), allocatable, dimension(:) :: HVB_TABLE
    ! momentum roughness length (m)
    real(kind=kind_noahmp), allocatable, dimension(:) :: Z0MVT_TABLE
    ! empirical canopy wind parameter
    real(kind=kind_noahmp), allocatable, dimension(:) :: CWPVT_TABLE
    ! tree density (no. of trunks per m2)
    real(kind=kind_noahmp), allocatable, dimension(:) :: DEN_TABLE
    ! leaf/stem orientation index
    real(kind=kind_noahmp), allocatable, dimension(:) :: XL_TABLE
    ! minimum leaf conductance (umol/m2/s)
    real(kind=kind_noahmp), allocatable, dimension(:) :: BP_TABLE
    ! co2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp), allocatable, dimension(:) :: KC25_TABLE
    ! o2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp), allocatable, dimension(:) :: KO25_TABLE
    ! q10 for kc25
    real(kind=kind_noahmp), allocatable, dimension(:) :: AKC_TABLE
    ! q10 for ko25
    real(kind=kind_noahmp), allocatable, dimension(:) :: AKO_TABLE
    ! Parameter used in radiation stress function
    real(kind=kind_noahmp), allocatable, dimension(:) :: RGL_TABLE
    ! Minimum stomatal resistance [s m-1]
    real(kind=kind_noahmp), allocatable, dimension(:) :: RS_TABLE
    ! Maximal stomatal resistance [s m-1]
    real(kind=kind_noahmp), allocatable, dimension(:) :: RSMAX_TABLE
    ! Optimum transpiration air temperature [K]
    real(kind=kind_noahmp), allocatable, dimension(:) :: TOPT_TABLE
    ! Parameter used in vapor pressure deficit function
    real(kind=kind_noahmp), allocatable, dimension(:) :: HS_TABLE
    ! characteristic leaf dimension (m)
    real(kind=kind_noahmp), allocatable, dimension(:) :: DLEAF_TABLE
    ! canopy biomass heat capacity parameter (m)
    real(kind=kind_noahmp), allocatable, dimension(:) :: CBIOM_TABLE
    ! albedo frozen lakes: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:) :: ALBLAK_TABLE
    ! two-stream parameter omega for snow
    real(kind=kind_noahmp), allocatable, dimension(:) :: OMEGAS_TABLE
    ! soil quartz content
    real(kind=kind_noahmp), allocatable, dimension(:) :: QUARTZ_TABLE
    ! Soil quartz content
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: quartz_3D

    ! Crop-specific table parameters
    ! minimum leaf conductance (umol/m2/s)
    real(kind=kind_noahmp), allocatable, dimension(:) :: BPI_TABLE
    ! co2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp), allocatable, dimension(:) :: KC25I_TABLE
    ! o2 michaelis-menten constant at 25c (pa)
    real(kind=kind_noahmp), allocatable, dimension(:) :: KO25I_TABLE
    ! q10 for kc25
    real(kind=kind_noahmp), allocatable, dimension(:) :: AKCI_TABLE
    ! q10 for ko25
    real(kind=kind_noahmp), allocatable, dimension(:) :: AKOI_TABLE

    ! 2D vegetation optical table parameters
    ! monthly leaf area index, one-sided
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: LAIM_TABLE
    ! monthly stem area index, one-sided
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SAIM_TABLE
    ! leaf reflectance: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RHOL_TABLE
    ! stem reflectance: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: RHOS_TABLE
    ! leaf transmittance: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TAUL_TABLE
    ! stem transmittance: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: TAUS_TABLE
    ! saturated soil albedos: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ALBSAT_TABLE
    ! dry soil albedos: 1=vis, 2=nir
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ALBDRY_TABLE
	
	! energy out transfer variables
	! 2D Energy Flux and State Variables
    ! surface radiative temperature [K]
	real(kind=kind_noahmp), allocatable, dimension(:) :: TSK
    ! sensible heat flux [W m-2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: HFX
    ! latent heat flux [W m-2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: LH
    ! ground/snow heat flux [W m-2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: GRDFLX
    ! total grid albedo []
	real(kind=kind_noahmp), allocatable, dimension(:) :: ALBEDO
    ! surface bulk emissivity
	real(kind=kind_noahmp), allocatable, dimension(:) :: EMISS
    ! roughness length output to host
	real(kind=kind_noahmp), allocatable, dimension(:) :: Z0
    ! roughness length output to host
	real(kind=kind_noahmp), allocatable, dimension(:) :: ZNT
    ! 2m temperature of vegetation part
	real(kind=kind_noahmp), allocatable, dimension(:) :: T2MVXY
    ! 2m temperature of bare ground part
	real(kind=kind_noahmp), allocatable, dimension(:) :: T2MBXY
    ! 2m mixing ratio of vegetation part
	real(kind=kind_noahmp), allocatable, dimension(:) :: Q2MVXY
    ! 2m mixing ratio of bare ground part
	real(kind=kind_noahmp), allocatable, dimension(:) :: Q2MBXY
    ! surface radiative temperature (k)
	real(kind=kind_noahmp), allocatable, dimension(:) :: TRADXY
    ! Noah-MP vegetation fraction [-]
	real(kind=kind_noahmp), allocatable, dimension(:) :: FVEGXY
    ! total absorbed solar radiation (w/m2)
	real(kind=kind_noahmp), allocatable, dimension(:) :: FSAXY
    ! total net longwave rad (w/m2) [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: FIRAXY
    ! photosyn active energy by canopy (w/m2)
	real(kind=kind_noahmp), allocatable, dimension(:) :: APARXY
    ! solar rad absorbed by veg. (w/m2)
	real(kind=kind_noahmp), allocatable, dimension(:) :: SAVXY
    ! solar rad absorbed by ground (w/m2)
	real(kind=kind_noahmp), allocatable, dimension(:) :: SAGXY
    ! sunlit leaf stomatal resistance (s/m)
	real(kind=kind_noahmp), allocatable, dimension(:) :: RSSUNXY
    ! shaded leaf stomatal resistance (s/m)
	real(kind=kind_noahmp), allocatable, dimension(:) :: RSSHAXY
    ! between gap fraction
	real(kind=kind_noahmp), allocatable, dimension(:) :: BGAPXY
    ! within gap fraction
	real(kind=kind_noahmp), allocatable, dimension(:) :: WGAPXY
    ! under canopy ground temperature[K]
	real(kind=kind_noahmp), allocatable, dimension(:) :: TGVXY
    ! bare ground temperature [K]
	real(kind=kind_noahmp), allocatable, dimension(:) :: TGBXY
    ! sensible heat exchange coefficient vegetated
	real(kind=kind_noahmp), allocatable, dimension(:) :: CHVXY
    ! sensible heat exchange coefficient bare-ground
	real(kind=kind_noahmp), allocatable, dimension(:) :: CHBXY
    ! veg ground sen. heat [w/m2]  [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: SHGXY
    ! canopy sen. heat [w/m2]  [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: SHCXY
    ! bare sensible heat [w/m2] [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: SHBXY
    ! veg ground evap. heat [w/m2] [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: EVGXY
    ! bare soil evaporation [w/m2] [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: EVBXY
    ! veg ground heat flux [w/m2] [+ to soil]
	real(kind=kind_noahmp), allocatable, dimension(:) :: GHVXY
    ! bare ground heat flux [w/m2] [+ to soil]
	real(kind=kind_noahmp), allocatable, dimension(:) :: GHBXY
    ! veg ground net LW rad. [w/m2] [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: IRGXY
    ! canopy net LW rad. [w/m2] [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: IRCXY
    ! bare net longwave rad. [w/m2] [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: IRBXY
    ! transpiration [w/m2] [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: TRXY
    ! canopy evaporation heat [w/m2] [+ to atm]
	real(kind=kind_noahmp), allocatable, dimension(:) :: EVCXY
    ! leaf exchange coefficient
	real(kind=kind_noahmp), allocatable, dimension(:) :: CHLEAFXY
    ! under canopy exchange coefficient
	real(kind=kind_noahmp), allocatable, dimension(:) :: CHUCXY
    ! veg 2m exchange coefficient
	real(kind=kind_noahmp), allocatable, dimension(:) :: CHV2XY
    ! bare 2m exchange coefficient
	real(kind=kind_noahmp), allocatable, dimension(:) :: CHB2XY
    ! Total stomatal resistance [s/m]
	real(kind=kind_noahmp), allocatable, dimension(:) :: RS
    ! latent heating from sprinkler evaporation (w/m2)
	real(kind=kind_noahmp), allocatable, dimension(:) :: IRRSPLH
    ! precipitation advected heat [W/m2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: PAHXY
    ! precipitation advected heat [W/m2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: PAHGXY
    ! precipitation advected heat [W/m2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: PAHBXY
    ! precipitation advected heat [W/m2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: PAHVXY
    ! accumulated heat flux through soil bottom per soil timestep [J/m2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: EFLXBXY
    ! energy content in soil relative to 273.16 [KJ/m2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: SOILENERGY
    ! energy content in snow relative to 273.16 [KJ/m2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: SNOWENERGY
    ! canopy heat storage change [W/m2]
	real(kind=kind_noahmp), allocatable, dimension(:) :: CANHSXY

	! 3D Multi-Layer Variables
    ! snow albedo (direct)
	real(kind=kind_noahmp), allocatable, dimension(:,:) :: ALBSNOWDIRXY
    ! snow albedo (diffuse)
	real(kind=kind_noahmp), allocatable, dimension(:,:) :: ALBSNOWDIFXY
    ! surface albedo (direct)
	real(kind=kind_noahmp), allocatable, dimension(:,:) :: ALBSFCDIRXY
    ! surface albedo (diffuse)
	real(kind=kind_noahmp), allocatable, dimension(:,:) :: ALBSFCDIFXY
	
!------------------------------------------------------------------------
! End 2D variables not used in WRF
!------------------------------------------------------------------------

    CHARACTER(LEN=256)                                     ::  MMINSL  = 'STAS'    ! soil classification

    ! SNICAR optics and spectral lookup data
    ! downward diffuse solar radiation spectral weights for wavelength band
    real(kind=kind_noahmp), allocatable, dimension(:)   :: flx_wgt_dif, flx_wgt_dir
    ! Mie single scatter albedos for hydrophillic BC
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_bc1, asm_prm_bc1, ext_cff_mss_bc1
    ! Mie single scatter albedos for hydrophobic BC
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_bc2, asm_prm_bc2, ext_cff_mss_bc2
    ! Mie single scatter albedos for hydrophillic OC
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_oc1, asm_prm_oc1, ext_cff_mss_oc1
    ! Mie single scatter albedos for hydrophobic OC
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_oc2, asm_prm_oc2, ext_cff_mss_oc2
    ! Mie single scatter albedos for dust species 1
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_dst1, asm_prm_dst1, ext_cff_mss_dst1
    ! Mie single scatter albedos for dust species 2
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_dst2, asm_prm_dst2, ext_cff_mss_dst2
    ! Mie single scatter albedos for dust species 3
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_dst3, asm_prm_dst3, ext_cff_mss_dst3
    ! Mie single scatter albedos for dust species 4
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_dst4, asm_prm_dst4, ext_cff_mss_dst4
    ! Mie single scatter albedos for dust species 5
    real(kind=kind_noahmp), allocatable, dimension(:)   :: ss_alb_dst5, asm_prm_dst5, ext_cff_mss_dst5
    ! Mie single scatter albedos for direct-beam ice
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ss_alb_snw_drc, asm_prm_snw_drc, ext_cff_mss_snw_drc
    ! Mie single scatter albedos for diffuse ice
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: ss_alb_snw_dfs, asm_prm_snw_dfs, ext_cff_mss_snw_dfs
	
	! water in transfer variables
	real(kind=kind_noahmp), allocatable, dimension(:)    ::  CANLIQXY            ! canopy-intercepted liquid water (mm)
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  CANICEXY            ! canopy-intercepted ice (mm)
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FWETXY              ! wetted or snowed fraction of the canopy (-)
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SNOW                ! snow water equivalent [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SNEQVOXY            ! snow mass at last time step(mm h2o)
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SNOWH               ! physical snow depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FIFRACT             ! flood irrigation fraction
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  IRWATFI             ! irrigation water amount [m] to be applied, Flood
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  MIFRACT             ! micro irrigation fraction
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  IRWATMI             ! irrigation water amount [m] to be applied, Micro
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SIFRACT             ! sprinkler irrigation fraction
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  IRWATSI             ! irrigation water amount [m] to be applied, Sprinkler
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ZWTXY               ! water table depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SMCWTDXY            ! groundwater storage [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  TD_FRACTION         ! tile drainage fraction
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  WAXY                ! water in the "aquifer" [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  WTXY                ! groundwater storage [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  WSLAKEXY            ! lake water storage [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  IRFRACT             ! irrigation fraction
    integer, allocatable, dimension(:)                   ::  IRNUMSI             ! irrigation event number, Sprinkler
    integer, allocatable, dimension(:)                   ::  IRNUMMI             ! irrigation event number, Micro
    integer, allocatable, dimension(:)                   ::  IRNUMFI             ! irrigation event number, Flood
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  SNICEXY             ! snow layer ice [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  SNLIQXY             ! snow layer liquid water [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  SH2O                ! volumetric liquid soil moisture [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  SMOIS               ! volumetric soil moisture [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  SMOISEQ             ! volumetric soil moisture [m3/m3]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FSATXY              ! saturated fraction of the grid (-)
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  WSURFXY             ! wetland water storage [mm]
	
	real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACC_QSEVAXY         ! accumulated soil surface evaporation [m/s * dt_soil/dt_main]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACC_QINSURXY        ! accumulated water flux into soil [m/s * dt_soil/dt_main]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACC_DWATERXY        ! accumulated snow,soil,canopy water change per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACC_PRCPXY          ! accumulated precipitation per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACC_ECANXY          ! accumulated net canopy evaporation per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACC_ETRANXY         ! accumulated transpiration per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACC_EDIRXY          ! accumulated net ground (soil/snow) evaporation per soil timestep [mm]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  ACC_ETRANIXY        ! accumualted transpiration rate within soil timestep [m/s * dt_soil/dt_main]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACC_GLAFLWXY        ! accumulated glacier excessive flow [mm] per soil timestep
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  SNFRXY              ! snow layer rate of snow freezing [mm/s]
	
    ! SNICAR snow-layer state. The CCPP block layout is (snow layer,column).
    ! snow layer effective grain radius [microns, m-6]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: SNRDSXY
    ! mass of hydrophillic Black Carbon in snow [kg/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: BCPHIXY, BCPHOXY, OCPHIXY, OCPHOXY
    ! mass of dust species 1 in snow [kg/m2]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: DUST1XY, DUST2XY, DUST3XY, DUST4XY, DUST5XY
    ! mass concentration of hydrophillic Black Carbon in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MassConcBCPHIXY, MassConcBCPHOXY
    ! mass concentration of hydrophillic Organic Carbon in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MassConcOCPHIXY, MassConcOCPHOXY
    ! mass concentration of dust species 1 in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MassConcDUST1XY, MassConcDUST2XY
    ! mass concentration of dust species 3 in snow [kg/kg]
    real(kind=kind_noahmp), allocatable, dimension(:,:) :: MassConcDUST3XY, MassConcDUST4XY, MassConcDUST5XY
    ! Snow aging parameters retrieved from lookup table [hour]
    real(kind=kind_noahmp), allocatable, dimension(:,:,:) :: snowage_tau, snowage_kappa, snowage_drdt0
    ! hydrophobic Black Carbon deposition [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:) :: DepBChydrophoXY, DepBChydrophiXY
    ! hydrophobic Organic Carbon deposition [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:) :: DepOChydrophoXY, DepOChydrophiXY
    ! dust species 1 deposition [kg m-2 s-1]
    real(kind=kind_noahmp), allocatable, dimension(:) :: DepDust1XY, DepDust2XY, DepDust3XY, DepDust4XY, DepDust5XY

	integer                                                ::  DRAIN_LAYER_OPT_TABLE     ! tile drainage layer
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  CH2OP_TABLE               ! maximum intercepted h2o per unit lai+sai (mm)
    real(kind=kind_noahmp)                                 ::  C2_SNOWCOMPACT_TABLE      ! overburden snow compaction parameter (m3/kg)
    real(kind=kind_noahmp)                                 ::  C3_SNOWCOMPACT_TABLE      ! snow desctructive metamorphism compaction parameter1 [1/s]
    real(kind=kind_noahmp)                                 ::  C4_SNOWCOMPACT_TABLE      ! snow desctructive metamorphism compaction parameter2 [1/k]
    real(kind=kind_noahmp)                                 ::  C5_SNOWCOMPACT_TABLE      ! snow desctructive metamorphism compaction parameter3
    real(kind=kind_noahmp)                                 ::  DM_SNOWCOMPACT_TABLE      ! upper Limit on destructive metamorphism compaction [kg/m3]
    real(kind=kind_noahmp)                                 ::  ETA0_SNOWCOMPACT_TABLE    ! snow viscosity coefficient [kg-s/m2]
    real(kind=kind_noahmp)                                 ::  SNOWCOMPACTm_AR24_TABLE   ! snow compaction m parameter for linear sfc temp fitting from AR24
    real(kind=kind_noahmp)                                 ::  SNOWCOMPACTb_AR24_TABLE   ! snow compaction b parameter for linear sfc temp fitting from AR24
    real(kind=kind_noahmp)                                 ::  SNOWCOMPACT_P1_AR24_TABLE ! lower constrain for SnowCompactBurdenFac for high pressure bin from AR24
    real(kind=kind_noahmp)                                 ::  SNOWCOMPACT_P2_AR24_TABLE ! lower constrain for SnowCompactBurdenFac for mid pressure bin from AR24
    real(kind=kind_noahmp)                                 ::  SNOWCOMPACT_P3_AR24_TABLE ! lower constrain for SnowCompactBurdenFac for low pressure bin from AR24
    real(kind=kind_noahmp)                                 ::  SCFm1_AR25_TABLE          ! m1 parameter for ground SCF from AR2025
    real(kind=kind_noahmp)                                 ::  SCFm2_AR25_TABLE          ! m2 parameter for ground SCF from AR2025
    real(kind=kind_noahmp)                                 ::  SCfac1_AR25_TABLE         ! SCfac1 parameter for ground SCF from AR2025
    real(kind=kind_noahmp)                                 ::  SCfac2_AR25_TABLE         ! SCfac2 parameter for ground SCF from AR2025
    real(kind=kind_noahmp)                                 ::  SNOWCOMPACT_Up_AR24_TABLE ! upper constraint on SnowCompactBurdenFac from AR24
    real(kind=kind_noahmp)                                 ::  SNLIQMAXFRAC_TABLE        ! maximum liquid water fraction in snow
    real(kind=kind_noahmp)                                 ::  SSI_TABLE                 ! liquid water holding capacity for snowpack (m3/m3) (0.03)
    real(kind=kind_noahmp)                                 ::  SNOW_RET_FAC_TABLE        ! snowpack water release timescale factor (1/s)
    real(kind=kind_noahmp)                                 ::  FIRTFAC_TABLE             ! flood application rate factor
    real(kind=kind_noahmp)                                 ::  MICIR_RATE_TABLE          ! mm/h, micro irrigation rate
    real(kind=kind_noahmp)                                 ::  REFDK_TABLE               ! Parameter in the surface runoff parameterization
    real(kind=kind_noahmp)                                 ::  REFKDT_TABLE              ! Parameter in the surface runoff parameterization
    real(kind=kind_noahmp)                                 ::  FRZK_TABLE                ! Frozen ground parameter
    real(kind=kind_noahmp)                                 ::  TIMEAN_TABLE              ! gridcell mean topgraphic index (global mean)
    real(kind=kind_noahmp)                                 ::  FSATMX_TABLE              ! maximum surface saturated fraction (global mean)
    real(kind=kind_noahmp)                                 ::  ROUS_TABLE                ! specific yield [-] for Niu et al. 2007 groundwater scheme
    real(kind=kind_noahmp)                                 ::  CMIC_TABLE                ! microprore content (0.0-1.0), 0.0: close to free drainage
    real(kind=kind_noahmp)                                 ::  WSLMAX_TABLE              ! maximum lake water storage (mm)
    real(kind=kind_noahmp)                                 ::  SWEMAXGLA_TABLE           ! Maximum SWE allowed at glaciers (mm)
    integer                                                ::  IRR_HAR_TABLE             ! number of days before harvest date to stop irrigation 
    real(kind=kind_noahmp)                                 ::  IRR_LAI_TABLE             ! Minimum lai to trigger irrigation
    real(kind=kind_noahmp)                                 ::  IRR_MAD_TABLE             ! management allowable deficit (0-1)
    real(kind=kind_noahmp)                                 ::  FILOSS_TABLE              ! factor of flood irrigation loss
    real(kind=kind_noahmp)                                 ::  SPRIR_RATE_TABLE          ! mm/h, sprinkler irrigation rate
    real(kind=kind_noahmp)                                 ::  IRR_FRAC_TABLE            ! irrigation Fraction
    real(kind=kind_noahmp)                                 ::  IR_RAIN_TABLE             ! maximum precipitation to stop irrigation trigger
    real(kind=kind_noahmp)                                 ::  SNOWDEN_MAX_TABLE         ! maximum fresh snowfall density (kg/m3)
    real(kind=kind_noahmp)                                 ::  SWEMX_TABLE               ! new snow mass to fully cover old snow (mm)
    real(kind=kind_noahmp)                                 ::  PSIWLT_TABLE              ! soil metric potential for wilting point (m)
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  MFSNO_TABLE               ! snowmelt curve parameter
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  SCFFAC_TABLE              ! snow cover factor (m) (replace original hard-coded 2.5*z0 in SCF formulation)
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  BVIC_TABLE                ! VIC model infiltration parameter (-) for opt_run=6
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  AXAJ_TABLE                ! Xinanjiang: Tension water distribution inflection parameter [-] for opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  BXAJ_TABLE                ! Xinanjiang: Tension water distribution shape parameter [-] for opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  XXAJ_TABLE                ! Xinanjiang: Free water distribution shape parameter [-] for opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  BBVIC_TABLE               ! heterogeniety parameter for DVIC infiltration [-]
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  GDVIC_TABLE               ! mean capilary drive (m)
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  BDVIC_TABLE               ! VIC model infiltration parameter (-)
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  TD_DC_TABLE               ! tile drainage coefficient [mm/d]
    integer               , allocatable, dimension(:)      ::  TD_DEPTH_TABLE            ! tile drainage depth (layer number) from soil surface
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  TDSMC_FAC_TABLE           ! tile drainage soil moisture factor
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  TD_DCOEF_TABLE            ! tile drainage coefficient [mm/d]
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  TD_ADEPTH_TABLE           ! actual depth of impervious layer from land surface [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  KLAT_FAC_TABLE            ! hydraulic conductivity mutiplification factor
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  TD_DDRAIN_TABLE           ! tile drainage depth [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  TD_SPAC_TABLE             ! distance between two drain tubes or tiles [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  TD_RADI_TABLE             ! effective radius of drain tubes [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  TD_D_TABLE                ! depth to impervious layer from drain water level [m]
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  NROOT_TABLE               ! number of soil layers with root present
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  SLOPE_TABLE               ! slope factor for soil drainage
    real(kind=kind_noahmp)                                 ::  WCAP_TABLE                ! maximum surface wetland capacity
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  SMCMAX_TABLE              ! porosity, saturated value of soil moisture (volumetric)
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  SMCWLT_TABLE              ! wilting point soil moisture (volumetric)
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  SMCREF_TABLE              ! reference soil moisture (field capacity) (volumetric)
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  SMCDRY_TABLE              ! dry soil moisture threshold
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  DWSAT_TABLE               ! saturated soil hydraulic diffusivity
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  DKSAT_TABLE               ! saturated soil hydraulic conductivity
	integer                                                ::  SLCATS_TABLE              ! number of soil categories
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  BEXP_TABLE                ! soil B parameter
    real(kind=kind_noahmp), allocatable, dimension(:)      ::  PSISAT_TABLE              ! saturated soil matric potential
	
	real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  bexp_3D             ! C-H B exponent  
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  smcdry_3D           ! Soil Moisture Limit: Dry
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  smcwlt_3D           ! Soil Moisture Limit: Wilt
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  smcref_3D           ! Soil Moisture Limit: Reference
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  smcmax_3D           ! Soil Moisture Limit: Max
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  dksat_3D            ! Saturated Soil Conductivity
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  dwsat_3D            ! Saturated Soil Diffusivity
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  psisat_3D           ! Saturated Matric Potential
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  refdk_2D            ! Reference Soil Conductivity
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  refkdt_2D           ! Soil Infiltration Parameter
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  bvic_2d             ! VIC model infiltration parameter [-] opt_run=6
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  axaj_2D             ! Tension water distribution inflection parameter [-] opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  bxaj_2D             ! Tension water distribution shape parameter [-] opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  xxaj_2D             ! Free water distribution shape parameter [-] opt_run=7
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  bdvic_2d            ! VIC model infiltration parameter [-] opt_run=8
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  gdvic_2d            ! Mean Capillary Drive (m) for infiltration models opt_run=8
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  bbvic_2d            ! DVIC heterogeniety parameter for infiltration [-] opt_run=8
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  irr_frac_2D         ! irrigation Fraction
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  irr_har_2D          ! number of days before harvest date to stop irrigation 
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  irr_lai_2D          ! Minimum lai to trigger irrigation
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  irr_mad_2D          ! management allowable deficit (0-1)
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  filoss_2D           ! fraction of flood irrigation loss (0-1) 
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  sprir_rate_2D       ! mm/h, sprinkler irrigation rate
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  micir_rate_2D       ! mm/h, micro irrigation rate
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  firtfac_2D          ! flood application rate factor
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ir_rain_2D          ! maximum precipitation to stop irrigation trigger
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  KLAT_FAC            ! factor multiplier to hydraulic conductivity
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  TDSMC_FAC           ! factor multiplier to field capacity
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  TD_DC               ! drainage coefficient for simple
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  TD_DCOEF            ! drainge coefficient for Hooghoudt 
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  TD_DDRAIN           ! depth of drain
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  TD_RADI             ! tile radius
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  TD_SPAC             ! tile spacing
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FSATMX              ! maximum saturated fraction
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  WCAP                ! maximum wetland capacity [m]
    real(kind=kind_noahmp), allocatable, dimension(:,:)  ::  soilcomp            ! Soil sand and clay content [fraction]
	
	! water out transfer variables
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QFX                 ! latent heat flux [kg s-1 m-2]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SMSTAV              ! soil moisture avail. [not used]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SMSTOT              ! total soil water [mm][not used]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SFCRUNOFF           ! accumulated surface runoff [m]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  UDRUNOFF            ! accumulated sub-surface runoff [m]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SNOWC               ! snow cover fraction []
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  CANWAT              ! total canopy water + ice [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACSNOM              ! accumulated snow melt leaving pack
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ACSNOW              ! accumulated snow on grid
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QTDRAIN             ! tile drain discharge [mm]

    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QSNOWXY             ! snowfall on the ground [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QRAINXY             ! rainfall on the ground [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  DEEPRECHXY          ! groundwater storage [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  RECHXY              ! groundwater storage [mm]

    real(kind=kind_noahmp), allocatable, dimension(:)    ::  RUNSFXY             ! surface runoff [mm per soil timestep]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  RUNSBXY             ! subsurface runoff [mm per soil timestep]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ECANXY              ! evaporation of intercepted water (mm/s)
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  EDIRXY              ! soil surface evaporation rate (mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  ETRANXY             ! transpiration rate (mm/s)

    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QINTSXY             ! canopy intercepted snow [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QINTRXY             ! canopy intercepted rain [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QDRIPSXY            ! canopy dripping snow [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QDRIPRXY            ! canopy dripping rain [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QTHROSXY            ! canopy throughfall snow [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QTHRORXY            ! canopy throughfall rain [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QSNSUBXY            ! snowpack sublimation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QMELTXY             ! snowpack melting rate due to phase change [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QSNFROXY            ! snowpack frost rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QSUBCXY             ! canopy snow sublimation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QFROCXY             ! canopy snow frost rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QEVACXY             ! canopy water evaporation rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QDEWCXY             ! canopy water dew rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QFRZCXY             ! canopy water freezing rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QMELTCXY            ! canopy snow melting rate [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  QSNBOTXY            ! total water (melt+rain through snow) out of snowpack bottom [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  PONDINGXY           ! total surface ponding [mm]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  FPICEXY             ! fraction of ice in total precipitation
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  RAINLSM             ! total rain rate at the surface [mm/s]
    real(kind=kind_noahmp), allocatable, dimension(:)    ::  SNOWLSM             ! total snow rate at the surface [mm/s]

    real(kind=kind_noahmp), allocatable, dimension(:)    :: IRELOSS              ! loss of irrigation water to evaporation,sprinkler [m/timestep]
    real(kind=kind_noahmp), allocatable, dimension(:)    :: IRSIVOL              ! amount of irrigation by sprinkler (mm)
    real(kind=kind_noahmp), allocatable, dimension(:)    :: IRMIVOL              ! amount of irrigation by micro (mm)
    real(kind=kind_noahmp), allocatable, dimension(:)    :: IRFIVOL              ! amount of irrigation by micro (mm)  

    ! wetland model (OPT_WETLAND=1 or 2)

    ! optional parameters
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_a_TABLE      ! sand coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_b_TABLE      ! clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_c_TABLE      ! orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_d_TABLE      ! sand*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_e_TABLE      ! clay*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_f_TABLE      ! sand*clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500t_g_TABLE      ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500_a_TABLE       ! theta_1500t coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_1500_b_TABLE       ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_a_TABLE        ! sand coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_b_TABLE        ! clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_c_TABLE        ! orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_d_TABLE        ! sand*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_e_TABLE        ! clay*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_f_TABLE        ! sand*clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33t_g_TABLE        ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_33_a_TABLE         ! theta_33t*theta_33t coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33_b_TABLE         ! theta_33t coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_33_c_TABLE         ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_a_TABLE       ! sand coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_b_TABLE       ! clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_c_TABLE       ! orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_d_TABLE       ! sand*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_e_TABLE       ! clay*orgm coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_f_TABLE       ! sand*clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33t_g_TABLE       ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33_a_TABLE        ! theta_s33t coefficient
    real(kind=kind_noahmp)                                 :: sr2006_theta_s33_b_TABLE        ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_a_TABLE           ! sand coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_b_TABLE           ! clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_c_TABLE           ! theta_s33 coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_d_TABLE           ! sand*theta_s33 coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_e_TABLE           ! clay*theta_s33 coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_f_TABLE           ! sand*clay coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_et_g_TABLE           ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_psi_e_a_TABLE            ! psi_et*psi_et coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_e_b_TABLE            ! psi_et coefficient
    real(kind=kind_noahmp)                                 :: sr2006_psi_e_c_TABLE            ! constant adjustment
    real(kind=kind_noahmp)                                 :: sr2006_smcmax_a_TABLE           ! sand adjustment
    real(kind=kind_noahmp)                                 :: sr2006_smcmax_b_TABLE           ! constant adjustment
	
	! ! SNICAR scheme parameters
    real(kind=kind_noahmp)                                 :: DepBChydropho_TABLE      ! hydrophobic Black Carbon deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: DepBChydrophi_TABLE      ! hydrophillic Black Carbon deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: DepOChydropho_TABLE      ! hydrophobic Organic Carbon deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: DepOChydrophi_TABLE      ! hydrophillic Organic Carbon deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: DepDust1_TABLE           ! dust species 1 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: DepDust2_TABLE           ! dust species 2 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: DepDust3_TABLE           ! dust species 3 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: DepDust4_TABLE           ! dust species 4 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: DepDust5_TABLE           ! dust species 5 deposition [kg m-2 s-1], assume constant read from table
    real(kind=kind_noahmp)                                 :: SnowRadiusMin_TABLE      ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
    real(kind=kind_noahmp)                                 :: FreshSnowRadiusMax_TABLE ! maximum warm fresh snow effective radius [microns]
    real(kind=kind_noahmp)                                 :: SnowRadiusRefrz_TABLE    ! effective radius of re-frozen snow [microns]
    real(kind=kind_noahmp)                                 :: ScavEffMeltScale_TABLE   ! Scaling factor modifying scavenging factors for aerosol in meltwater (-)
    real(kind=kind_noahmp)                                 :: ScavEffMeltBCphi_TABLE   ! scavenging factor for hydrophillic BC inclusion in meltwater [frc]
    real(kind=kind_noahmp)                                 :: ScavEffMeltBCpho_TABLE   ! scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
    real(kind=kind_noahmp)                                 :: ScavEffMeltOCphi_TABLE   ! scavenging factor for hydrophillic OC inclusion in meltwater [frc]
    real(kind=kind_noahmp)                                 :: ScavEffMeltOCpho_TABLE   ! scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
    real(kind=kind_noahmp)                                 :: ScavEffMeltDust1_TABLE   ! scavenging factor for dust species 1 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)                                 :: ScavEffMeltDust2_TABLE   ! scavenging factor for dust species 2 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)                                 :: ScavEffMeltDust3_TABLE   ! scavenging factor for dust species 3 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)                                 :: ScavEffMeltDust4_TABLE   ! scavenging factor for dust species 4 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)                                 :: ScavEffMeltDust5_TABLE   ! scavenging factor for dust species 5 inclusion in meltwater  [frc]
    real(kind=kind_noahmp)                                 :: SnowRadiusMax_TABLE      ! maximum allowed snow effective radius [microns]
    real(kind=kind_noahmp)                                 :: SnowWetAgeC1Brun89_TABLE ! constant for liquid water grain growth [m3 s-1], from Brun89
    real(kind=kind_noahmp)                                 :: SnowWetAgeC2Brun89_TABLE ! Constant for liquid water grain growth [m3 s-1], from Brun89: corrected for LWC 
    real(kind=kind_noahmp)                                 :: SnowAgeScaleFac_TABLE    ! Arbitrary scaling factor applied to snow aging rate (-)
	
	!Biochem in transfers
    ! plant growth stage
	integer, allocatable, dimension(:)                   :: PGSXY
    ! planting day
	real(kind=kind_noahmp), allocatable, dimension(:)    :: PLANTING
    ! harvest day
	real(kind=kind_noahmp), allocatable, dimension(:)    :: HARVEST
    ! seasonal GDD
	real(kind=kind_noahmp), allocatable, dimension(:)    :: SEASON_GDD
	
    ! leaf mass [g/m2]
	real(kind=kind_noahmp), allocatable, dimension(:)    :: LFMASSXY
    ! mass of fine roots [g/m2]
	real(kind=kind_noahmp), allocatable, dimension(:)    :: RTMASSXY
    ! stem mass [g/m2]
	real(kind=kind_noahmp), allocatable, dimension(:)    :: STMASSXY
    ! mass of wood (incl. woody roots) [g/m2]
	real(kind=kind_noahmp), allocatable, dimension(:)    :: WOODXY
    ! stable carbon in deep soil [g/m2]
	real(kind=kind_noahmp), allocatable, dimension(:)    :: STBLCPXY
    ! short-lived carbon, shallow soil [g/m2]
	real(kind=kind_noahmp), allocatable, dimension(:)    :: FASTCPXY
    ! XING mass of grain!THREE
	real(kind=kind_noahmp), allocatable, dimension(:)    :: GRAINXY
    ! XINGgrowingdegressday
	real(kind=kind_noahmp), allocatable, dimension(:)    :: GDDXY
	
    ! Plant density [per ha] - used?
	real(kind=kind_noahmp), allocatable, dimension(:)    :: PLANTPOP_TABLE
    ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
	real(kind=kind_noahmp), allocatable, dimension(:)    :: IRRI_TABLE
    ! q10 for qe25
	real(kind=kind_noahmp), allocatable, dimension(:)      :: AQE_TABLE
    ! foliage nitrogen concentration when f(n)=1 (%)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: FOLNMX_TABLE
    ! quantum efficiency at 25c (umol co2 / umol photon)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: QE25_TABLE
    ! maximum rate of carboxylation at 25c (umol co2/m2/s)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: VCMX25_TABLE
    ! q10 for vcmx25
	real(kind=kind_noahmp), allocatable, dimension(:)      :: AVCMX_TABLE
    ! photosynthetic pathway: 0. = c4, 1. = c3
	real(kind=kind_noahmp), allocatable, dimension(:)      :: C3PSN_TABLE
    ! slope of conductance-to-photosynthesis relationship
	real(kind=kind_noahmp), allocatable, dimension(:)      :: MP_TABLE
    ! q10 for maintenance respiration
	real(kind=kind_noahmp), allocatable, dimension(:)      :: ARM_TABLE
    ! leaf maintenance respiration at 25c (umol co2/m2/s)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: RMF25_TABLE
    ! stem maintenance respiration at 25c (umol co2/kg bio/s)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: RMS25_TABLE
    ! root maintenance respiration at 25c (umol co2/kg bio/s)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: RMR25_TABLE
    ! wood to non-wood ratio
	real(kind=kind_noahmp), allocatable, dimension(:)      :: WRRAT_TABLE
    ! wood pool (switch 1 or 0) depending on woody or not [-]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: WDPOOL_TABLE
    ! leaf turnover [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: LTOVRC_TABLE
    ! characteristic T for leaf freezing [K]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: TDLEF_TABLE
    ! coeficient for leaf stress death [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: DILEFW_TABLE
    ! coeficient for leaf stress death [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: DILEFC_TABLE
    ! fraction of growth respiration !original was 0.3
	real(kind=kind_noahmp), allocatable, dimension(:)      :: FRAGR_TABLE
    ! microbial respiration parameter (umol co2 /kg c/ s)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: MRP_TABLE
    ! minimum temperature for photosynthesis (k)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: TMIN_TABLE
    ! single-side leaf area per Kg [m2/kg]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: SLA_TABLE
    ! minimum stem area index [m2/m2]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: XSAMIN_TABLE
    ! parameter for present wood allocation [-]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: BF_TABLE
    ! water stress coeficient [-]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: WSTRC_TABLE
    ! minimum leaf area index [m2/m2]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: LAIMIN_TABLE
    ! root turnover coefficient [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: RTOVRC_TABLE
    ! wood respiration coeficient [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: RSWOODC_TABLE
	
    ! Planting date
	integer, allocatable, dimension(:)                     :: PLTDAY_TABLE
    ! Harvest date
	integer, allocatable, dimension(:)                     :: HSDAY_TABLE
	
    ! foliage nitrogen concentration when f(n)=1 (%)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: FOLNMXI_TABLE
    ! quantum efficiency at 25c (umol co2 / umol photon)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: QE25I_TABLE
    ! maximum rate of carboxylation at 25c (umol co2/m2/s)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: VCMX25I_TABLE
    ! q10 for vcmx25
	real(kind=kind_noahmp), allocatable, dimension(:)      :: AVCMXI_TABLE
    ! photosynthetic pathway: 0. = c4, 1. = c3 ! Zhe Zhang 2020-07-03
	real(kind=kind_noahmp), allocatable, dimension(:)      :: C3PSNI_TABLE
    ! slope of conductance-to-photosynthesis relationship
	real(kind=kind_noahmp), allocatable, dimension(:)      :: MPI_TABLE

    ! q10 for maintainance respiration
	real(kind=kind_noahmp), allocatable, dimension(:)      :: Q10MR_TABLE
    ! leaf maintenance respiration at 25C [umol CO2/m2/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: LFMR25_TABLE
    ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: STMR25_TABLE
    ! root maintenance respiration at 25C [umol CO2/kg bio/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: RTMR25_TABLE

    ! fraction of growth respiration
	real(kind=kind_noahmp), allocatable, dimension(:)      :: FRA_GR_TABLE
    ! characteristic T for leaf freezing [K]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: LEFREEZ_TABLE
    ! leaf area per living leaf biomass [m2/kg]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: BIO2LAI_TABLE

    ! Base temperature for GDD accumulation [C]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDTBASE_TABLE
    ! Upper temperature for GDD accumulation [C]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDTCUT_TABLE
    ! GDD from seeding to emergence
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS1_TABLE
    ! GDD from seeding to initial vegetative
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS2_TABLE
    ! GDD from seeding to post vegetative
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS3_TABLE
    ! GDD from seeding to intial reproductive
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS4_TABLE
    ! GDD from seeding to pysical maturity
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GDDS5_TABLE

    ! Fraction of incoming solar radiation to photosynthetically active radiation
	real(kind=kind_noahmp), allocatable, dimension(:)      :: I2PAR_TABLE
    ! Minimum temperature for CO2 assimulation [C]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: TASSIM0_TABLE
    ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: TASSIM1_TABLE
    ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: TASSIM2_TABLE
    ! reference maximum CO2 assimulation rate
	real(kind=kind_noahmp), allocatable, dimension(:)      :: AREF_TABLE
    ! light extinction coefficient
	real(kind=kind_noahmp), allocatable, dimension(:)      :: K_TABLE
    ! initial light use efficiency
	real(kind=kind_noahmp), allocatable, dimension(:)      :: EPSI_TABLE
    ! CO2 assimulation reduction factor(0-1) (caused by non-modeled part, pest,weeds)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: PSNRF_TABLE

    ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GRAINMR25_TABLE
	
    ! coeficient for temperature leaf stress death [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: DILE_FC_TABLE
    ! coeficient for water leaf stress death [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: DILE_FW_TABLE

    ! fraction of carbohydrate translocation from leaf to grain
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: LFCT_TABLE
    ! fraction of carbohydrate translocation from stem to grain
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: STCT_TABLE
    ! fraction of carbohydrate translocation from root to grain
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: RTCT_TABLE

    ! fraction of carbohydrate flux to leaf
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: LFPT_TABLE
    ! fraction of carbohydrate flux to stem
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: STPT_TABLE
    ! fraction of carbohydrate flux to root
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: RTPT_TABLE
    ! fraction of carbohydrate flux to grain
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: GRAINPT_TABLE

    ! fraction of leaf turnover [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: LF_OVRC_TABLE
    ! fraction of stem turnover [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: ST_OVRC_TABLE
    ! fraction of root tunrover [1/s]
	real(kind=kind_noahmp), allocatable, dimension(:,:)    :: RT_OVRC_TABLE
	
	!Biochem out transfers
    ! net ecosys exchange (g/m2/s CO2)
	real(kind=kind_noahmp), allocatable, dimension(:)      :: NEEXY
    ! gross primary assimilation [g/m2/s C]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: GPPXY
    ! net primary productivity [g/m2/s C]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: NPPXY
    ! total photosynthesis (umol co2/m2/s) [+]
	real(kind=kind_noahmp), allocatable, dimension(:)      :: PSNXY
	
	! Read namelist variables
    ! soil time step (s) (default=0: same as main NoahMP timstep)
	real(kind=kind_noahmp)                                 ::  soiltstep
    ! Legacy HRLDAS LAI update setting; unused by CCPP.
    logical                                                ::  update_lai
    ! Legacy HRLDAS vegetation update setting; unused by CCPP.
    logical                                                ::  update_veg

    ! Legacy HRLDAS forcing directory; unused by CCPP.
    character(len=256)                                     ::  indir
    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  forcing_timestep
    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  noah_timestep
    ! Legacy HRLDAS simulation start-time component; unused by CCPP.
    integer                                                ::  start_year
    ! Legacy HRLDAS simulation start-time component; unused by CCPP.
    integer                                                ::  start_month
    ! Legacy HRLDAS simulation start-time component; unused by CCPP.
    integer                                                ::  start_day
    ! Legacy HRLDAS simulation start-time component; unused by CCPP.
    integer                                                ::  start_hour
    ! Legacy HRLDAS simulation start-time component; unused by CCPP.
    integer                                                ::  start_min
    ! Legacy HRLDAS output directory; unused by CCPP.
    character(len=256)                                     ::  outdir
    ! Legacy HRLDAS driver setting; unused by CCPP.
    character(len=256)                                     ::  restart_filename_requested
    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  restart_frequency_hours
    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  output_timestep
    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  spinup_loops

    ! =0: default output; >0 include additional output
    integer                                                ::  noahmp_output
    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  split_output_count
    ! Legacy HRLDAS driver setting; unused by CCPP.
    logical                                                ::  skip_first_output
    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  khour
    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  kday
    ! Legacy HRLDAS driver setting; unused by CCPP.
    real(kind=kind_noahmp)                                 ::  zlvl

    ! Legacy HRLDAS driver setting; unused by CCPP.
    integer                                                ::  use_wudapt_lcz
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_ndm
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_ng
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_nwr
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_ngb
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_nf
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_nz
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_nbui
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_hi
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    integer                                                ::  num_urban_ngr
    ! Legacy HRLDAS urban-grid dimension; unused by CCPP.
    real(kind=kind_noahmp)                                 ::  urban_atmosphere_thickness
    ! atmospheric levels including ZLVL for BEP/BEM models
    integer                                                ::  num_urban_atmosphere

    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_T
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_Q
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_U
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_V
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_P
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_LW
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_SW
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_PR
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_SN
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_DirFrac
    ! Legacy HRLDAS forcing variable name; unused by CCPP.
    character(len=256)                                     ::  forcing_name_VisFrac

    ! Legacy HRLDAS setup filename; unused by CCPP.
    character(len=256)                                     ::  hrldas_setup_file
    ! Legacy HRLDAS spatial input filename; unused by CCPP.
    character(len=256)                                     ::  spatial_filename
    ! Legacy HRLDAS vegetation input template; unused by CCPP.
    character(len=256)                                     ::  external_veg_filename_template
    ! Legacy HRLDAS LAI input template; unused by CCPP.
    character(len=256)                                     ::  external_lai_filename_template
    ! Legacy HRLDAS agricultural input filename; unused by CCPP.
    character(len=256)                                     ::  agdata_flnm
    ! Legacy HRLDAS tile-drainage input filename; unused by CCPP.
    character(len=256)                                     ::  tdinput_flnm
    ! SNICAR filename for optics parameters
    character(len=256)                                     ::  snicar_optic_flnm
    ! SNICAR filename for snow aging parameters
    character(len=256)                                     ::  snicar_age_flnm

    ! Legacy HRLDAS domain X lower bound; unused by CCPP.
    integer                                                ::  xstart
    ! Legacy HRLDAS domain Y lower bound; unused by CCPP.
    integer                                                ::  ystart
    ! Legacy HRLDAS domain X upper bound; unused by CCPP.
    integer                                                ::  xend
    ! Legacy HRLDAS domain Y upper bound; unused by CCPP.
    integer                                                ::  yend
    ! Legacy HRLDAS maximum soil-layer count; unused by CCPP.
    integer                                                ::  MAX_SOIL_LEVELS

    ! forcing variable for hydrophilic black carbon deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_BCPHI
    ! forcing variable for hydrophobic black carbon deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_BCPHO
    ! forcing variable for hydrophilic organic carbon deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_OCPHI
    ! forcing variable for hydrophobic organic carbon deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_OCPHO
    ! forcing variable for dust size bin 1 deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_DUST1
    ! forcing variable for dust size bin 2 deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_DUST2
    ! forcing variable for dust size bin 3 deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_DUST3
    ! forcing variable for dust size bin 4 deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_DUST4
    ! forcing variable for dust size bin 5 deposition flux [kg/m2/s]
    character(len=256)                                     ::  forcing_name_DUST5
	
  end type NoahmpIO_type

end module NoahmpIOVarType
