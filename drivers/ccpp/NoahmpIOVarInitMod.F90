! Harsh Kamath (NOAA GSL/CIRES), September 2026


module NoahmpIOVarInitMod

  use Machine
  use NoahmpIOVarType, only : NoahmpIO_type

  implicit none

contains

!=== initialize with default values

  subroutine NoahmpIOVarInitDefault(NoahmpIO, capacity)

    implicit none

    type(NoahmpIO_type), intent(inout) :: NoahmpIO
	
	integer, intent(in) :: capacity
	
	associate(                               &
              KMS     =>  NoahmpIO%KMS      ,&
              KME     =>  NoahmpIO%KME      ,&
              NSOIL   =>  NoahmpIO%NSOIL    ,&
              NSNOW   =>  NoahmpIO%NSNOW    ,&
              NUMRAD  =>  NoahmpIO%NUMRAD    &
        )
	
	
	NoahmpIO%capacity = capacity
    NoahmpIO%ncol     = 0
    if (.not. allocated(NoahmpIO%source_i)) allocate(NoahmpIO%source_i(1:capacity))
    if (.not. allocated(NoahmpIO%source_j)) allocate(NoahmpIO%source_j(1:capacity))
   
	! Configuration and column geometry (following the WRF driver order)
	if ( .not. allocated (NoahmpIO%DX)        ) allocate ( NoahmpIO%DX         (1:capacity) )
	if ( .not. allocated (NoahmpIO%DY)        ) allocate ( NoahmpIO%DY         (1:capacity) )
	if ( .not. allocated (NoahmpIO%COSZEN)    ) allocate ( NoahmpIO%COSZEN     (1:capacity) ) ! cosine zenith angle
	if ( .not. allocated (NoahmpIO%XLAT)      ) allocate ( NoahmpIO%XLAT       (1:capacity) ) ! latitude [radians] 
	if ( .not. allocated (NoahmpIO%ZSOIL)     ) allocate ( NoahmpIO%ZSOIL      (1:NSOIL)                         ) ! depth to soil interfaces [m] 
	if ( .not. allocated (NoahmpIO%IVGTYP)    ) allocate ( NoahmpIO%IVGTYP     (1:capacity) ) ! vegetation type
	if ( .not. allocated (NoahmpIO%ISLTYP)    ) allocate ( NoahmpIO%ISLTYP     (1:capacity) ) ! soil type
	if ( .not. allocated (NoahmpIO%SOILCOL)   ) allocate ( NoahmpIO%SOILCOL    (1:capacity) ) ! soil color type
	if ( .not. allocated (NoahmpIO%SLOPETYP)  ) allocate ( NoahmpIO%SLOPETYP   (1:capacity) ) ! slope type
	if ( .not. allocated (NoahmpIO%XICE)      ) allocate ( NoahmpIO%XICE       (1:capacity) ) ! fraction of grid that is seaice
	if ( .not. allocated (NoahmpIO%ICE)       ) allocate ( NoahmpIO%ICE        (1:capacity) ) ! sea-ice/glacier flag
	if ( .not. allocated (NoahmpIO%DZ8W)      ) allocate ( NoahmpIO%DZ8W       (kms:kme,1:capacity) ) ! thickness of atmo layers [m]

	! Allocate all four soil-class inputs; the selected soil option uses what it needs.
	   if ( .not. allocated (NoahmpIO%soilcl1) ) allocate ( NoahmpIO%soilcl1  (1:capacity) ) ! Soil texture class with depth
	   if ( .not. allocated (NoahmpIO%soilcl2) ) allocate ( NoahmpIO%soilcl2  (1:capacity) ) ! Soil texture class with depth
	   if ( .not. allocated (NoahmpIO%soilcl3) ) allocate ( NoahmpIO%soilcl3  (1:capacity) ) ! Soil texture class with depth
	   if ( .not. allocated (NoahmpIO%soilcl4) ) allocate ( NoahmpIO%soilcl4  (1:capacity) ) ! Soil texture class with depth

	if ( .not. allocated (NoahmpIO%ISNOWXY)   ) allocate ( NoahmpIO%ISNOWXY    (1:capacity) ) ! actual no. of snow layers
	if ( .not. allocated (NoahmpIO%ZSNSOXY)   ) allocate ( NoahmpIO%ZSNSOXY    (-NSNOW+1:NSOIL,1:capacity) ) ! snow layer depth [m]
	if ( .not. allocated (NoahmpIO%CROPCAT)   ) allocate ( NoahmpIO%CROPCAT    (1:capacity) )
	
	NoahmpIO%DX              = undefined_real
	NoahmpIO%DY              = undefined_real
	NoahmpIO%ICE             = undefined_int
	NoahmpIO%IVGTYP          = undefined_int
	NoahmpIO%ISLTYP          = undefined_int
	NoahmpIO%SOILCOL         = undefined_int
	NoahmpIO%ISNOWXY         = undefined_int
	NoahmpIO%COSZEN          = undefined_real
	NoahmpIO%XLAT            = undefined_real
	NoahmpIO%DZ8W            = undefined_real
	NoahmpIO%ZSOIL           = undefined_real
	NoahmpIO%XICE            = undefined_real

	NoahmpIO%ZSNSOXY         = undefined_real
	NoahmpIO%CROPCAT         = undefined_int

	   NoahmpIO%SOILCL1      = undefined_real
	   NoahmpIO%SOILCL2      = undefined_real
	   NoahmpIO%SOILCL3      = undefined_real
	   NoahmpIO%SOILCL4      = undefined_real

	NoahmpIO%XICE              = 0.0      ! fraction of grid that is seaice
	NoahmpIO%SLOPETYP          = 1        ! soil parameter slope type
	NoahmpIO%soil_update_steps = 1        ! number of model time step to update soil proces
	NoahmpIO%calculate_soil    = .false.  ! index for if do soil process
	

   ! Forcing inputs
	
	if ( .not. allocated (NoahmpIO%T_PHY)     ) allocate ( NoahmpIO%T_PHY      (kms:kme,1:capacity))  ! 3D atmospheric temperature valid at mid-levels [K]
    ! CCPP specific humidity at the lowest model layer [kg/kg].
    if ( .not. allocated (NoahmpIO%QV_CURR)   ) &
       allocate ( NoahmpIO%QV_CURR(kms:kme,1:capacity))
    if ( .not. allocated (NoahmpIO%U_PHY)     ) allocate ( NoahmpIO%U_PHY      (kms:kme,1:capacity))  ! 3D U wind component [m/s]
    if ( .not. allocated (NoahmpIO%V_PHY)     ) allocate ( NoahmpIO%V_PHY      (kms:kme,1:capacity))  ! 3D V wind component [m/s]
	
	if ( .not. allocated (NoahmpIO%SWDOWN)    ) allocate ( NoahmpIO%SWDOWN     (1:capacity) ) ! solar down at surface [W m-2]
	if ( .not. allocated (NoahmpIO%GLW)       ) allocate ( NoahmpIO%GLW        (1:capacity) ) ! longwave down at surface [W m-2]
	if ( .not. allocated (NoahmpIO%PS)        ) allocate ( NoahmpIO%PS         (1:capacity) ) ! surface air pressure [Pa]
	if ( .not. allocated (NoahmpIO%PRSL1)     ) allocate ( NoahmpIO%PRSL1      (1:capacity) ) ! lowest-layer pressure [Pa]
    if ( .not. allocated (NoahmpIO%PBLH)      ) allocate ( NoahmpIO%PBLH       (1:capacity) ) ! PBL height [m]
    if ( .not. allocated (NoahmpIO%RAINBL)    ) &
       allocate ( NoahmpIO%RAINBL(1:capacity) ) ! total rate [mm/s]
    if ( .not. allocated (NoahmpIO%SR)        ) &
       allocate ( NoahmpIO%SR(1:capacity) ) ! frozen fraction/type [-]
    ! CCPP convective precipitation rate [mm/s].
    if ( .not. allocated (NoahmpIO%MP_RAINC)  ) &
       allocate ( NoahmpIO%MP_RAINC(1:capacity) )
    ! CCPP explicit/nonconvective precipitation rate [mm/s].
    if ( .not. allocated (NoahmpIO%MP_RAINNC) ) &
       allocate ( NoahmpIO%MP_RAINNC(1:capacity) )
    ! Transitional shallow-convective precipitation rate [mm/s].
    if ( .not. allocated (NoahmpIO%MP_SHCV)   ) &
       allocate ( NoahmpIO%MP_SHCV(1:capacity) )
    if ( .not. allocated (NoahmpIO%MP_SNOW)   ) &
       allocate ( NoahmpIO%MP_SNOW(1:capacity) )
    if ( .not. allocated (NoahmpIO%MP_GRAUP)  ) &
       allocate ( NoahmpIO%MP_GRAUP(1:capacity) )
    if ( .not. allocated (NoahmpIO%MP_HAIL)   ) &
       allocate ( NoahmpIO%MP_HAIL(1:capacity) )
	
	if ( .not. allocated (NoahmpIO%TMN)       ) allocate ( NoahmpIO%TMN        (1:capacity) )
	if ( .not. allocated (NoahmpIO%RadSwDirFrac)) allocate ( NoahmpIO%RadSwDirFrac    (1:capacity) )
	if ( .not. allocated (NoahmpIO%RadSwVisFrac)) allocate ( NoahmpIO%RadSwVisFrac    (1:capacity) )
	
	! forcing out transfers
	if ( .not. allocated (NoahmpIO%FORCTLSM)    ) allocate ( NoahmpIO%FORCTLSM     (1:capacity) )
    if ( .not. allocated (NoahmpIO%FORCQLSM)    ) allocate ( NoahmpIO%FORCQLSM     (1:capacity) )
    if ( .not. allocated (NoahmpIO%FORCPLSM)    ) allocate ( NoahmpIO%FORCPLSM     (1:capacity) )
    if ( .not. allocated (NoahmpIO%FORCZLSM)    ) allocate ( NoahmpIO%FORCZLSM     (1:capacity) )
    if ( .not. allocated (NoahmpIO%FORCWLSM)    ) allocate ( NoahmpIO%FORCWLSM     (1:capacity) )
	
	
    !-------------------------------------------------------------------
    ! Initialize variables with default values 
    !-------------------------------------------------------------------

	NoahmpIO%source_i             = undefined_int
	NoahmpIO%source_j             = undefined_int
	
	NoahmpIO%T_PHY           = undefined_real
    NoahmpIO%QV_CURR         = undefined_real
    NoahmpIO%U_PHY           = undefined_real
    NoahmpIO%V_PHY           = undefined_real

	NoahmpIO%SWDOWN          = undefined_real
	NoahmpIO%GLW             = undefined_real
	NoahmpIO%PS              = undefined_real
	NoahmpIO%PRSL1           = undefined_real
    NoahmpIO%PBLH            = 1000.0_kind_noahmp
	NoahmpIO%RAINBL          = undefined_real
    NoahmpIO%SR              = undefined_real
	
	NoahmpIO%MP_RAINC        = 0.0
    NoahmpIO%MP_RAINNC       = 0.0
    NoahmpIO%MP_SHCV         = 0.0
    NoahmpIO%MP_SNOW         = 0.0
    NoahmpIO%MP_GRAUP        = 0.0
    NoahmpIO%MP_HAIL         = 0.0
	

	NoahmpIO%FORCTLSM        = undefined_real
    NoahmpIO%FORCQLSM        = undefined_real
    NoahmpIO%FORCPLSM        = undefined_real
    NoahmpIO%FORCZLSM        = undefined_real
    NoahmpIO%FORCWLSM        = undefined_real

	! energy in
	! 2D/3D state, flux, and spatial arrays
    if ( .not. allocated (NoahmpIO%LAI)          ) allocate ( NoahmpIO%LAI          (1:capacity) )
    if ( .not. allocated (NoahmpIO%XSAIXY)       ) allocate ( NoahmpIO%XSAIXY       (1:capacity) )
    if ( .not. allocated (NoahmpIO%QSFC)         ) allocate ( NoahmpIO%QSFC         (1:capacity) )
    if ( .not. allocated (NoahmpIO%TGXY)         ) allocate ( NoahmpIO%TGXY         (1:capacity) )
    if ( .not. allocated (NoahmpIO%TVXY)         ) allocate ( NoahmpIO%TVXY         (1:capacity) )
    if ( .not. allocated (NoahmpIO%TAUSSXY)      ) allocate ( NoahmpIO%TAUSSXY      (1:capacity) )
    if ( .not. allocated (NoahmpIO%ALBOLDXY)     ) allocate ( NoahmpIO%ALBOLDXY     (1:capacity) )
    if ( .not. allocated (NoahmpIO%EAHXY)        ) allocate ( NoahmpIO%EAHXY        (1:capacity) )
    if ( .not. allocated (NoahmpIO%TAHXY)        ) allocate ( NoahmpIO%TAHXY        (1:capacity) )
    if ( .not. allocated (NoahmpIO%CHXY)         ) allocate ( NoahmpIO%CHXY         (1:capacity) )
    if ( .not. allocated (NoahmpIO%CMXY)         ) allocate ( NoahmpIO%CMXY         (1:capacity) )
    if ( .not. allocated (NoahmpIO%TSNOXY)       ) allocate ( NoahmpIO%TSNOXY       (-NSNOW+1:0,1:capacity) )
    if ( .not. allocated (NoahmpIO%TSLB)         ) allocate ( NoahmpIO%TSLB         (1:NSOIL,1:capacity) )
    if ( .not. allocated (NoahmpIO%ALBSOILDIRXY) ) allocate ( NoahmpIO%ALBSOILDIRXY (1:NUMRAD,1:capacity) )
    if ( .not. allocated (NoahmpIO%ALBSOILDIFXY) ) allocate ( NoahmpIO%ALBSOILDIFXY (1:NUMRAD,1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_SSOILXY)  ) allocate ( NoahmpIO%ACC_SSOILXY  (1:capacity) )
    if ( .not. allocated (NoahmpIO%GVFMAX)       ) allocate ( NoahmpIO%GVFMAX       (1:capacity) )
    if ( .not. allocated (NoahmpIO%VEGFRA)       ) allocate ( NoahmpIO%VEGFRA       (1:capacity) )

    ! Conditional spatial parameter allocation (OptSoilProperty = 4)
    !if ( NoahmpIO%IOPT_SOIL == 4 ) then
       if ( .not. allocated (NoahmpIO%quartz_3D) ) allocate ( NoahmpIO%quartz_3D (1:NSOIL,1:capacity) )
    !endif

    ! SNICAR scheme array allocations (OptSnowAlbedo = 3)
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       if ( .not. allocated (NoahmpIO%ss_alb_snw_drc)      ) allocate ( NoahmpIO%ss_alb_snw_drc      (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_snw_drc)     ) allocate ( NoahmpIO%asm_prm_snw_drc     (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_snw_drc) ) allocate ( NoahmpIO%ext_cff_mss_snw_drc (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_snw_dfs)      ) allocate ( NoahmpIO%ss_alb_snw_dfs      (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_snw_dfs)     ) allocate ( NoahmpIO%asm_prm_snw_dfs     (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_snw_dfs) ) allocate ( NoahmpIO%ext_cff_mss_snw_dfs (NoahmpIO%idx_Mie_snw_mx,NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_bc1)          ) allocate ( NoahmpIO%ss_alb_bc1          (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_bc1)         ) allocate ( NoahmpIO%asm_prm_bc1         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_bc1)     ) allocate ( NoahmpIO%ext_cff_mss_bc1     (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_bc2)          ) allocate ( NoahmpIO%ss_alb_bc2          (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_bc2)         ) allocate ( NoahmpIO%asm_prm_bc2         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_bc2)     ) allocate ( NoahmpIO%ext_cff_mss_bc2     (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_oc1)          ) allocate ( NoahmpIO%ss_alb_oc1          (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_oc1)         ) allocate ( NoahmpIO%asm_prm_oc1         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_oc1)     ) allocate ( NoahmpIO%ext_cff_mss_oc1     (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_oc2)          ) allocate ( NoahmpIO%ss_alb_oc2          (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_oc2)         ) allocate ( NoahmpIO%asm_prm_oc2         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_oc2)     ) allocate ( NoahmpIO%ext_cff_mss_oc2     (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst1)         ) allocate ( NoahmpIO%ss_alb_dst1         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst1)        ) allocate ( NoahmpIO%asm_prm_dst1        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst1)    ) allocate ( NoahmpIO%ext_cff_mss_dst1    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst2)         ) allocate ( NoahmpIO%ss_alb_dst2         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst2)        ) allocate ( NoahmpIO%asm_prm_dst2        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst2)    ) allocate ( NoahmpIO%ext_cff_mss_dst2    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst3)         ) allocate ( NoahmpIO%ss_alb_dst3         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst3)        ) allocate ( NoahmpIO%asm_prm_dst3        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst3)    ) allocate ( NoahmpIO%ext_cff_mss_dst3    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst4)         ) allocate ( NoahmpIO%ss_alb_dst4         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst4)        ) allocate ( NoahmpIO%asm_prm_dst4        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst4)    ) allocate ( NoahmpIO%ext_cff_mss_dst4    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ss_alb_dst5)         ) allocate ( NoahmpIO%ss_alb_dst5         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%asm_prm_dst5)        ) allocate ( NoahmpIO%asm_prm_dst5        (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%ext_cff_mss_dst5)    ) allocate ( NoahmpIO%ext_cff_mss_dst5    (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%flx_wgt_dir)         ) allocate ( NoahmpIO%flx_wgt_dir         (NoahmpIO%snicar_numrad_snw) )
       if ( .not. allocated (NoahmpIO%flx_wgt_dif)         ) allocate ( NoahmpIO%flx_wgt_dif         (NoahmpIO%snicar_numrad_snw) )
    endif

    ! Default value assignments
    NoahmpIO%VEGFRA        = undefined_real
    NoahmpIO%QSFC          = undefined_real
    NoahmpIO%TSLB          = undefined_real
    NoahmpIO%TVXY          = undefined_real
    NoahmpIO%TGXY          = undefined_real
    NoahmpIO%EAHXY         = undefined_real
    NoahmpIO%TAHXY         = undefined_real
    NoahmpIO%CMXY          = undefined_real
    NoahmpIO%CHXY          = undefined_real
    NoahmpIO%ALBOLDXY      = undefined_real
    NoahmpIO%TSNOXY        = undefined_real
    NoahmpIO%LAI           = undefined_real
    NoahmpIO%XSAIXY        = undefined_real
    NoahmpIO%ALBSOILDIRXY  = 0.0
    NoahmpIO%ALBSOILDIFXY  = 0.0
    NoahmpIO%TAUSSXY       = 0.0
    NoahmpIO%ACC_SSOILXY   = 0.0
    NoahmpIO%GVFMAX        = undefined_real

    if ( NoahmpIO%IOPT_ALB == 3 ) then
       NoahmpIO%ss_alb_snw_drc = undefined_real
       NoahmpIO%asm_prm_snw_drc = undefined_real
       NoahmpIO%ext_cff_mss_snw_drc = undefined_real
       NoahmpIO%ss_alb_snw_dfs = undefined_real
       NoahmpIO%asm_prm_snw_dfs = undefined_real
       NoahmpIO%ext_cff_mss_snw_dfs = undefined_real
       NoahmpIO%ss_alb_bc1 = undefined_real
       NoahmpIO%asm_prm_bc1 = undefined_real
       NoahmpIO%ext_cff_mss_bc1 = undefined_real
       NoahmpIO%ss_alb_bc2 = undefined_real
       NoahmpIO%asm_prm_bc2 = undefined_real
       NoahmpIO%ext_cff_mss_bc2 = undefined_real
       NoahmpIO%ss_alb_oc1 = undefined_real
       NoahmpIO%asm_prm_oc1 = undefined_real
       NoahmpIO%ext_cff_mss_oc1 = undefined_real
       NoahmpIO%ss_alb_oc2 = undefined_real
       NoahmpIO%asm_prm_oc2 = undefined_real
       NoahmpIO%ext_cff_mss_oc2 = undefined_real
       NoahmpIO%ss_alb_dst1 = undefined_real
       NoahmpIO%asm_prm_dst1 = undefined_real
       NoahmpIO%ext_cff_mss_dst1 = undefined_real
       NoahmpIO%ss_alb_dst2 = undefined_real
       NoahmpIO%asm_prm_dst2 = undefined_real
       NoahmpIO%ext_cff_mss_dst2 = undefined_real
       NoahmpIO%ss_alb_dst3 = undefined_real
       NoahmpIO%asm_prm_dst3 = undefined_real
       NoahmpIO%ext_cff_mss_dst3 = undefined_real
       NoahmpIO%ss_alb_dst4 = undefined_real
       NoahmpIO%asm_prm_dst4 = undefined_real
       NoahmpIO%ext_cff_mss_dst4 = undefined_real
       NoahmpIO%ss_alb_dst5 = undefined_real
       NoahmpIO%asm_prm_dst5 = undefined_real
       NoahmpIO%ext_cff_mss_dst5 = undefined_real
       NoahmpIO%flx_wgt_dir = undefined_real
       NoahmpIO%flx_wgt_dif = undefined_real
    endif
	
	! energy out transfer variables
	
	! 2D Variables (XSTART:XEND, YSTART:YEND)
	if ( .not. allocated (NoahmpIO%TSK)          ) allocate ( NoahmpIO%TSK          (1:capacity) )
	if ( .not. allocated (NoahmpIO%HFX)          ) allocate ( NoahmpIO%HFX          (1:capacity) )
	if ( .not. allocated (NoahmpIO%LH)           ) allocate ( NoahmpIO%LH           (1:capacity) )
	if ( .not. allocated (NoahmpIO%GRDFLX)       ) allocate ( NoahmpIO%GRDFLX       (1:capacity) )
	if ( .not. allocated (NoahmpIO%ALBEDO)       ) allocate ( NoahmpIO%ALBEDO       (1:capacity) )
	if ( .not. allocated (NoahmpIO%EMISS)        ) allocate ( NoahmpIO%EMISS        (1:capacity) )
	if ( .not. allocated (NoahmpIO%Z0)           ) allocate ( NoahmpIO%Z0           (1:capacity) )
	if ( .not. allocated (NoahmpIO%ZNT)          ) allocate ( NoahmpIO%ZNT          (1:capacity) )
	if ( .not. allocated (NoahmpIO%T2MVXY)       ) allocate ( NoahmpIO%T2MVXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%T2MBXY)       ) allocate ( NoahmpIO%T2MBXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%Q2MVXY)       ) allocate ( NoahmpIO%Q2MVXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%Q2MBXY)       ) allocate ( NoahmpIO%Q2MBXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%TRADXY)       ) allocate ( NoahmpIO%TRADXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%FVEGXY)       ) allocate ( NoahmpIO%FVEGXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%FSAXY)        ) allocate ( NoahmpIO%FSAXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%FIRAXY)       ) allocate ( NoahmpIO%FIRAXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%APARXY)       ) allocate ( NoahmpIO%APARXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%SAVXY)        ) allocate ( NoahmpIO%SAVXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%SAGXY)        ) allocate ( NoahmpIO%SAGXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%RSSUNXY)      ) allocate ( NoahmpIO%RSSUNXY      (1:capacity) )
	if ( .not. allocated (NoahmpIO%RSSHAXY)      ) allocate ( NoahmpIO%RSSHAXY      (1:capacity) )
	if ( .not. allocated (NoahmpIO%BGAPXY)       ) allocate ( NoahmpIO%BGAPXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%WGAPXY)       ) allocate ( NoahmpIO%WGAPXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%TGVXY)        ) allocate ( NoahmpIO%TGVXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%TGBXY)        ) allocate ( NoahmpIO%TGBXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%CHVXY)        ) allocate ( NoahmpIO%CHVXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%CHBXY)        ) allocate ( NoahmpIO%CHBXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%SHGXY)        ) allocate ( NoahmpIO%SHGXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%SHCXY)        ) allocate ( NoahmpIO%SHCXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%SHBXY)        ) allocate ( NoahmpIO%SHBXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%EVGXY)        ) allocate ( NoahmpIO%EVGXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%EVBXY)        ) allocate ( NoahmpIO%EVBXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%GHVXY)        ) allocate ( NoahmpIO%GHVXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%GHBXY)        ) allocate ( NoahmpIO%GHBXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%IRGXY)        ) allocate ( NoahmpIO%IRGXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%IRCXY)        ) allocate ( NoahmpIO%IRCXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%IRBXY)        ) allocate ( NoahmpIO%IRBXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%TRXY)         ) allocate ( NoahmpIO%TRXY         (1:capacity) )
	if ( .not. allocated (NoahmpIO%EVCXY)        ) allocate ( NoahmpIO%EVCXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%CHLEAFXY)     ) allocate ( NoahmpIO%CHLEAFXY     (1:capacity) )
	if ( .not. allocated (NoahmpIO%CHUCXY)       ) allocate ( NoahmpIO%CHUCXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%CHV2XY)       ) allocate ( NoahmpIO%CHV2XY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%CHB2XY)       ) allocate ( NoahmpIO%CHB2XY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%RS)           ) allocate ( NoahmpIO%RS           (1:capacity) )
	if ( .not. allocated (NoahmpIO%IRRSPLH)      ) allocate ( NoahmpIO%IRRSPLH      (1:capacity) )
	if ( .not. allocated (NoahmpIO%PAHXY)        ) allocate ( NoahmpIO%PAHXY        (1:capacity) )
	if ( .not. allocated (NoahmpIO%PAHGXY)       ) allocate ( NoahmpIO%PAHGXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%PAHBXY)       ) allocate ( NoahmpIO%PAHBXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%PAHVXY)       ) allocate ( NoahmpIO%PAHVXY       (1:capacity) )
	if ( .not. allocated (NoahmpIO%EFLXBXY)      ) allocate ( NoahmpIO%EFLXBXY      (1:capacity) )
	if ( .not. allocated (NoahmpIO%SOILENERGY)   ) allocate ( NoahmpIO%SOILENERGY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%SNOWENERGY)   ) allocate ( NoahmpIO%SNOWENERGY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%CANHSXY)      ) allocate ( NoahmpIO%CANHSXY      (1:capacity) )

	! 3D Multi-Layer Variables
	if ( .not. allocated (NoahmpIO%ALBSNOWDIRXY) ) allocate ( NoahmpIO%ALBSNOWDIRXY (1:NUMRAD, 1:capacity) )
	if ( .not. allocated (NoahmpIO%ALBSNOWDIFXY) ) allocate ( NoahmpIO%ALBSNOWDIFXY (1:NUMRAD, 1:capacity) )
	if ( .not. allocated (NoahmpIO%ALBSFCDIRXY)  ) allocate ( NoahmpIO%ALBSFCDIRXY  (1:NUMRAD, 1:capacity) )
	if ( .not. allocated (NoahmpIO%ALBSFCDIFXY)  ) allocate ( NoahmpIO%ALBSFCDIFXY  (1:NUMRAD, 1:capacity) )


	!=======================================================================
	! INITIAL VALUE ASSIGNMENTS
	!=======================================================================

	NoahmpIO%QSFC          = undefined_real
	NoahmpIO%TSK           = undefined_real
	NoahmpIO%TSLB          = undefined_real
	NoahmpIO%ALBEDO        = undefined_real
	NoahmpIO%TVXY          = undefined_real
	NoahmpIO%TGXY          = undefined_real
	NoahmpIO%EAHXY         = undefined_real
	NoahmpIO%TAHXY         = undefined_real
	NoahmpIO%CMXY          = undefined_real
	NoahmpIO%CHXY          = undefined_real
	NoahmpIO%ALBOLDXY      = undefined_real
	NoahmpIO%TSNOXY        = undefined_real
	NoahmpIO%LAI           = undefined_real
	NoahmpIO%XSAIXY        = undefined_real
	NoahmpIO%GRDFLX        = undefined_real
	NoahmpIO%HFX           = undefined_real
	NoahmpIO%LH            = undefined_real
	NoahmpIO%EMISS         = undefined_real
	NoahmpIO%T2MVXY        = undefined_real
	NoahmpIO%T2MBXY        = undefined_real
	NoahmpIO%Q2MVXY        = undefined_real
	NoahmpIO%Q2MBXY        = undefined_real
	NoahmpIO%TRADXY        = undefined_real
	NoahmpIO%FVEGXY        = undefined_real
	NoahmpIO%FSAXY         = undefined_real
	NoahmpIO%FIRAXY        = undefined_real
	NoahmpIO%APARXY        = undefined_real
	NoahmpIO%SAVXY         = undefined_real
	NoahmpIO%SAGXY         = undefined_real
	NoahmpIO%RSSUNXY       = undefined_real
	NoahmpIO%RSSHAXY       = undefined_real
	NoahmpIO%BGAPXY        = undefined_real
	NoahmpIO%WGAPXY        = undefined_real
	NoahmpIO%TGVXY         = undefined_real
	NoahmpIO%TGBXY         = undefined_real
	NoahmpIO%CHVXY         = undefined_real
	NoahmpIO%CHBXY         = undefined_real
	NoahmpIO%SHGXY         = undefined_real
	NoahmpIO%SHCXY         = undefined_real
	NoahmpIO%SHBXY         = undefined_real
	NoahmpIO%EVGXY         = undefined_real
	NoahmpIO%EVBXY         = undefined_real
	NoahmpIO%GHVXY         = undefined_real
	NoahmpIO%GHBXY         = undefined_real
	NoahmpIO%IRGXY         = undefined_real
	NoahmpIO%IRCXY         = undefined_real
	NoahmpIO%IRBXY         = undefined_real
	NoahmpIO%TRXY          = undefined_real
	NoahmpIO%EVCXY         = undefined_real
	NoahmpIO%CHLEAFXY      = undefined_real
	NoahmpIO%CHUCXY        = undefined_real
	NoahmpIO%CHV2XY        = undefined_real
	NoahmpIO%CHB2XY        = undefined_real
	NoahmpIO%RS            = undefined_real
	NoahmpIO%CANHSXY       = undefined_real
	NoahmpIO%Z0            = undefined_real
	NoahmpIO%ZNT           = undefined_real
	NoahmpIO%ALBSNOWDIRXY  = undefined_real
	NoahmpIO%ALBSNOWDIFXY  = undefined_real
	NoahmpIO%ALBSFCDIRXY   = undefined_real
	NoahmpIO%ALBSFCDIFXY   = undefined_real
	NoahmpIO%PAHXY         = undefined_real
	NoahmpIO%PAHGXY        = undefined_real
	NoahmpIO%PAHBXY        = undefined_real
	NoahmpIO%PAHVXY        = undefined_real
	NoahmpIO%EFLXBXY       = undefined_real
	NoahmpIO%SOILENERGY    = undefined_real
	NoahmpIO%SNOWENERGY    = undefined_real

	NoahmpIO%ALBSOILDIRXY  = 0.0
	NoahmpIO%ALBSOILDIFXY  = 0.0
	NoahmpIO%TAUSSXY       = 0.0
	NoahmpIO%ACC_SSOILXY   = 0.0
	NoahmpIO%IRRSPLH       = 0.0
	
	!water in transfers
	!if ( NoahmpIO%IOPT_SOIL > 1 ) then
       if ( .not. allocated (NoahmpIO%soilcomp)) allocate ( NoahmpIO%soilcomp (1:2*NSOIL,1:capacity) )
    !endif
    !if ( NoahmpIO%IOPT_SOIL == 4 ) then
       if ( .not. allocated (NoahmpIO%bexp_3d)      ) allocate ( NoahmpIO%bexp_3d       (1:NSOIL,1:capacity) )
       if ( .not. allocated (NoahmpIO%smcdry_3D)    ) allocate ( NoahmpIO%smcdry_3D     (1:NSOIL,1:capacity) )
       if ( .not. allocated (NoahmpIO%smcwlt_3D)    ) allocate ( NoahmpIO%smcwlt_3D     (1:NSOIL,1:capacity) )
       if ( .not. allocated (NoahmpIO%smcref_3D)    ) allocate ( NoahmpIO%smcref_3D     (1:NSOIL,1:capacity) )
       if ( .not. allocated (NoahmpIO%smcmax_3D)    ) allocate ( NoahmpIO%smcmax_3D     (1:NSOIL,1:capacity) )
       if ( .not. allocated (NoahmpIO%dksat_3D)     ) allocate ( NoahmpIO%dksat_3D      (1:NSOIL,1:capacity) )
       if ( .not. allocated (NoahmpIO%dwsat_3D)     ) allocate ( NoahmpIO%dwsat_3D      (1:NSOIL,1:capacity) )
       if ( .not. allocated (NoahmpIO%psisat_3D)    ) allocate ( NoahmpIO%psisat_3D     (1:NSOIL,1:capacity) )
	   if ( .not. allocated (NoahmpIO%refdk_2D)     ) allocate ( NoahmpIO%refdk_2D      (1:capacity) )
       if ( .not. allocated (NoahmpIO%refkdt_2D)    ) allocate ( NoahmpIO%refkdt_2D     (1:capacity) )
       if ( .not. allocated (NoahmpIO%irr_frac_2D)  ) allocate ( NoahmpIO%irr_frac_2D   (1:capacity) )
       if ( .not. allocated (NoahmpIO%irr_har_2D)   ) allocate ( NoahmpIO%irr_har_2D    (1:capacity) )
       if ( .not. allocated (NoahmpIO%irr_lai_2D)   ) allocate ( NoahmpIO%irr_lai_2D    (1:capacity) )
       if ( .not. allocated (NoahmpIO%irr_mad_2D)   ) allocate ( NoahmpIO%irr_mad_2D    (1:capacity) )
       if ( .not. allocated (NoahmpIO%filoss_2D)    ) allocate ( NoahmpIO%filoss_2D     (1:capacity) )
       if ( .not. allocated (NoahmpIO%sprir_rate_2D)) allocate ( NoahmpIO%sprir_rate_2D (1:capacity) )
       if ( .not. allocated (NoahmpIO%micir_rate_2D)) allocate ( NoahmpIO%micir_rate_2D (1:capacity) )
       if ( .not. allocated (NoahmpIO%firtfac_2D)   ) allocate ( NoahmpIO%firtfac_2D    (1:capacity) )
       if ( .not. allocated (NoahmpIO%ir_rain_2D)   ) allocate ( NoahmpIO%ir_rain_2D    (1:capacity) )
       if ( .not. allocated (NoahmpIO%bvic_2D)      ) allocate ( NoahmpIO%bvic_2D       (1:capacity) )
       if ( .not. allocated (NoahmpIO%axaj_2D)      ) allocate ( NoahmpIO%axaj_2D       (1:capacity) )
       if ( .not. allocated (NoahmpIO%bxaj_2D)      ) allocate ( NoahmpIO%bxaj_2D       (1:capacity) )
       if ( .not. allocated (NoahmpIO%xxaj_2D)      ) allocate ( NoahmpIO%xxaj_2D       (1:capacity) )
       if ( .not. allocated (NoahmpIO%bdvic_2D)     ) allocate ( NoahmpIO%bdvic_2D      (1:capacity) )
       if ( .not. allocated (NoahmpIO%gdvic_2D)     ) allocate ( NoahmpIO%gdvic_2D      (1:capacity) )
       if ( .not. allocated (NoahmpIO%bbvic_2D)     ) allocate ( NoahmpIO%bbvic_2D      (1:capacity) )
       if ( .not. allocated (NoahmpIO%KLAT_FAC)     ) allocate ( NoahmpIO%KLAT_FAC      (1:capacity) )
       if ( .not. allocated (NoahmpIO%TDSMC_FAC)    ) allocate ( NoahmpIO%TDSMC_FAC     (1:capacity) )
       if ( .not. allocated (NoahmpIO%TD_DC)        ) allocate ( NoahmpIO%TD_DC         (1:capacity) )
       if ( .not. allocated (NoahmpIO%TD_DCOEF)     ) allocate ( NoahmpIO%TD_DCOEF      (1:capacity) )
       if ( .not. allocated (NoahmpIO%TD_DDRAIN)    ) allocate ( NoahmpIO%TD_DDRAIN     (1:capacity) )
       if ( .not. allocated (NoahmpIO%TD_RADI)      ) allocate ( NoahmpIO%TD_RADI       (1:capacity) )
       if ( .not. allocated (NoahmpIO%TD_SPAC)      ) allocate ( NoahmpIO%TD_SPAC       (1:capacity) )
    !endif

    ! INOUT state and flux variables
    if ( .not. allocated (NoahmpIO%SNOW)     ) allocate ( NoahmpIO%SNOW      (1:capacity) )
    if ( .not. allocated (NoahmpIO%SNOWH)    ) allocate ( NoahmpIO%SNOWH     (1:capacity) )
    if ( .not. allocated (NoahmpIO%SMOISEQ)  ) allocate ( NoahmpIO%SMOISEQ   (1:NSOIL,1:capacity) )
    if ( .not. allocated (NoahmpIO%SMOIS)    ) allocate ( NoahmpIO%SMOIS     (1:NSOIL,1:capacity) )
    if ( .not. allocated (NoahmpIO%SH2O)     ) allocate ( NoahmpIO%SH2O      (1:NSOIL,1:capacity) )
	if ( .not. allocated (NoahmpIO%CANICEXY)  ) allocate ( NoahmpIO%CANICEXY  (1:capacity) )
    if ( .not. allocated (NoahmpIO%CANLIQXY)  ) allocate ( NoahmpIO%CANLIQXY  (1:capacity) )
    if ( .not. allocated (NoahmpIO%FWETXY)    ) allocate ( NoahmpIO%FWETXY    (1:capacity) )
    if ( .not. allocated (NoahmpIO%SNEQVOXY)  ) allocate ( NoahmpIO%SNEQVOXY  (1:capacity) )
    if ( .not. allocated (NoahmpIO%WSLAKEXY)  ) allocate ( NoahmpIO%WSLAKEXY  (1:capacity) )
    if ( .not. allocated (NoahmpIO%ZWTXY)     ) allocate ( NoahmpIO%ZWTXY     (1:capacity) )
    if ( .not. allocated (NoahmpIO%WAXY)      ) allocate ( NoahmpIO%WAXY      (1:capacity) )
    if ( .not. allocated (NoahmpIO%WTXY)      ) allocate ( NoahmpIO%WTXY      (1:capacity) )
    if ( .not. allocated (NoahmpIO%SMCWTDXY)  ) allocate ( NoahmpIO%SMCWTDXY  (1:capacity) )
    if ( .not. allocated (NoahmpIO%SNICEXY)   ) allocate ( NoahmpIO%SNICEXY   (-NSNOW+1:0,    1:capacity) )
    if ( .not. allocated (NoahmpIO%SNLIQXY)   ) allocate ( NoahmpIO%SNLIQXY   (-NSNOW+1:0,    1:capacity) )

! irrigation
    if ( .not. allocated (NoahmpIO%IRFRACT) ) allocate ( NoahmpIO%IRFRACT (1:capacity) )
    if ( .not. allocated (NoahmpIO%SIFRACT) ) allocate ( NoahmpIO%SIFRACT (1:capacity) )
    if ( .not. allocated (NoahmpIO%MIFRACT) ) allocate ( NoahmpIO%MIFRACT (1:capacity) )
    if ( .not. allocated (NoahmpIO%FIFRACT) ) allocate ( NoahmpIO%FIFRACT (1:capacity) )
    if ( .not. allocated (NoahmpIO%IRNUMSI) ) allocate ( NoahmpIO%IRNUMSI (1:capacity) )
    if ( .not. allocated (NoahmpIO%IRNUMMI) ) allocate ( NoahmpIO%IRNUMMI (1:capacity) )
    if ( .not. allocated (NoahmpIO%IRNUMFI) ) allocate ( NoahmpIO%IRNUMFI (1:capacity) )
    if ( .not. allocated (NoahmpIO%IRWATSI) ) allocate ( NoahmpIO%IRWATSI (1:capacity) )
    if ( .not. allocated (NoahmpIO%IRWATMI) ) allocate ( NoahmpIO%IRWATMI (1:capacity) )
    if ( .not. allocated (NoahmpIO%IRWATFI) ) allocate ( NoahmpIO%IRWATFI (1:capacity) )

    ! tile drainage
    if ( .not. allocated (NoahmpIO%TD_FRACTION)) allocate ( NoahmpIO%TD_FRACTION (1:capacity) )

    ! additional output
    if ( .not. allocated (NoahmpIO%ACC_DWATERXY)) allocate ( NoahmpIO%ACC_DWATERXY (1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_PRCPXY)  ) allocate ( NoahmpIO%ACC_PRCPXY   (1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_ECANXY)  ) allocate ( NoahmpIO%ACC_ECANXY   (1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_ETRANXY) ) allocate ( NoahmpIO%ACC_ETRANXY  (1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_EDIRXY)  ) allocate ( NoahmpIO%ACC_EDIRXY   (1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_QINSURXY)) allocate ( NoahmpIO%ACC_QINSURXY (1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_QSEVAXY) ) allocate ( NoahmpIO%ACC_QSEVAXY  (1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_ETRANIXY)) allocate ( NoahmpIO%ACC_ETRANIXY (1:NSOIL,1:capacity) )
    if ( .not. allocated (NoahmpIO%ACC_GLAFLWXY)) allocate ( NoahmpIO%ACC_GLAFLWXY (1:capacity) )

    ! Needed for SNICAR SNOW ALBEDO (IOPT_ALB = 3)
    if ( NoahmpIO%IOPT_ALB == 3 ) then
       allocate (NoahmpIO%snowage_tau(NoahmpIO%idx_rhos_max,NoahmpIO%idx_Tgrd_max,NoahmpIO%idx_T_max))
       allocate (NoahmpIO%snowage_kappa(NoahmpIO%idx_rhos_max,NoahmpIO%idx_Tgrd_max,NoahmpIO%idx_T_max))
       allocate (NoahmpIO%snowage_drdt0(NoahmpIO%idx_rhos_max,NoahmpIO%idx_Tgrd_max,NoahmpIO%idx_T_max))
       allocate (NoahmpIO%SNRDSXY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%SNFRXY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%BCPHIXY(-NSNOW+1:0,1:capacity), NoahmpIO%BCPHOXY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%OCPHIXY(-NSNOW+1:0,1:capacity), NoahmpIO%OCPHOXY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%DUST1XY(-NSNOW+1:0,1:capacity), NoahmpIO%DUST2XY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%DUST3XY(-NSNOW+1:0,1:capacity), NoahmpIO%DUST4XY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%DUST5XY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%MassConcBCPHIXY(-NSNOW+1:0,1:capacity), NoahmpIO%MassConcBCPHOXY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%MassConcOCPHIXY(-NSNOW+1:0,1:capacity), NoahmpIO%MassConcOCPHOXY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%MassConcDUST1XY(-NSNOW+1:0,1:capacity), NoahmpIO%MassConcDUST2XY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%MassConcDUST3XY(-NSNOW+1:0,1:capacity), NoahmpIO%MassConcDUST4XY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%MassConcDUST5XY(-NSNOW+1:0,1:capacity))
       allocate (NoahmpIO%DepBChydrophoXY(1:capacity), NoahmpIO%DepBChydrophiXY(1:capacity))
       allocate (NoahmpIO%DepOChydrophoXY(1:capacity), NoahmpIO%DepOChydrophiXY(1:capacity))
       allocate (NoahmpIO%DepDust1XY(1:capacity), NoahmpIO%DepDust2XY(1:capacity), &
                 NoahmpIO%DepDust3XY(1:capacity), NoahmpIO%DepDust4XY(1:capacity), &
                 NoahmpIO%DepDust5XY(1:capacity))
       NoahmpIO%SNRDSXY = 0.0_kind_noahmp
       NoahmpIO%SNFRXY = 0.0_kind_noahmp
       NoahmpIO%BCPHIXY = 0.0_kind_noahmp
       NoahmpIO%BCPHOXY = 0.0_kind_noahmp
       NoahmpIO%OCPHIXY = 0.0_kind_noahmp
       NoahmpIO%OCPHOXY = 0.0_kind_noahmp
       NoahmpIO%DUST1XY = 0.0_kind_noahmp
       NoahmpIO%DUST2XY = 0.0_kind_noahmp
       NoahmpIO%DUST3XY = 0.0_kind_noahmp
       NoahmpIO%DUST4XY = 0.0_kind_noahmp
       NoahmpIO%DUST5XY = 0.0_kind_noahmp
       NoahmpIO%MassConcBCPHIXY = 0.0_kind_noahmp
       NoahmpIO%MassConcBCPHOXY = 0.0_kind_noahmp
       NoahmpIO%MassConcOCPHIXY = 0.0_kind_noahmp
       NoahmpIO%MassConcOCPHOXY = 0.0_kind_noahmp
       NoahmpIO%MassConcDUST1XY = 0.0_kind_noahmp
       NoahmpIO%MassConcDUST2XY = 0.0_kind_noahmp
       NoahmpIO%MassConcDUST3XY = 0.0_kind_noahmp
       NoahmpIO%MassConcDUST4XY = 0.0_kind_noahmp
       NoahmpIO%MassConcDUST5XY = 0.0_kind_noahmp
       NoahmpIO%DepBChydrophoXY = 0.0_kind_noahmp
       NoahmpIO%DepBChydrophiXY = 0.0_kind_noahmp
       NoahmpIO%DepOChydrophoXY = 0.0_kind_noahmp
       NoahmpIO%DepOChydrophiXY = 0.0_kind_noahmp
       NoahmpIO%DepDust1XY = 0.0_kind_noahmp
       NoahmpIO%DepDust2XY = 0.0_kind_noahmp
       NoahmpIO%DepDust3XY = 0.0_kind_noahmp
       NoahmpIO%DepDust4XY = 0.0_kind_noahmp
       NoahmpIO%DepDust5XY = 0.0_kind_noahmp
    endif

! Needed for Zhang et al. 2022 wetland model (OPT_WETLAND=1 or 2)
    !if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       if ( .not. allocated (NoahmpIO%FSATXY) ) allocate ( NoahmpIO%FSATXY     (1:capacity) )
       if ( .not. allocated (NoahmpIO%WSURFXY)) allocate ( NoahmpIO%WSURFXY    (1:capacity) )
    !endif
    !if ( NoahmpIO%IOPT_WETLAND == 2 ) then
       if ( .not. allocated (NoahmpIO%FSATMX) ) allocate ( NoahmpIO%FSATMX     (1:capacity) )
       if ( .not. allocated (NoahmpIO%WCAP)   ) allocate ( NoahmpIO%WCAP       (1:capacity) )
    !endif


! Default value initializations
    NoahmpIO%SMOIS             = undefined_real
    NoahmpIO%SH2O              = undefined_real
    NoahmpIO%SNOW              = undefined_real
    NoahmpIO%SNOWH             = undefined_real
    NoahmpIO%SMOISEQ           = undefined_real
    NoahmpIO%CANICEXY          = undefined_real
    NoahmpIO%CANLIQXY          = undefined_real
    NoahmpIO%FWETXY            = undefined_real
    NoahmpIO%SNEQVOXY          = undefined_real
    NoahmpIO%WSLAKEXY          = undefined_real
    NoahmpIO%ZWTXY             = undefined_real
    NoahmpIO%WAXY              = undefined_real
    NoahmpIO%WTXY              = undefined_real
    NoahmpIO%SNICEXY           = undefined_real
    NoahmpIO%SNLIQXY           = undefined_real
    NoahmpIO%SMCWTDXY          = undefined_real

    NoahmpIO%ACC_DWATERXY      = 0.0
    NoahmpIO%ACC_PRCPXY        = 0.0
    NoahmpIO%ACC_ECANXY        = 0.0
    NoahmpIO%ACC_ETRANXY       = 0.0
    NoahmpIO%ACC_EDIRXY        = 0.0
    NoahmpIO%ACC_QINSURXY      = 0.0
    NoahmpIO%ACC_QSEVAXY       = 0.0
    NoahmpIO%ACC_ETRANIXY      = 0.0
    NoahmpIO%ACC_GLAFLWXY      = 0.0

    ! tile drainage
    NoahmpIO%TD_FRACTION       = undefined_real

    ! irrigation
    NoahmpIO%IRFRACT           = 0.0
    NoahmpIO%SIFRACT           = 0.0
    NoahmpIO%MIFRACT           = 0.0
    NoahmpIO%FIFRACT           = 0.0
    NoahmpIO%IRNUMSI           = 0
    NoahmpIO%IRNUMMI           = 0
    NoahmpIO%IRNUMFI           = 0
    NoahmpIO%IRWATSI           = 0.0
    NoahmpIO%IRWATMI           = 0.0
    NoahmpIO%IRWATFI           = 0.0

    ! wetland model (Zhang et al. 2022)
    !if ( NoahmpIO%IOPT_WETLAND > 0 ) then
       NoahmpIO%FSATXY       = undefined_real
       NoahmpIO%WSURFXY      = undefined_real
    !endif
    !if ( NoahmpIO%IOPT_WETLAND == 2 ) then
       NoahmpIO%FSATMX       = undefined_real
       NoahmpIO%WCAP         = undefined_real
    !endif

    ! spatial varying soil texture
    !if ( NoahmpIO%IOPT_SOIL > 1 ) then
       NoahmpIO%SOILCOMP     = undefined_real
    !endif


	!water out transfers
	if ( .not. allocated (NoahmpIO%QFX)      ) allocate ( NoahmpIO%QFX       (1:capacity) ) ! latent heat flux [kg s-1 m-2]
	if ( .not. allocated (NoahmpIO%SMSTAV)   ) allocate ( NoahmpIO%SMSTAV    (1:capacity) ) ! soil moisture avail. [not used]
	if ( .not. allocated (NoahmpIO%SMSTOT)   ) allocate ( NoahmpIO%SMSTOT    (1:capacity) ) ! total soil water [mm][not used]
	if ( .not. allocated (NoahmpIO%SFCRUNOFF)) allocate ( NoahmpIO%SFCRUNOFF (1:capacity) ) ! accumulated surface runoff [m]
	if ( .not. allocated (NoahmpIO%UDRUNOFF) ) allocate ( NoahmpIO%UDRUNOFF  (1:capacity) ) ! accumulated sub-surface runoff [m]
	if ( .not. allocated (NoahmpIO%SNOWC)    ) allocate ( NoahmpIO%SNOWC     (1:capacity) ) ! snow cover fraction []
	if ( .not. allocated (NoahmpIO%CANWAT)   ) allocate ( NoahmpIO%CANWAT    (1:capacity) ) ! total canopy water + ice [mm]
	if ( .not. allocated (NoahmpIO%ACSNOM)   ) allocate ( NoahmpIO%ACSNOM    (1:capacity) ) ! accumulated snow melt leaving pack
	if ( .not. allocated (NoahmpIO%ACSNOW)   ) allocate ( NoahmpIO%ACSNOW    (1:capacity) ) ! accumulated snow on grid

	! INOUT (with no Noah LSM equivalent) (as defined in WRF)
	if ( .not. allocated (NoahmpIO%QSNOWXY)   ) allocate ( NoahmpIO%QSNOWXY    (1:capacity) ) ! snowfall on the ground [mm/s]
	if ( .not. allocated (NoahmpIO%QRAINXY)   ) allocate ( NoahmpIO%QRAINXY    (1:capacity) ) ! rainfall on the ground [mm/s]
	if ( .not. allocated (NoahmpIO%DEEPRECHXY)) allocate ( NoahmpIO%DEEPRECHXY (1:capacity) ) ! recharge to the water table when deep (m)
	if ( .not. allocated (NoahmpIO%RECHXY)    ) allocate ( NoahmpIO%RECHXY     (1:capacity) ) ! recharge to the water table (diagnostic) (mm)

	! irrigation
	if ( .not. allocated (NoahmpIO%IRELOSS) ) allocate ( NoahmpIO%IRELOSS (1:capacity) ) ! loss of irrigation water to evaporation,sprinkler [mm]
	if ( .not. allocated (NoahmpIO%IRSIVOL) ) allocate ( NoahmpIO%IRSIVOL (1:capacity) ) ! amount of irrigation by sprinkler (mm)
	if ( .not. allocated (NoahmpIO%IRMIVOL) ) allocate ( NoahmpIO%IRMIVOL (1:capacity) ) ! amount of irrigation by micro (mm)
	if ( .not. allocated (NoahmpIO%IRFIVOL) ) allocate ( NoahmpIO%IRFIVOL (1:capacity) ) ! amount of irrigation by micro (mm)

	! OUT (with no Noah LSM equivalent) (as defined in WRF)
	if ( .not. allocated (NoahmpIO%RUNSFXY) ) allocate ( NoahmpIO%RUNSFXY (1:capacity) ) ! surface runoff [mm per soil timestep]
	if ( .not. allocated (NoahmpIO%RUNSBXY) ) allocate ( NoahmpIO%RUNSBXY (1:capacity) ) ! subsurface runoff [mm per soil timestep]
	if ( .not. allocated (NoahmpIO%ECANXY)  ) allocate ( NoahmpIO%ECANXY  (1:capacity) ) ! evaporation of intercepted water (mm/s)
	if ( .not. allocated (NoahmpIO%EDIRXY)  ) allocate ( NoahmpIO%EDIRXY  (1:capacity) ) ! soil surface evaporation rate (mm/s]
	if ( .not. allocated (NoahmpIO%ETRANXY) ) allocate ( NoahmpIO%ETRANXY (1:capacity) ) ! transpiration rate (mm/s)
	if ( .not. allocated (NoahmpIO%QTDRAIN) ) allocate ( NoahmpIO%QTDRAIN (1:capacity) ) ! tile drainage (mm)

	! additional output variables
	if ( .not. allocated (NoahmpIO%QINTSXY)   ) allocate ( NoahmpIO%QINTSXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%QINTRXY)   ) allocate ( NoahmpIO%QINTRXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%QDRIPSXY)  ) allocate ( NoahmpIO%QDRIPSXY  (1:capacity) )
	if ( .not. allocated (NoahmpIO%QDRIPRXY)  ) allocate ( NoahmpIO%QDRIPRXY  (1:capacity) )
	if ( .not. allocated (NoahmpIO%QTHROSXY)  ) allocate ( NoahmpIO%QTHROSXY  (1:capacity) )
	if ( .not. allocated (NoahmpIO%QTHRORXY)  ) allocate ( NoahmpIO%QTHRORXY  (1:capacity) )
	if ( .not. allocated (NoahmpIO%QSNSUBXY)  ) allocate ( NoahmpIO%QSNSUBXY  (1:capacity) )
	if ( .not. allocated (NoahmpIO%QSNFROXY)  ) allocate ( NoahmpIO%QSNFROXY  (1:capacity) )
	if ( .not. allocated (NoahmpIO%QSUBCXY)   ) allocate ( NoahmpIO%QSUBCXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%QFROCXY)   ) allocate ( NoahmpIO%QFROCXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%QEVACXY)   ) allocate ( NoahmpIO%QEVACXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%QDEWCXY)   ) allocate ( NoahmpIO%QDEWCXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%QFRZCXY)   ) allocate ( NoahmpIO%QFRZCXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%QMELTCXY)  ) allocate ( NoahmpIO%QMELTCXY  (1:capacity) )
	if ( .not. allocated (NoahmpIO%QSNBOTXY)  ) allocate ( NoahmpIO%QSNBOTXY  (1:capacity) )
	if ( .not. allocated (NoahmpIO%QMELTXY)   ) allocate ( NoahmpIO%QMELTXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%PONDINGXY) ) allocate ( NoahmpIO%PONDINGXY (1:capacity) )
	if ( .not. allocated (NoahmpIO%FPICEXY)   ) allocate ( NoahmpIO%FPICEXY   (1:capacity) )
	if ( .not. allocated (NoahmpIO%RAINLSM)   ) allocate ( NoahmpIO%RAINLSM   (1:capacity) )
	if ( .not. allocated (NoahmpIO%SNOWLSM)   ) allocate ( NoahmpIO%SNOWLSM   (1:capacity) )


	NoahmpIO%QFX             = undefined_real
	NoahmpIO%SMSTAV          = undefined_real
	NoahmpIO%SMSTOT          = undefined_real
	NoahmpIO%CANWAT          = undefined_real
	NoahmpIO%SNOWC           = undefined_real

	NoahmpIO%QSNOWXY         = undefined_real
	NoahmpIO%QRAINXY         = undefined_real

	NoahmpIO%RUNSFXY         = undefined_real
	NoahmpIO%RUNSBXY         = undefined_real
	NoahmpIO%ECANXY          = undefined_real
	NoahmpIO%EDIRXY          = undefined_real
	NoahmpIO%ETRANXY         = undefined_real

	NoahmpIO%DEEPRECHXY      = 0.0
	NoahmpIO%RECHXY          = 0.0
	NoahmpIO%ACSNOM          = 0.0
	NoahmpIO%ACSNOW          = 0.0
	NoahmpIO%SFCRUNOFF       = 0.0
	NoahmpIO%UDRUNOFF        = 0.0

	! additional output
	NoahmpIO%QINTSXY         = undefined_real
	NoahmpIO%QINTRXY         = undefined_real
	NoahmpIO%QDRIPSXY        = undefined_real
	NoahmpIO%QDRIPRXY        = undefined_real
	NoahmpIO%QTHROSXY        = undefined_real
	NoahmpIO%QTHRORXY        = undefined_real
	NoahmpIO%QSNSUBXY        = undefined_real
	NoahmpIO%QSNFROXY        = undefined_real
	NoahmpIO%QSUBCXY         = undefined_real
	NoahmpIO%QFROCXY         = undefined_real
	NoahmpIO%QEVACXY         = undefined_real
	NoahmpIO%QDEWCXY         = undefined_real
	NoahmpIO%QFRZCXY         = undefined_real
	NoahmpIO%QMELTCXY        = undefined_real
	NoahmpIO%QSNBOTXY        = undefined_real
	NoahmpIO%QMELTXY         = undefined_real
	NoahmpIO%FPICEXY         = undefined_real
	NoahmpIO%RAINLSM         = undefined_real
	NoahmpIO%SNOWLSM         = undefined_real
	NoahmpIO%PONDINGXY       = 0.0

	! tile drainage
	NoahmpIO%QTDRAIN         = 0.0

	! irrigation
	NoahmpIO%IRELOSS         = 0.0
	NoahmpIO%IRSIVOL         = 0.0
	NoahmpIO%IRMIVOL         = 0.0
	NoahmpIO%IRFIVOL         = 0.0
	
	!biochem in transfers
	if ( .not. allocated(NoahmpIO%PGSXY) ) allocate(NoahmpIO%PGSXY(1:capacity))
	if ( .not. allocated(NoahmpIO%PLANTING) ) allocate(NoahmpIO%PLANTING(1:capacity))
	if ( .not. allocated(NoahmpIO%HARVEST) ) allocate(NoahmpIO%HARVEST(1:capacity))
	if ( .not. allocated(NoahmpIO%SEASON_GDD) ) allocate(NoahmpIO%SEASON_GDD(1:capacity))	   
	if ( .not. allocated(NoahmpIO%LFMASSXY) ) allocate(NoahmpIO%LFMASSXY(1:capacity))
	if ( .not. allocated(NoahmpIO%RTMASSXY) ) allocate(NoahmpIO%RTMASSXY(1:capacity))
	if ( .not. allocated(NoahmpIO%STMASSXY) ) allocate(NoahmpIO%STMASSXY(1:capacity))
	if ( .not. allocated(NoahmpIO%WOODXY) ) allocate(NoahmpIO%WOODXY(1:capacity))
	if ( .not. allocated(NoahmpIO%GRAINXY) ) allocate(NoahmpIO%GRAINXY(1:capacity))
	if ( .not. allocated(NoahmpIO%GDDXY) ) allocate(NoahmpIO%GDDXY(1:capacity))
	if ( .not. allocated(NoahmpIO%STBLCPXY) ) allocate(NoahmpIO%STBLCPXY(1:capacity))
	if ( .not. allocated(NoahmpIO%FASTCPXY) ) allocate(NoahmpIO%FASTCPXY(1:capacity))
   
	NoahmpIO%LFMASSXY = undefined_real
	NoahmpIO%RTMASSXY = undefined_real
	NoahmpIO%STMASSXY = undefined_real
	NoahmpIO%WOODXY   = undefined_real

	NoahmpIO%STBLCPXY = undefined_real
	NoahmpIO%FASTCPXY = undefined_real
	
	NoahmpIO%GRAINXY = undefined_real
	NoahmpIO%GDDXY   = undefined_real

	NoahmpIO%PGSXY      = undefined_int
	NoahmpIO%PLANTING   = undefined_real
	NoahmpIO%HARVEST    = undefined_real
	NoahmpIO%SEASON_GDD = undefined_real
	
	!Biochem out transfers
	if ( .not. allocated(NoahmpIO%NEEXY) ) allocate(NoahmpIO%NEEXY(1:capacity))
	if ( .not. allocated(NoahmpIO%GPPXY) ) allocate(NoahmpIO%GPPXY(1:capacity))
	if ( .not. allocated(NoahmpIO%NPPXY) ) allocate(NoahmpIO%NPPXY(1:capacity))
	if ( .not. allocated(NoahmpIO%PSNXY) ) allocate(NoahmpIO%PSNXY(1:capacity))
	
	NoahmpIO%NEEXY(1:capacity) = undefined_real
	NoahmpIO%GPPXY(1:capacity) = undefined_real
	NoahmpIO%NPPXY(1:capacity) = undefined_real
	NoahmpIO%PSNXY(1:capacity) = undefined_real
	
	end associate
	   
  end subroutine NoahmpIOVarInitDefault

end module NoahmpIOVarInitMod
