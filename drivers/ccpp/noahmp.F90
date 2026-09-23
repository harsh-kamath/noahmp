! Harsh Kamath (NOAA GSL/CIRES), September 2026


module noahmp

  use Machine,               only : kind_phys
  use NoahmpIOVarType, only : NoahmpIO_type
  use NoahmpIOVarInitMod, only : NoahmpIOVarInitDefault
  use NoahmpReadTableMod, only : NoahmpReadTable
  use NoahmpDriverMainMod, only : NoahmpDriverMain
  use SnowInputSnicarMod, only : SnowInputSnicar
  use StabilityFunctionsMod, only : InitializeStabilityFunctions
  use ConstantDefineMod, only : ConstGasDryAir, ConstHeatCapacAir, ConstLatHeatEvap, &
                                ConstLatHeatSublim, ConstLatHeatFusion

  implicit none
  private



  integer, parameter :: noahmp_shortwave_radiation_bands = 2
  character(len=*), parameter :: legacy_ccpp_land_use_data_name = &
								"MODIFIED_IGBP_MODIS_NOAH"
								
  type :: NoahmpCcppContext_type
     type(NoahmpIO_type), allocatable :: io
  end type NoahmpCcppContext_type

  ! Each CCPP instance owns an independent driver context and workspace.
  type(NoahmpCcppContext_type), allocatable :: NoahmpContexts(:)

  public :: noahmp_init
  public :: noahmp_run
  public :: noahmp_final

contains

  !> \section arg_table_noahmp_init Argument Table
  !! \htmlinclude noahmp_init.html
  !!
  subroutine noahmp_init(instance_number, number_of_instances, horizontal_dimension, &
    vertical_dimension_of_soil, vertical_dimension_of_surface_snow, &
    dynamic_vegetation_option, rain_snow_partition_option, soil_water_transpiration_option, &
    ground_resistance_evaporation_option, surface_drag_option, surface_stability_function_option, &
    surface_thermal_roughness_option, stomata_resistance_option, snow_albedo_option, &
    canopy_radiation_transfer_option, snow_soil_temperature_time_option, snow_thermal_conductivity_option, &
    soil_temperature_bottom_option, soil_supercooled_water_option, frozen_soil_permeability_option, &
    dynamic_vic_infiltration_option, tile_drainage_option, irrigation_option, &
    irrigation_method_option, crop_model_option, soil_property_option, &
    pedotransfer_option, surface_runoff_option, subsurface_runoff_option, &
    glacier_treatment_option, snow_compaction_option, wetland_model_option, &
    snow_cover_fraction_option, snicar_snow_shape_option, snicar_rt_solver_option, &
    snicar_band_number_option, snicar_solar_spectrum_option, snicar_snow_optics_option, &
    snicar_dust_optics_option, snicar_snow_bc_internal_mixing, snicar_snow_dust_internal_mixing, &
    snicar_use_aerosol, snicar_use_organic_carbon, snicar_read_aerosol_table, &
    urban_physics_option, errmsg, errflg)
	
    integer, intent(in) :: instance_number
    integer, intent(in) :: number_of_instances
    integer, intent(in) :: horizontal_dimension
    integer, intent(in) :: vertical_dimension_of_soil
    integer, intent(in) :: vertical_dimension_of_surface_snow
    integer, parameter  :: number_of_shortwave_radiation_bands = &
							noahmp_shortwave_radiation_bands
    integer, intent(in) :: dynamic_vegetation_option
    integer, intent(in) :: rain_snow_partition_option
    integer, intent(in) :: soil_water_transpiration_option
    integer, intent(in) :: ground_resistance_evaporation_option
    integer, intent(in) :: surface_drag_option
    integer, intent(in) :: surface_stability_function_option
    integer, intent(in) :: surface_thermal_roughness_option
    integer, intent(in) :: stomata_resistance_option
    integer, intent(in) :: snow_albedo_option
    integer, intent(in) :: canopy_radiation_transfer_option
    integer, intent(in) :: snow_soil_temperature_time_option
    integer, intent(in) :: snow_thermal_conductivity_option
    integer, intent(in) :: soil_temperature_bottom_option
    integer, intent(in) :: soil_supercooled_water_option
    integer, intent(in) :: frozen_soil_permeability_option
    integer, intent(in) :: dynamic_vic_infiltration_option
    integer, intent(in) :: tile_drainage_option
    integer, intent(in) :: irrigation_option
    integer, intent(in) :: irrigation_method_option
    integer, intent(in) :: crop_model_option
    integer, intent(in) :: soil_property_option
    integer, intent(in) :: pedotransfer_option
    integer, intent(in) :: surface_runoff_option
    integer, intent(in) :: subsurface_runoff_option
    integer, intent(in) :: glacier_treatment_option
    integer, intent(in) :: snow_compaction_option
    integer, intent(in) :: wetland_model_option
    integer, intent(in) :: snow_cover_fraction_option
    integer, intent(in) :: snicar_snow_shape_option
    integer, intent(in) :: snicar_rt_solver_option
    integer, intent(in) :: snicar_band_number_option
    integer, intent(in) :: snicar_solar_spectrum_option
    integer, intent(in) :: snicar_snow_optics_option
    integer, intent(in) :: snicar_dust_optics_option
    logical, intent(in) :: snicar_snow_bc_internal_mixing
    logical, intent(in) :: snicar_snow_dust_internal_mixing
    logical, intent(in) :: snicar_use_aerosol
    logical, intent(in) :: snicar_use_organic_carbon
    logical, intent(in) :: snicar_read_aerosol_table
    integer, intent(in) :: urban_physics_option

    character(len=*), intent(out) :: errmsg
    integer, intent(out)          :: errflg

    integer :: allocation_status
    logical :: context_registry_created

    errmsg = ''
    errflg = 0

    if (number_of_instances <= 0) then
       errflg = 1
       errmsg = 'noahmp_init: number_of_instances must be positive'
       return
    end if
    if (instance_number < 1 .or. instance_number > number_of_instances) then
       errflg = 1
       errmsg = 'noahmp_init: instance_number is outside number_of_instances'
       return
    end if

    context_registry_created = .false.
    if (allocated(NoahmpContexts)) then
       if (size(NoahmpContexts) /= number_of_instances) then
          errflg = 1
          errmsg = 'noahmp_init: number_of_instances changed after context allocation'
          return
       end if
       if (allocated(NoahmpContexts(instance_number)%io)) return
    end if

    if (horizontal_dimension <= 0 .or. vertical_dimension_of_soil <= 0 .or. &
        vertical_dimension_of_surface_snow <= 0 .or. &
        number_of_shortwave_radiation_bands <= 0) then
       errflg = 1
       errmsg = 'noahmp_init: horizontal, soil, snow, and radiation dimensions must be positive'
       return
    end if

    if (surface_drag_option < 1 .or. surface_drag_option > 4) then
       errflg = 1
       errmsg = 'noahmp_init: surface_drag_option must be in the range 1 through 4'
       return
    end if
    if (surface_thermal_roughness_option < 0 .or. &
        surface_thermal_roughness_option > 3) then
       errflg = 1
       errmsg = 'noahmp_init: surface_thermal_roughness_option must be in the range 0 through 3'
       return
    end if

    if (snow_albedo_option == 3 .and. snicar_band_number_option /= 1 .and. &
        snicar_band_number_option /= 2) then
       errflg = 1
       errmsg = 'noahmp_init: snicar_band_number_option must be 1 or 2'
       return
    end if

    call InitializeStabilityFunctions(surface_stability_function_option, errmsg, errflg)
    if (errflg /= 0) return

    if (.not. allocated(NoahmpContexts)) then
       allocate(NoahmpContexts(number_of_instances), stat=allocation_status)
       if (allocation_status /= 0) then
          errflg = 1
          write(errmsg, '(a,i0)') &
               'noahmp_init: unable to allocate the context registry; stat=', &
               allocation_status
          return
       end if
       context_registry_created = .true.
    end if

    allocate(NoahmpContexts(instance_number)%io, stat=allocation_status)
    if (allocation_status /= 0) then
       errflg = 1
       write(errmsg, '(a,i0,a,i0)') &
            'noahmp_init: unable to allocate context ', instance_number, &
            '; stat=', allocation_status
       if (context_registry_created) deallocate(NoahmpContexts)
       return
    end if
    associate (NoahmpIO => NoahmpContexts(instance_number)%io)
	
    NoahmpIO%KMS = 1
    NoahmpIO%KME = 2
    NoahmpIO%NSOIL = vertical_dimension_of_soil
    NoahmpIO%NSNOW = vertical_dimension_of_surface_snow
    NoahmpIO%NUMRAD = number_of_shortwave_radiation_bands

    ! The selected instance context owns configuration after initialization.
    NoahmpIO%IOPT_DVEG = dynamic_vegetation_option
    NoahmpIO%IOPT_SNF = rain_snow_partition_option
    NoahmpIO%IOPT_BTR = soil_water_transpiration_option
    NoahmpIO%IOPT_RSF = ground_resistance_evaporation_option
    NoahmpIO%IOPT_SFC = surface_drag_option
    NoahmpIO%PSI_OPT = surface_stability_function_option
    NoahmpIO%IZ0TLND = surface_thermal_roughness_option
    NoahmpIO%IOPT_CRS = stomata_resistance_option
    NoahmpIO%IOPT_ALB = snow_albedo_option
    NoahmpIO%IOPT_RAD = canopy_radiation_transfer_option
    NoahmpIO%IOPT_STC = snow_soil_temperature_time_option
    NoahmpIO%IOPT_TKSNO = snow_thermal_conductivity_option
    NoahmpIO%IOPT_TBOT = soil_temperature_bottom_option
    NoahmpIO%IOPT_FRZ = soil_supercooled_water_option
    NoahmpIO%IOPT_INF = frozen_soil_permeability_option
    NoahmpIO%IOPT_INFDV = dynamic_vic_infiltration_option
    NoahmpIO%IOPT_TDRN = tile_drainage_option
    NoahmpIO%IOPT_IRR = irrigation_option
    NoahmpIO%IOPT_IRRM = irrigation_method_option
    NoahmpIO%IOPT_CROP = crop_model_option
    NoahmpIO%IOPT_SOIL = soil_property_option
    NoahmpIO%IOPT_PEDO = pedotransfer_option
    NoahmpIO%IOPT_RUNSRF = surface_runoff_option
    NoahmpIO%IOPT_RUNSUB = subsurface_runoff_option
    NoahmpIO%IOPT_GLA = glacier_treatment_option
    NoahmpIO%IOPT_COMPACT = snow_compaction_option
    NoahmpIO%IOPT_WETLAND = wetland_model_option
    NoahmpIO%IOPT_SCF = snow_cover_fraction_option
    NoahmpIO%SNICAR_SNOWSHAPE_OPT = snicar_snow_shape_option
    NoahmpIO%SNICAR_RTSOLVER_OPT = snicar_rt_solver_option
    NoahmpIO%SNICAR_BANDNUMBER_OPT = snicar_band_number_option
    NoahmpIO%SNICAR_SOLARSPEC_OPT = snicar_solar_spectrum_option
    NoahmpIO%SNICAR_SNOWOPTICS_OPT = snicar_snow_optics_option
    NoahmpIO%SNICAR_DUSTOPTICS_OPT = snicar_dust_optics_option
    NoahmpIO%SNICAR_SNOWBC_INTMIX = snicar_snow_bc_internal_mixing
    NoahmpIO%SNICAR_SNOWDUST_INTMIX = snicar_snow_dust_internal_mixing
    NoahmpIO%SNICAR_USE_AEROSOL = snicar_use_aerosol
    NoahmpIO%SNICAR_USE_OC = snicar_use_organic_carbon
    NoahmpIO%SNICAR_AEROSOL_READTABLE = snicar_read_aerosol_table
    NoahmpIO%SF_URBAN_PHYSICS = urban_physics_option
	
    if (snow_albedo_option == 3) then
       select case (snicar_band_number_option)
       case (1)
          NoahmpIO%snicar_numrad_snw = 5
          NoahmpIO%snicar_optic_flnm = 'snicar_optics_5bnd_c013122.nc'
       case (2)
          NoahmpIO%snicar_numrad_snw = 480
          NoahmpIO%snicar_optic_flnm = 'snicar_optics_480bnd_c012422.nc'
       end select
       NoahmpIO%snicar_age_flnm = 'snicar_drdt_bst_fit_60_c070416.nc'
    endif
	
    call NoahmpIOVarInitDefault(NoahmpIO, horizontal_dimension)

    NoahmpIO%LLANDUSE = legacy_ccpp_land_use_data_name
	
    call NoahmpReadTable(NoahmpIO)
	
    if (snow_albedo_option == 3) then
       call SnowInputSnicar(NoahmpIO)
       NoahmpIO%SNRDSXY = NoahmpIO%SnowRadiusMin_TABLE
    endif

    end associate

  end subroutine noahmp_init


  !> \section arg_table_noahmp_run Argument Table
  !! \htmlinclude noahmp_run.html
  !!
  subroutine noahmp_run(instance_number, horizontal_loop_extent, timestep_index, timestep_for_physics, &
    forecast_julian_day, number_of_days_in_current_year, &
    characteristic_grid_lengthscale, vertical_dimension_of_soil, &
    vertical_dimension_of_surface_snow, &
    lower_bound_of_vertical_dimension_of_surface_snow, &
    depth_of_soil_layer_interfaces, vegetation_category, soil_category, &
    soil_color_classification, surface_slope_classification, &
    active_snow_layer_lower_index, depth_of_snow_soil_layer_interfaces, &
    atmospheric_forcing_height, atmosphere_boundary_layer_thickness, &
    flag_for_reference_pressure_theta, leaf_area_index, stem_area_index, &
    specific_humidity_at_surface, ground_temperature, canopy_temperature, &
    dimensionless_snow_age, snow_albedo_on_previous_timestep, &
    canopy_air_vapor_pressure, canopy_air_temperature, &
    surface_exchange_coefficient_for_heat, &
    surface_exchange_coefficient_for_momentum, &
    maximum_vegetation_fraction, vegetation_fraction, soil_temperature, &
    snow_temperature, canopy_liquid_water, canopy_ice_water, &
    canopy_wet_fraction, snow_water_equivalent, &
    snow_water_equivalent_on_previous_timestep, &
    snow_depth, flood_irrigation_fraction, flood_irrigation_water_amount, &
    micro_irrigation_fraction, micro_irrigation_water_amount, &
    sprinkler_irrigation_fraction, sprinkler_irrigation_water_amount, &
    water_table_depth, soil_moisture_below_soil_column, &
    tile_drainage_fraction, aquifer_water_storage, &
    soil_aquifer_water_storage, lake_water_storage, &
    grid_irrigation_fraction, sprinkler_irrigation_event_count, &
    micro_irrigation_event_count, flood_irrigation_event_count, &
    soil_saturated_fraction, wetland_water_storage, &
    snow_layer_ice, snow_layer_liquid_water, soil_liquid_water, &
    soil_moisture, equilibrium_soil_moisture, leaf_mass, &
    root_mass, stem_mass, wood_mass, &
    deep_soil_carbon_mass, shallow_soil_carbon_mass, &
    accumulated_ground_heat_flux, sprinkler_heat_accumulation, &
    accumulated_soil_surface_evaporation, accumulated_soil_surface_inflow, &
    accumulated_surface_water_change, accumulated_precipitation, &
    accumulated_canopy_evaporation, accumulated_transpiration, &
    accumulated_ground_evaporation, accumulated_soil_layer_transpiration, &
    accumulated_glacier_excess_flow, accumulated_surface_runoff, &
    accumulated_subsurface_runoff, accumulated_tile_drainage, &
    accumulated_snowfall, accumulated_snowmelt, &
    accumulated_shallow_groundwater_recharge, &
    accumulated_deep_groundwater_recharge, accumulated_sprinkler_irrigation, &
    accumulated_micro_irrigation, accumulated_flood_irrigation, &
    accumulated_sprinkler_evaporation_loss, direct_soil_albedo, &
    diffuse_soil_albedo, plant_growth_stage, grain_mass, &
    growing_degree_day, &
    wilting_soil_moisture, reference_soil_moisture, &
    surface_radiative_temperature, surface_emissivity, &
    surface_roughness_length, temperature_at_2m_from_noahmp, &
    specific_humidity_at_2m_from_noahmp, &
    surface_drag_wind_speed_for_momentum, &
    surface_drag_mass_flux_for_heat_and_moisture, &
    surface_albedo, canopy_resistance, latent_heat_flux, &
    sensible_heat_flux, ground_heat_flux, soil_upward_latent_heat_flux, &
    transpiration_latent_heat_flux, canopy_upward_latent_heat_flux, &
    precipitation_advected_heat_flux, soil_moisture_content, &
    normalized_soil_wetness, snow_cover_fraction, total_canopy_water, &
    ground_snowfall_rate, surface_runoff, surface_runoff_flux, &
    subsurface_runoff, net_canopy_evaporation, net_ground_evaporation, &
    transpiration, snow_deposition_sublimation_upward_latent_heat_flux, &
    snow_freezing_rain_upward_latent_heat_flux, &
    direct_visible_surface_albedo, direct_nir_surface_albedo, &
    diffuse_visible_surface_albedo, diffuse_nir_surface_albedo, &
    latitude, instantaneous_cosine_of_zenith_angle, &
    surface_downwelling_shortwave_flux, surface_downwelling_longwave_flux, &
    direct_nir_shortwave_flux_at_surface, &
    diffuse_nir_shortwave_flux_at_surface, &
    direct_vis_shortwave_flux_at_surface, &
    diffuse_vis_shortwave_flux_at_surface, &
    snicar_hydrophobic_black_carbon_deposition_flux, &
    snicar_hydrophilic_black_carbon_deposition_flux, &
    snicar_hydrophobic_organic_carbon_deposition_flux, &
    snicar_hydrophilic_organic_carbon_deposition_flux, &
    snicar_dust_species_1_deposition_flux, &
    snicar_dust_species_2_deposition_flux, &
    snicar_dust_species_3_deposition_flux, &
    snicar_dust_species_4_deposition_flux, &
    snicar_dust_species_5_deposition_flux, &
    surface_air_pressure, air_pressure_at_surface_adjacent_layer, &
    deep_soil_temperature, air_temperature_at_lowest_model_layer, &
    specific_humidity_at_lowest_model_layer, &
    eastward_wind_at_lowest_model_layer, &
    northward_wind_at_lowest_model_layer, &
    convective_precipitation_rate_on_previous_timestep, &
    explicit_precipitation_rate_on_previous_timestep, &
    precipitation_amount_on_dynamics_timestep, &
    snowfall_rate_on_previous_timestep, &
    graupel_precipitation_rate_on_previous_timestep, &
    ice_precipitation_rate_on_previous_timestep, precipitation_type, &
    flag_for_iteration, flag_nonzero_land_surface_fraction, &
    apply_urban_irrigation, errmsg, errflg)
	
    integer, intent(in)                    :: instance_number
    integer, intent(in)                    :: horizontal_loop_extent
    integer, intent(in)                    :: timestep_index
    real(kind=kind_phys), intent(in)        :: timestep_for_physics
    real(kind=kind_phys), intent(in)     :: forecast_julian_day
    integer, intent(in)                    :: number_of_days_in_current_year
    real(kind=kind_phys), intent(in)       :: characteristic_grid_lengthscale(:)
    integer, intent(in)                    :: vertical_dimension_of_soil
    integer, intent(in)                    :: vertical_dimension_of_surface_snow
    integer, intent(in)                    :: lower_bound_of_vertical_dimension_of_surface_snow
    integer, parameter                     :: upper_bound_of_vertical_dimension_of_surface_snow = 0
    integer, parameter                     :: number_of_shortwave_radiation_bands = &
         noahmp_shortwave_radiation_bands
    real(kind=kind_phys), intent(in)     :: depth_of_soil_layer_interfaces(:)
    integer, intent(in)                    :: vegetation_category(:)
    integer, intent(in)                    :: soil_category(:)
    integer, intent(in)                    :: soil_color_classification(:)
    integer, intent(in)                    :: surface_slope_classification(:)
    real(kind=kind_phys), intent(inout), optional :: active_snow_layer_lower_index(:)
    real(kind=kind_phys), intent(inout), optional :: depth_of_snow_soil_layer_interfaces( &
                                              :,lower_bound_of_vertical_dimension_of_surface_snow:)
    real(kind=kind_phys), intent(in)       :: atmospheric_forcing_height(:)
    real(kind=kind_phys), intent(in)       :: atmosphere_boundary_layer_thickness(:)
    logical, intent(in)                    :: flag_for_reference_pressure_theta
    real(kind=kind_phys), intent(inout), optional :: leaf_area_index(:)
    real(kind=kind_phys), intent(inout), optional :: stem_area_index(:)
    real(kind=kind_phys), intent(inout)  :: specific_humidity_at_surface(:)
    real(kind=kind_phys), intent(inout), optional :: ground_temperature(:)
    real(kind=kind_phys), intent(inout), optional :: canopy_temperature(:)
    real(kind=kind_phys), intent(inout), optional :: dimensionless_snow_age(:)
    real(kind=kind_phys), intent(inout), optional :: snow_albedo_on_previous_timestep(:)
    real(kind=kind_phys), intent(inout), optional :: canopy_air_vapor_pressure(:)
    real(kind=kind_phys), intent(inout), optional :: canopy_air_temperature(:)
    real(kind=kind_phys), intent(inout), optional :: surface_exchange_coefficient_for_heat(:)
    real(kind=kind_phys), intent(inout), optional :: surface_exchange_coefficient_for_momentum(:)
    real(kind=kind_phys), intent(inout)  :: maximum_vegetation_fraction(:)
    real(kind=kind_phys), intent(inout)  :: vegetation_fraction(:)
    real(kind=kind_phys), intent(inout)  :: soil_temperature(:,:)
    real(kind=kind_phys), intent(inout), optional :: snow_temperature( &
                                              :,lower_bound_of_vertical_dimension_of_surface_snow:)
    real(kind=kind_phys), intent(inout), optional :: canopy_liquid_water(:)
    real(kind=kind_phys), intent(inout), optional :: canopy_ice_water(:)
    real(kind=kind_phys), intent(inout), optional :: canopy_wet_fraction(:)
    real(kind=kind_phys), intent(inout)  :: snow_water_equivalent(:)
    real(kind=kind_phys), intent(inout), optional :: snow_water_equivalent_on_previous_timestep(:)
    real(kind=kind_phys), intent(inout)  :: snow_depth(:)
    real(kind=kind_phys), intent(in), optional    :: flood_irrigation_fraction(:)
    real(kind=kind_phys), intent(inout), optional :: flood_irrigation_water_amount(:)
    real(kind=kind_phys), intent(in), optional    :: micro_irrigation_fraction(:)
    real(kind=kind_phys), intent(inout), optional :: micro_irrigation_water_amount(:)
    real(kind=kind_phys), intent(in), optional    :: sprinkler_irrigation_fraction(:)
    real(kind=kind_phys), intent(inout), optional :: sprinkler_irrigation_water_amount(:)
    real(kind=kind_phys), intent(inout), optional :: water_table_depth(:)
    real(kind=kind_phys), intent(inout), optional :: soil_moisture_below_soil_column(:)
    real(kind=kind_phys), intent(in), optional :: tile_drainage_fraction(:)
    real(kind=kind_phys), intent(inout), optional :: aquifer_water_storage(:)
    real(kind=kind_phys), intent(inout), optional :: soil_aquifer_water_storage(:)
    real(kind=kind_phys), intent(inout), optional :: lake_water_storage(:)
    real(kind=kind_phys), intent(in), optional    :: grid_irrigation_fraction(:)
    integer, intent(inout), optional              :: sprinkler_irrigation_event_count(:)
    integer, intent(inout), optional              :: micro_irrigation_event_count(:)
    integer, intent(inout), optional              :: flood_irrigation_event_count(:)
    real(kind=kind_phys), intent(inout), optional :: soil_saturated_fraction(:)
    real(kind=kind_phys), intent(inout), optional :: wetland_water_storage(:)
    real(kind=kind_phys), intent(inout), optional :: snow_layer_ice( &
                                              :,lower_bound_of_vertical_dimension_of_surface_snow:)
    real(kind=kind_phys), intent(inout), optional :: snow_layer_liquid_water( &
                                              :,lower_bound_of_vertical_dimension_of_surface_snow:)
    real(kind=kind_phys), intent(inout)  :: soil_liquid_water(:,:)
    real(kind=kind_phys), intent(inout)  :: soil_moisture(:,:)
    real(kind=kind_phys), intent(in), optional :: equilibrium_soil_moisture(:,:)
    real(kind=kind_phys), intent(inout), optional :: leaf_mass(:)
    real(kind=kind_phys), intent(inout), optional :: root_mass(:)
    real(kind=kind_phys), intent(inout), optional :: stem_mass(:)
    real(kind=kind_phys), intent(inout), optional :: wood_mass(:)
    real(kind=kind_phys), intent(inout), optional :: deep_soil_carbon_mass(:)
    real(kind=kind_phys), intent(inout), optional :: shallow_soil_carbon_mass(:)
    real(kind=kind_phys), intent(inout) :: accumulated_ground_heat_flux(:)
    real(kind=kind_phys), intent(inout) :: sprinkler_heat_accumulation(:)
    real(kind=kind_phys), intent(inout) :: accumulated_soil_surface_evaporation(:)
    real(kind=kind_phys), intent(inout) :: accumulated_soil_surface_inflow(:)
    real(kind=kind_phys), intent(inout) :: accumulated_surface_water_change(:)
    real(kind=kind_phys), intent(inout) :: accumulated_precipitation(:)
    real(kind=kind_phys), intent(inout) :: accumulated_canopy_evaporation(:)
    real(kind=kind_phys), intent(inout) :: accumulated_transpiration(:)
    real(kind=kind_phys), intent(inout) :: accumulated_ground_evaporation(:)
    real(kind=kind_phys), intent(inout) :: accumulated_soil_layer_transpiration(:,:)
    real(kind=kind_phys), intent(inout) :: accumulated_glacier_excess_flow(:)
    real(kind=kind_phys), intent(inout) :: accumulated_surface_runoff(:)
    real(kind=kind_phys), intent(inout) :: accumulated_subsurface_runoff(:)
    real(kind=kind_phys), intent(inout) :: accumulated_tile_drainage(:)
    real(kind=kind_phys), intent(inout) :: accumulated_snowfall(:)
    real(kind=kind_phys), intent(inout) :: accumulated_snowmelt(:)
    real(kind=kind_phys), intent(inout) :: accumulated_shallow_groundwater_recharge(:)
    real(kind=kind_phys), intent(inout) :: accumulated_deep_groundwater_recharge(:)
    real(kind=kind_phys), intent(inout) :: accumulated_sprinkler_irrigation(:)
    real(kind=kind_phys), intent(inout) :: accumulated_micro_irrigation(:)
    real(kind=kind_phys), intent(inout) :: accumulated_flood_irrigation(:)
    real(kind=kind_phys), intent(inout) :: accumulated_sprinkler_evaporation_loss(:)
    real(kind=kind_phys), intent(inout) :: direct_soil_albedo(:,:)
    real(kind=kind_phys), intent(inout) :: diffuse_soil_albedo(:,:)
    integer, intent(inout) :: plant_growth_stage(:)
    real(kind=kind_phys), intent(inout) :: grain_mass(:)
    real(kind=kind_phys), intent(inout) :: growing_degree_day(:)
    real(kind=kind_phys), intent(inout)  :: wilting_soil_moisture(:)
    real(kind=kind_phys), intent(inout)  :: reference_soil_moisture(:)
    real(kind=kind_phys), intent(inout)    :: surface_radiative_temperature(:)
    real(kind=kind_phys), intent(inout)    :: surface_emissivity(:)
    real(kind=kind_phys), intent(inout)    :: surface_roughness_length(:)
    real(kind=kind_phys), intent(inout), optional :: temperature_at_2m_from_noahmp(:)
    real(kind=kind_phys), intent(inout), optional :: specific_humidity_at_2m_from_noahmp(:)
    real(kind=kind_phys), intent(inout) :: surface_drag_wind_speed_for_momentum(:)
    real(kind=kind_phys), intent(inout) :: surface_drag_mass_flux_for_heat_and_moisture(:)
    real(kind=kind_phys), intent(inout)    :: surface_albedo(:)
    real(kind=kind_phys), intent(inout), optional :: canopy_resistance(:)
    real(kind=kind_phys), intent(inout)    :: latent_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: sensible_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: ground_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: soil_upward_latent_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: transpiration_latent_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: canopy_upward_latent_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: precipitation_advected_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: soil_moisture_content(:)
    real(kind=kind_phys), intent(inout), optional :: normalized_soil_wetness(:)
    real(kind=kind_phys), intent(inout)    :: snow_cover_fraction(:)
    real(kind=kind_phys), intent(inout)    :: total_canopy_water(:)
    real(kind=kind_phys), intent(inout), optional :: ground_snowfall_rate(:)
    real(kind=kind_phys), intent(inout)    :: surface_runoff(:)
    real(kind=kind_phys), intent(inout)    :: surface_runoff_flux(:)
    real(kind=kind_phys), intent(inout)    :: subsurface_runoff(:)
    real(kind=kind_phys), intent(inout)    :: net_canopy_evaporation(:)
    real(kind=kind_phys), intent(inout)    :: net_ground_evaporation(:)
    real(kind=kind_phys), intent(inout)    :: transpiration(:)
    real(kind=kind_phys), intent(inout)    :: snow_deposition_sublimation_upward_latent_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: snow_freezing_rain_upward_latent_heat_flux(:)
    real(kind=kind_phys), intent(inout)    :: direct_visible_surface_albedo(:)
    real(kind=kind_phys), intent(inout)    :: direct_nir_surface_albedo(:)
    real(kind=kind_phys), intent(inout)    :: diffuse_visible_surface_albedo(:)
    real(kind=kind_phys), intent(inout)    :: diffuse_nir_surface_albedo(:)
    real(kind=kind_phys), intent(in)     :: latitude(:)
    real(kind=kind_phys), intent(in)     :: instantaneous_cosine_of_zenith_angle(:)
    real(kind=kind_phys), intent(in)     :: surface_downwelling_shortwave_flux(:)
    real(kind=kind_phys), intent(in)     :: surface_downwelling_longwave_flux(:)
    real(kind=kind_phys), intent(in) :: direct_nir_shortwave_flux_at_surface(:)
    real(kind=kind_phys), intent(in) :: diffuse_nir_shortwave_flux_at_surface(:)
    real(kind=kind_phys), intent(in) :: direct_vis_shortwave_flux_at_surface(:)
    real(kind=kind_phys), intent(in) :: diffuse_vis_shortwave_flux_at_surface(:)
    real(kind=kind_phys), intent(in) :: snicar_hydrophobic_black_carbon_deposition_flux(:)
    real(kind=kind_phys), intent(in) :: snicar_hydrophilic_black_carbon_deposition_flux(:)
    real(kind=kind_phys), intent(in) :: snicar_hydrophobic_organic_carbon_deposition_flux(:)
    real(kind=kind_phys), intent(in) :: snicar_hydrophilic_organic_carbon_deposition_flux(:)
    real(kind=kind_phys), intent(in) :: snicar_dust_species_1_deposition_flux(:)
    real(kind=kind_phys), intent(in) :: snicar_dust_species_2_deposition_flux(:)
    real(kind=kind_phys), intent(in) :: snicar_dust_species_3_deposition_flux(:)
    real(kind=kind_phys), intent(in) :: snicar_dust_species_4_deposition_flux(:)
    real(kind=kind_phys), intent(in) :: snicar_dust_species_5_deposition_flux(:)
    real(kind=kind_phys), intent(in)     :: surface_air_pressure(:)
    real(kind=kind_phys), intent(in)     :: air_pressure_at_surface_adjacent_layer(:)
    real(kind=kind_phys), intent(in)     :: deep_soil_temperature(:)
    real(kind=kind_phys), intent(in)     :: air_temperature_at_lowest_model_layer(:)
    real(kind=kind_phys), intent(in)     :: specific_humidity_at_lowest_model_layer(:)
    real(kind=kind_phys), intent(in)     :: eastward_wind_at_lowest_model_layer(:)
    real(kind=kind_phys), intent(in)     :: northward_wind_at_lowest_model_layer(:)
    real(kind=kind_phys), intent(in), optional :: convective_precipitation_rate_on_previous_timestep(:)
    real(kind=kind_phys), intent(in), optional :: explicit_precipitation_rate_on_previous_timestep(:)
    real(kind=kind_phys), intent(in) :: precipitation_amount_on_dynamics_timestep(:)
    real(kind=kind_phys), intent(in), optional :: snowfall_rate_on_previous_timestep(:)
    real(kind=kind_phys), intent(in), optional :: graupel_precipitation_rate_on_previous_timestep(:)
    real(kind=kind_phys), intent(in), optional :: ice_precipitation_rate_on_previous_timestep(:)
    real(kind=kind_phys), intent(in)     :: precipitation_type(:)
    logical, intent(in)                    :: flag_for_iteration(:)
    logical, intent(in)                    :: flag_nonzero_land_surface_fraction(:)
    logical, intent(in), optional          :: apply_urban_irrigation(:)
    character(len=*), intent(out)          :: errmsg
    integer, intent(out)                   :: errflg
    integer                                :: n
    integer                                :: soil_layer
    integer                                :: number_of_active_columns
    integer                                :: active_source_column(horizontal_loop_extent)
    logical                                :: packed_apply_urban_irrigation(horizontal_loop_extent)
    real(kind=kind_phys)                   :: shortwave_component_total(horizontal_loop_extent)
    real(kind=kind_phys)                   :: air_density(horizontal_loop_extent)
    real(kind=kind_phys)                   :: runoff_accumulation_time(horizontal_loop_extent)

    errmsg = ''
    errflg = 0

    if (.not. allocated(NoahmpContexts)) then
       errflg = 1
       errmsg = 'noahmp_run called before noahmp_init'
       return
    end if
    if (instance_number < 1 .or. instance_number > size(NoahmpContexts)) then
       errflg = 1
       errmsg = 'noahmp_run: instance_number is outside the allocated context range'
       return
    end if
    if (.not. allocated(NoahmpContexts(instance_number)%io)) then
       errflg = 1
       errmsg = 'noahmp_run called before noahmp_init for this instance'
       return
    end if

    associate (NoahmpIO => NoahmpContexts(instance_number)%io)

    if (horizontal_loop_extent < 0) then
       errflg = 1
       errmsg = 'NoahmpModular_run: horizontal_loop_extent must not be negative'
       return
    end if

    if (horizontal_loop_extent > NoahmpIO%capacity) then
       errflg = 1
       write(errmsg, '(a,i0,a,i0)') &
             'NoahmpModular_run: horizontal_loop_extent=', &
             horizontal_loop_extent, ' exceeds internal NoahmpIO capacity=', NoahmpIO%capacity
       return
    end if

    if (size(flag_for_iteration) < horizontal_loop_extent .or. &
        size(flag_nonzero_land_surface_fraction) < horizontal_loop_extent .or. &
        size(apply_urban_irrigation) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: execution-control masks are too small'
       return
    end if

    number_of_active_columns = 0
    do n = 1, horizontal_loop_extent
       if (flag_for_iteration(n) .and. &
           flag_nonzero_land_surface_fraction(n)) then
          number_of_active_columns = number_of_active_columns + 1
          active_source_column(number_of_active_columns) = n
          packed_apply_urban_irrigation(number_of_active_columns) = &
               apply_urban_irrigation(n)
       end if
    end do

    NoahmpIO%ncol = number_of_active_columns
    if (number_of_active_columns == 0) return

    if (timestep_for_physics <= 0) then
       errflg = 1
       errmsg = 'NoahmpModular_run: timestep_for_physics must be positive'
       return
    end if

    if (number_of_days_in_current_year <= 0) then
       errflg = 1
       errmsg = 'NoahmpModular_run: number_of_days_in_current_year must be positive'
       return
    end if

    if (vertical_dimension_of_soil <= 0 .or. &
        vertical_dimension_of_surface_snow <= 0 .or. &
        number_of_shortwave_radiation_bands <= 0) then
       errflg = 1
       errmsg = 'NoahmpModular_run: soil, snow, and radiation dimensions must be positive'
       return
    end if

    if (vertical_dimension_of_soil /= NoahmpIO%NSOIL .or. &
        vertical_dimension_of_surface_snow /= NoahmpIO%NSNOW .or. &
        number_of_shortwave_radiation_bands /= NoahmpIO%NUMRAD) then
       errflg = 1
       errmsg = 'NoahmpModular_run: CCPP dimensions disagree with NoahmpIO allocation dimensions'
       return
    end if

    if (lower_bound_of_vertical_dimension_of_surface_snow /= &
             1-vertical_dimension_of_surface_snow .or. &
        upper_bound_of_vertical_dimension_of_surface_snow /= 0) then
       errflg = 1
       errmsg = 'NoahmpModular_run: CCPP snow bounds must be 1-NSNOW:0'
       return
    end if

    if (size(depth_of_soil_layer_interfaces) < vertical_dimension_of_soil) then
       errflg = 1
       errmsg = 'NoahmpModular_run: soil-interface depth array is too small'
       return
    end if

    if (size(characteristic_grid_lengthscale) < horizontal_loop_extent .or. &
        size(vegetation_category) < horizontal_loop_extent .or. &
        size(soil_category) < horizontal_loop_extent .or. &
        size(soil_color_classification) < horizontal_loop_extent .or. &
        size(surface_slope_classification) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: CCPP geometry or surface-category arrays are too small'
       return
    end if

    if (.not. present(active_snow_layer_lower_index) .or. &
        .not. present(depth_of_snow_soil_layer_interfaces) .or. &
        .not. present(leaf_area_index) .or. &
        .not. present(stem_area_index) .or. &
        .not. present(ground_temperature) .or. &
        .not. present(canopy_temperature) .or. &
        .not. present(dimensionless_snow_age) .or. &
        .not. present(snow_albedo_on_previous_timestep) .or. &
        .not. present(canopy_air_vapor_pressure) .or. &
        .not. present(canopy_air_temperature) .or. &
        .not. present(snow_temperature) .or. &
        .not. present(canopy_liquid_water) .or. &
        .not. present(canopy_ice_water) .or. &
        .not. present(canopy_wet_fraction) .or. &
        .not. present(snow_water_equivalent_on_previous_timestep) .or. &
        .not. present(snow_layer_ice) .or. &
        .not. present(snow_layer_liquid_water)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: required conditional Noah-MP state is unavailable'
       return
    end if

    if (.not. present(surface_exchange_coefficient_for_heat) .or. &
        .not. present(surface_exchange_coefficient_for_momentum) .or. &
        .not. present(water_table_depth) .or. &
        .not. present(soil_moisture_below_soil_column) .or. &
        .not. present(aquifer_water_storage) .or. &
        .not. present(soil_aquifer_water_storage) .or. &
        .not. present(lake_water_storage) .or. &
        .not. present(leaf_mass) .or. &
        .not. present(root_mass) .or. &
        .not. present(stem_mass) .or. &
        .not. present(wood_mass) .or. &
        .not. present(deep_soil_carbon_mass) .or. &
        .not. present(shallow_soil_carbon_mass) .or. &
        .not. present(ground_snowfall_rate) .or. &
        .not. present(canopy_resistance) .or. &
        .not. present(temperature_at_2m_from_noahmp) .or. &
        .not. present(specific_humidity_at_2m_from_noahmp) .or. &
        .not. present(normalized_soil_wetness)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: required legacy Noah-MP coupling state is unavailable'
       return
    end if

    if (.not. present(flood_irrigation_fraction) .or. &
        .not. present(flood_irrigation_water_amount) .or. &
        .not. present(micro_irrigation_fraction) .or. &
        .not. present(micro_irrigation_water_amount) .or. &
        .not. present(sprinkler_irrigation_fraction) .or. &
        .not. present(sprinkler_irrigation_water_amount) .or. &
        .not. present(tile_drainage_fraction) .or. &
        .not. present(grid_irrigation_fraction) .or. &
        .not. present(sprinkler_irrigation_event_count) .or. &
        .not. present(micro_irrigation_event_count) .or. &
        .not. present(flood_irrigation_event_count) .or. &
        .not. present(soil_saturated_fraction) .or. &
        .not. present(wetland_water_storage) .or. &
        .not. present(equilibrium_soil_moisture) .or. &
        .not. present(apply_urban_irrigation)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: required SCM water-management state is unavailable'
       return
    end if

    if (size(active_snow_layer_lower_index) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: active snow-layer-index array is too small'
       return
    end if

    if (size(depth_of_snow_soil_layer_interfaces, 1) < horizontal_loop_extent .or. &
        lbound(depth_of_snow_soil_layer_interfaces, 2) > &
             lower_bound_of_vertical_dimension_of_surface_snow .or. &
        ubound(depth_of_snow_soil_layer_interfaces, 2) < vertical_dimension_of_soil) then
       errflg = 1
       errmsg = 'NoahmpModular_run: snow-soil interface-depth array is too small'
       return
    end if

    if (size(atmospheric_forcing_height) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: atmospheric-forcing-height output is too small'
       return
    end if

    if (size(leaf_area_index) < horizontal_loop_extent .or. &
        size(stem_area_index) < horizontal_loop_extent .or. &
        size(specific_humidity_at_surface) < horizontal_loop_extent .or. &
        size(ground_temperature) < horizontal_loop_extent .or. &
        size(canopy_temperature) < horizontal_loop_extent .or. &
        size(dimensionless_snow_age) < horizontal_loop_extent .or. &
        size(snow_albedo_on_previous_timestep) < horizontal_loop_extent .or. &
        size(canopy_air_vapor_pressure) < horizontal_loop_extent .or. &
        size(canopy_air_temperature) < horizontal_loop_extent .or. &
        size(surface_exchange_coefficient_for_heat) < horizontal_loop_extent .or. &
        size(surface_exchange_coefficient_for_momentum) < horizontal_loop_extent .or. &
        size(maximum_vegetation_fraction) < horizontal_loop_extent .or. &
        size(vegetation_fraction) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: scalar energy-state arrays are too small'
       return
    end if

    if (size(soil_temperature, 1) < horizontal_loop_extent .or. &
        size(soil_temperature, 2) < vertical_dimension_of_soil .or. &
        size(snow_temperature, 1) < horizontal_loop_extent .or. &
        lbound(snow_temperature, 2) > lower_bound_of_vertical_dimension_of_surface_snow .or. &
        ubound(snow_temperature, 2) < upper_bound_of_vertical_dimension_of_surface_snow) then
       errflg = 1
       errmsg = 'NoahmpModular_run: layered energy-state arrays are too small'
       return
    end if

    if (size(canopy_liquid_water) < horizontal_loop_extent .or. &
        size(canopy_ice_water) < horizontal_loop_extent .or. &
        size(canopy_wet_fraction) < horizontal_loop_extent .or. &
        size(snow_water_equivalent) < horizontal_loop_extent .or. &
        size(snow_water_equivalent_on_previous_timestep) < horizontal_loop_extent .or. &
        size(snow_depth) < horizontal_loop_extent .or. &
        size(flood_irrigation_fraction) < horizontal_loop_extent .or. &
        size(flood_irrigation_water_amount) < horizontal_loop_extent .or. &
        size(micro_irrigation_fraction) < horizontal_loop_extent .or. &
        size(micro_irrigation_water_amount) < horizontal_loop_extent .or. &
        size(sprinkler_irrigation_fraction) < horizontal_loop_extent .or. &
        size(sprinkler_irrigation_water_amount) < horizontal_loop_extent .or. &
        size(water_table_depth) < horizontal_loop_extent .or. &
        size(soil_moisture_below_soil_column) < horizontal_loop_extent .or. &
        size(tile_drainage_fraction) < horizontal_loop_extent .or. &
        size(aquifer_water_storage) < horizontal_loop_extent .or. &
        size(soil_aquifer_water_storage) < horizontal_loop_extent .or. &
        size(lake_water_storage) < horizontal_loop_extent .or. &
        size(grid_irrigation_fraction) < horizontal_loop_extent .or. &
        size(sprinkler_irrigation_event_count) < horizontal_loop_extent .or. &
        size(micro_irrigation_event_count) < horizontal_loop_extent .or. &
        size(flood_irrigation_event_count) < horizontal_loop_extent .or. &
        size(soil_saturated_fraction) < horizontal_loop_extent .or. &
        size(wetland_water_storage) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: scalar water-state arrays are too small'
       return
    end if

    if (size(snow_layer_ice, 1) < horizontal_loop_extent .or. &
        lbound(snow_layer_ice, 2) > lower_bound_of_vertical_dimension_of_surface_snow .or. &
        ubound(snow_layer_ice, 2) < upper_bound_of_vertical_dimension_of_surface_snow .or. &
        size(snow_layer_liquid_water, 1) < horizontal_loop_extent .or. &
        lbound(snow_layer_liquid_water, 2) > &
             lower_bound_of_vertical_dimension_of_surface_snow .or. &
        ubound(snow_layer_liquid_water, 2) < &
             upper_bound_of_vertical_dimension_of_surface_snow .or. &
        size(soil_liquid_water, 1) < horizontal_loop_extent .or. &
        size(soil_liquid_water, 2) < vertical_dimension_of_soil .or. &
        size(soil_moisture, 1) < horizontal_loop_extent .or. &
        size(soil_moisture, 2) < vertical_dimension_of_soil .or. &
        size(equilibrium_soil_moisture, 1) < horizontal_loop_extent .or. &
        size(equilibrium_soil_moisture, 2) < vertical_dimension_of_soil) then
       errflg = 1
       errmsg = 'NoahmpModular_run: layered water-state arrays are too small'
       return
    end if


    if (size(leaf_mass) < horizontal_loop_extent .or. &
        size(root_mass) < horizontal_loop_extent .or. &
        size(stem_mass) < horizontal_loop_extent .or. &
        size(wood_mass) < horizontal_loop_extent .or. &
        size(deep_soil_carbon_mass) < horizontal_loop_extent .or. &
        size(shallow_soil_carbon_mass) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: biochemistry and crop arrays are too small'
       return
    end if



    if ( size(surface_radiative_temperature) < horizontal_loop_extent .or. &
        size(surface_emissivity) < horizontal_loop_extent .or. &
        size(surface_roughness_length) < horizontal_loop_extent .or. &
        size(temperature_at_2m_from_noahmp) < horizontal_loop_extent .or. &
        size(specific_humidity_at_2m_from_noahmp) < horizontal_loop_extent .or. &
        size(surface_drag_wind_speed_for_momentum) < horizontal_loop_extent .or. &
        size(surface_drag_mass_flux_for_heat_and_moisture) < horizontal_loop_extent .or. &
        size(surface_albedo) < horizontal_loop_extent .or. &
        size(canopy_resistance) < horizontal_loop_extent .or. &
        size(latent_heat_flux) < horizontal_loop_extent .or. &
        size(sensible_heat_flux) < horizontal_loop_extent .or. &
        size(ground_heat_flux) < horizontal_loop_extent .or. &
        size(soil_upward_latent_heat_flux) < horizontal_loop_extent .or. &
        size(transpiration_latent_heat_flux) < horizontal_loop_extent .or. &
        size(canopy_upward_latent_heat_flux) < horizontal_loop_extent .or. &
        size(precipitation_advected_heat_flux) < horizontal_loop_extent .or. &
        size(soil_moisture_content) < horizontal_loop_extent .or. &
        size(normalized_soil_wetness) < horizontal_loop_extent .or. &
        size(wilting_soil_moisture) < horizontal_loop_extent .or. &
        size(reference_soil_moisture) < horizontal_loop_extent .or. &
        size(snow_cover_fraction) < horizontal_loop_extent .or. &
        size(total_canopy_water) < horizontal_loop_extent .or. &
        size(ground_snowfall_rate) < horizontal_loop_extent .or. &
        size(surface_runoff) < horizontal_loop_extent .or. &
        size(surface_runoff_flux) < horizontal_loop_extent .or. &
        size(subsurface_runoff) < horizontal_loop_extent .or. &
        size(net_canopy_evaporation) < horizontal_loop_extent .or. &
        size(net_ground_evaporation) < horizontal_loop_extent .or. &
        size(transpiration) < horizontal_loop_extent .or. &
        size(snow_deposition_sublimation_upward_latent_heat_flux) < horizontal_loop_extent .or. &
        size(snow_freezing_rain_upward_latent_heat_flux) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: scalar diagnostic output arrays are too small'
       return
    end if

    if (size(direct_visible_surface_albedo) < horizontal_loop_extent .or. &
        size(direct_nir_surface_albedo) < horizontal_loop_extent .or. &
        size(diffuse_visible_surface_albedo) < horizontal_loop_extent .or. &
        size(diffuse_nir_surface_albedo) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: surface albedo output arrays are too small'
       return
    end if


    if (size(latitude) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: latitude is smaller than horizontal_loop_extent'
       return
    end if

    if (size(instantaneous_cosine_of_zenith_angle) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: cosine of zenith angle is too small'
       return
    end if


    if (size(surface_downwelling_shortwave_flux) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: shortwave forcing is smaller than horizontal_loop_extent'
       return
    end if

    if (size(surface_downwelling_longwave_flux) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: longwave forcing is smaller than horizontal_loop_extent'
       return
    end if

    if (size(direct_nir_shortwave_flux_at_surface) < &
             horizontal_loop_extent .or. &
        size(diffuse_nir_shortwave_flux_at_surface) < &
             horizontal_loop_extent .or. &
        size(direct_vis_shortwave_flux_at_surface) < &
             horizontal_loop_extent .or. &
        size(diffuse_vis_shortwave_flux_at_surface) < &
             horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: partitioned shortwave forcing is too small'
       return
    end if

    if (size(surface_air_pressure) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: surface pressure is smaller than horizontal_loop_extent'
       return
    end if

    if (size(air_pressure_at_surface_adjacent_layer) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: lowest-layer pressure is smaller than horizontal_loop_extent'
       return
    end if

    if (size(deep_soil_temperature) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: deep-soil temperature is smaller than horizontal_loop_extent'
       return
    end if

    if (size(air_temperature_at_lowest_model_layer) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: air temperature is smaller than horizontal_loop_extent'
       return
    end if

    if (size(specific_humidity_at_lowest_model_layer) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: specific humidity is smaller than horizontal_loop_extent'
       return
    end if

    if (size(eastward_wind_at_lowest_model_layer) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: eastward wind is smaller than horizontal_loop_extent'
       return
    end if

    if (size(northward_wind_at_lowest_model_layer) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: northward wind is smaller than horizontal_loop_extent'
       return
    end if

    if (.not. present(convective_precipitation_rate_on_previous_timestep) .or. &
        .not. present(explicit_precipitation_rate_on_previous_timestep) .or. &
        .not. present(snowfall_rate_on_previous_timestep) .or. &
        .not. present(graupel_precipitation_rate_on_previous_timestep) .or. &
        .not. present(ice_precipitation_rate_on_previous_timestep)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: required conditional precipitation state is unavailable'
       return
    end if

    if (size(convective_precipitation_rate_on_previous_timestep) < &
        horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: convective precipitation is smaller than horizontal_loop_extent'
       return
    end if

    if (size(explicit_precipitation_rate_on_previous_timestep) < &
        horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: nonconvective precipitation is smaller than horizontal_loop_extent'
       return
    end if


    if (size(precipitation_amount_on_dynamics_timestep) < &
             horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: precipitation amount is too small'
       return
    end if

    if (size(snowfall_rate_on_previous_timestep) < &
        horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: snow precipitation rate is too small'
       return
    end if

    if (size(graupel_precipitation_rate_on_previous_timestep) < &
        horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: graupel precipitation rate is too small'
       return
    end if

    if (size(ice_precipitation_rate_on_previous_timestep) < &
         horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: ice precipitation rate is too small'
       return
    end if

    if (size(precipitation_type) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: precipitation type array is too small'
       return
    end if

    if (.not. allocated(NoahmpIO%SWDOWN) .or. .not. allocated(NoahmpIO%GLW)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal radiation forcing arrays are not allocated'
       return
    end if

    if (size(NoahmpIO%SWDOWN) < horizontal_loop_extent .or. &
        size(NoahmpIO%GLW) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal radiation forcing arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%RadSwDirFrac) .or. &
        .not. allocated(NoahmpIO%RadSwVisFrac)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal shortwave-fraction arrays are not allocated'
       return
    end if

    if (size(NoahmpIO%RadSwDirFrac) < horizontal_loop_extent .or. &
        size(NoahmpIO%RadSwVisFrac) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal shortwave-fraction arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%PS) .or. .not. allocated(NoahmpIO%PRSL1)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal pressure forcing arrays are not allocated'
       return
    end if

    if (size(NoahmpIO%PS) < horizontal_loop_extent .or. &
        size(NoahmpIO%PRSL1) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal pressure forcing arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%TMN)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal deep-soil temperature array is not allocated'
       return
    end if

    if (size(NoahmpIO%TMN) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal deep-soil temperature array is too small'
       return
    end if

    if (.not. allocated(NoahmpIO%MP_RAINC)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal convective-precipitation array is not allocated'
       return
    end if

    if (size(NoahmpIO%MP_RAINC) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal convective-precipitation array is too small'
       return
    end if

    if (.not. allocated(NoahmpIO%MP_RAINNC)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal nonconvective-precipitation array is not allocated'
       return
    end if

    if (size(NoahmpIO%MP_RAINNC) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal nonconvective-precipitation array is too small'
       return
    end if

    if (.not. allocated(NoahmpIO%MP_SHCV)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal shallow-precipitation array is not allocated'
       return
    end if

    if (size(NoahmpIO%MP_SHCV) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal shallow-precipitation array is too small'
       return
    end if

    if (.not. allocated(NoahmpIO%RAINBL) .or. &
        .not. allocated(NoahmpIO%MP_SNOW) .or. &
        .not. allocated(NoahmpIO%MP_GRAUP) .or. &
        .not. allocated(NoahmpIO%MP_HAIL) .or. &
        .not. allocated(NoahmpIO%SR)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal precipitation arrays are not allocated'
       return
    end if

    if (size(NoahmpIO%RAINBL) < horizontal_loop_extent .or. &
        size(NoahmpIO%MP_SNOW) < horizontal_loop_extent .or. &
        size(NoahmpIO%MP_GRAUP) < horizontal_loop_extent .or. &
        size(NoahmpIO%MP_HAIL) < horizontal_loop_extent .or. &
        size(NoahmpIO%SR) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal precipitation arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%T_PHY) .or. &
        .not. allocated(NoahmpIO%QV_CURR) .or. &
        .not. allocated(NoahmpIO%U_PHY) .or. &
        .not. allocated(NoahmpIO%V_PHY)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal temperature/humidity/wind arrays are not allocated'
       return
    end if

    NoahmpIO%KTS = 1
    if (NoahmpIO%KTS < lbound(NoahmpIO%T_PHY, 1) .or. &
        NoahmpIO%KTS > ubound(NoahmpIO%T_PHY, 1) .or. &
        NoahmpIO%KTS < lbound(NoahmpIO%QV_CURR, 1) .or. &
        NoahmpIO%KTS > ubound(NoahmpIO%QV_CURR, 1) .or. &
        NoahmpIO%KTS < lbound(NoahmpIO%U_PHY, 1) .or. &
        NoahmpIO%KTS > ubound(NoahmpIO%U_PHY, 1) .or. &
        NoahmpIO%KTS < lbound(NoahmpIO%V_PHY, 1) .or. &
        NoahmpIO%KTS > ubound(NoahmpIO%V_PHY, 1)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: KTS is outside internal forcing-array bounds'
       return
    end if

    if (size(NoahmpIO%T_PHY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%QV_CURR, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%U_PHY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%V_PHY, 2) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal temperature/humidity/wind arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%PBLH) .or. &
        size(NoahmpIO%PBLH) < horizontal_loop_extent .or. &
        size(atmosphere_boundary_layer_thickness) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: boundary-layer-height array is unavailable or too small'
       return
    end if

    if (.not. allocated(NoahmpIO%XLAT) .or. &
        .not. allocated(NoahmpIO%COSZEN) .or. &
        .not. allocated(NoahmpIO%DZ8W)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal geometry arrays are not allocated'
       return
    end if

    if (NoahmpIO%KTS < lbound(NoahmpIO%DZ8W, 1) .or. &
        NoahmpIO%KTS > ubound(NoahmpIO%DZ8W, 1)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: KTS is outside DZ8W bounds'
       return
    end if

    if (size(NoahmpIO%XLAT) < horizontal_loop_extent .or. &
        size(NoahmpIO%COSZEN) < horizontal_loop_extent .or. &
        size(NoahmpIO%DZ8W, 2) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal geometry arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%DX) .or. &
        .not. allocated(NoahmpIO%DY) .or. &
        .not. allocated(NoahmpIO%ZSOIL) .or. &
        .not. allocated(NoahmpIO%IVGTYP) .or. &
        .not. allocated(NoahmpIO%CROPCAT) .or. &
        .not. allocated(NoahmpIO%ISLTYP) .or. &
        .not. allocated(NoahmpIO%SOILCOL) .or. &
        .not. allocated(NoahmpIO%ICE) .or. &
        .not. allocated(NoahmpIO%XICE) .or. &
        .not. allocated(NoahmpIO%SOILCL1) .or. &
        .not. allocated(NoahmpIO%SOILCL2) .or. &
        .not. allocated(NoahmpIO%SOILCL3) .or. &
        .not. allocated(NoahmpIO%SOILCL4)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal soil and surface-category arrays are not allocated'
       return
    end if

    if (size(NoahmpIO%DX) < horizontal_loop_extent .or. &
        size(NoahmpIO%DY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ZSOIL) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%IVGTYP) < horizontal_loop_extent .or. &
        size(NoahmpIO%CROPCAT) < horizontal_loop_extent .or. &
        size(NoahmpIO%ISLTYP) < horizontal_loop_extent .or. &
        size(NoahmpIO%SOILCOL) < horizontal_loop_extent .or. &
        size(NoahmpIO%ICE) < horizontal_loop_extent .or. &
        size(NoahmpIO%XICE) < horizontal_loop_extent .or. &
        size(NoahmpIO%SOILCL1) < horizontal_loop_extent .or. &
        size(NoahmpIO%SOILCL2) < horizontal_loop_extent .or. &
        size(NoahmpIO%SOILCL3) < horizontal_loop_extent .or. &
        size(NoahmpIO%SOILCL4) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal soil and surface-category arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%source_i) .or. &
        .not. allocated(NoahmpIO%source_j) .or. &
        .not. allocated(NoahmpIO%ISNOWXY) .or. &
        .not. allocated(NoahmpIO%ZSNSOXY) .or. &
        .not. allocated(NoahmpIO%FORCZLSM)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal grid and snow-configuration arrays are not allocated'
       return
    end if

    if (size(NoahmpIO%source_i) < horizontal_loop_extent .or. &
        size(NoahmpIO%source_j) < horizontal_loop_extent .or. &
        size(NoahmpIO%ISNOWXY) < horizontal_loop_extent .or. &
        lbound(NoahmpIO%ZSNSOXY, 1) > -vertical_dimension_of_surface_snow+1 .or. &
        ubound(NoahmpIO%ZSNSOXY, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%ZSNSOXY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%FORCZLSM) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal grid and snow-configuration arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%LAI) .or. &
        .not. allocated(NoahmpIO%XSAIXY) .or. &
        .not. allocated(NoahmpIO%QSFC) .or. &
        .not. allocated(NoahmpIO%TGXY) .or. &
        .not. allocated(NoahmpIO%TVXY) .or. &
        .not. allocated(NoahmpIO%TAUSSXY) .or. &
        .not. allocated(NoahmpIO%ALBOLDXY) .or. &
        .not. allocated(NoahmpIO%EAHXY) .or. &
        .not. allocated(NoahmpIO%TAHXY) .or. &
        .not. allocated(NoahmpIO%CHXY) .or. &
        .not. allocated(NoahmpIO%CMXY) .or. &
        .not. allocated(NoahmpIO%ACC_SSOILXY) .or. &
        .not. allocated(NoahmpIO%GVFMAX) .or. &
        .not. allocated(NoahmpIO%VEGFRA) .or. &
        .not. allocated(NoahmpIO%IRRSPLH) .or. &
        .not. allocated(NoahmpIO%TSLB) .or. &
        .not. allocated(NoahmpIO%TSNOXY) .or. &
        .not. allocated(NoahmpIO%ALBSOILDIRXY) .or. &
        .not. allocated(NoahmpIO%ALBSOILDIFXY)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal energy-state arrays are not allocated'
       return
    end if

    if (size(NoahmpIO%LAI) < horizontal_loop_extent .or. &
        size(NoahmpIO%XSAIXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QSFC) < horizontal_loop_extent .or. &
        size(NoahmpIO%TGXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%TVXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%TAUSSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ALBOLDXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%EAHXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%TAHXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CHXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CMXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_SSOILXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%GVFMAX) < horizontal_loop_extent .or. &
        size(NoahmpIO%VEGFRA) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRRSPLH) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal scalar energy-state arrays are too small'
       return
    end if

    if (size(NoahmpIO%TSLB, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%TSLB, 2) < horizontal_loop_extent .or. &
        lbound(NoahmpIO%TSNOXY, 1) > -vertical_dimension_of_surface_snow+1 .or. &
        ubound(NoahmpIO%TSNOXY, 1) < 0 .or. &
        size(NoahmpIO%TSNOXY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%ALBSOILDIRXY, 1) < number_of_shortwave_radiation_bands .or. &
        size(NoahmpIO%ALBSOILDIRXY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%ALBSOILDIFXY, 1) < number_of_shortwave_radiation_bands .or. &
        size(NoahmpIO%ALBSOILDIFXY, 2) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal layered energy-state arrays are too small'
       return
    end if

    if (.not. allocated(NoahmpIO%CANLIQXY) .or. &
        .not. allocated(NoahmpIO%CANICEXY) .or. &
        .not. allocated(NoahmpIO%FWETXY) .or. &
        .not. allocated(NoahmpIO%SNOW) .or. &
        .not. allocated(NoahmpIO%SNEQVOXY) .or. &
        .not. allocated(NoahmpIO%SNOWH) .or. &
        .not. allocated(NoahmpIO%FIFRACT) .or. &
        .not. allocated(NoahmpIO%IRWATFI) .or. &
        .not. allocated(NoahmpIO%MIFRACT) .or. &
        .not. allocated(NoahmpIO%IRWATMI) .or. &
        .not. allocated(NoahmpIO%SIFRACT) .or. &
        .not. allocated(NoahmpIO%IRWATSI) .or. &
        .not. allocated(NoahmpIO%ZWTXY) .or. &
        .not. allocated(NoahmpIO%SMCWTDXY) .or. &
        .not. allocated(NoahmpIO%TD_FRACTION) .or. &
        .not. allocated(NoahmpIO%WAXY) .or. &
        .not. allocated(NoahmpIO%WTXY) .or. &
        .not. allocated(NoahmpIO%WSLAKEXY) .or. &
        .not. allocated(NoahmpIO%IRFRACT) .or. &
        .not. allocated(NoahmpIO%IRNUMSI) .or. &
        .not. allocated(NoahmpIO%IRNUMMI) .or. &
        .not. allocated(NoahmpIO%IRNUMFI) .or. &
        .not. allocated(NoahmpIO%FSATXY) .or. &
        .not. allocated(NoahmpIO%WSURFXY) .or. &
        .not. allocated(NoahmpIO%SNICEXY) .or. &
        .not. allocated(NoahmpIO%SNLIQXY) .or. &
        .not. allocated(NoahmpIO%SH2O) .or. &
        .not. allocated(NoahmpIO%SMOIS) .or. &
        .not. allocated(NoahmpIO%SMOISEQ) .or. &
        .not. allocated(NoahmpIO%ACC_QSEVAXY) .or. &
        .not. allocated(NoahmpIO%ACC_QINSURXY) .or. &
        .not. allocated(NoahmpIO%ACC_DWATERXY) .or. &
        .not. allocated(NoahmpIO%ACC_PRCPXY) .or. &
        .not. allocated(NoahmpIO%ACC_ECANXY) .or. &
        .not. allocated(NoahmpIO%ACC_ETRANXY) .or. &
        .not. allocated(NoahmpIO%ACC_EDIRXY) .or. &
        .not. allocated(NoahmpIO%ACC_ETRANIXY) .or. &
        .not. allocated(NoahmpIO%ACC_GLAFLWXY) .or. &
        .not. allocated(NoahmpIO%SFCRUNOFF) .or. &
        .not. allocated(NoahmpIO%UDRUNOFF) .or. &
        .not. allocated(NoahmpIO%QTDRAIN) .or. &
        .not. allocated(NoahmpIO%ACSNOW) .or. &
        .not. allocated(NoahmpIO%ACSNOM) .or. &
        .not. allocated(NoahmpIO%RECHXY) .or. &
        .not. allocated(NoahmpIO%DEEPRECHXY) .or. &
        .not. allocated(NoahmpIO%IRSIVOL) .or. &
        .not. allocated(NoahmpIO%IRMIVOL) .or. &
        .not. allocated(NoahmpIO%IRFIVOL) .or. &
        .not. allocated(NoahmpIO%IRELOSS)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal persistent water arrays are not allocated'
       return
    end if

    if (size(NoahmpIO%CANLIQXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CANICEXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%FWETXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SNOW) < horizontal_loop_extent .or. &
        size(NoahmpIO%SNEQVOXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SNOWH) < horizontal_loop_extent .or. &
        size(NoahmpIO%FIFRACT) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRWATFI) < horizontal_loop_extent .or. &
        size(NoahmpIO%MIFRACT) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRWATMI) < horizontal_loop_extent .or. &
        size(NoahmpIO%SIFRACT) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRWATSI) < horizontal_loop_extent .or. &
        size(NoahmpIO%ZWTXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMCWTDXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%TD_FRACTION) < horizontal_loop_extent .or. &
        size(NoahmpIO%WAXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%WTXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%WSLAKEXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRFRACT) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRNUMSI) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRNUMMI) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRNUMFI) < horizontal_loop_extent .or. &
        size(NoahmpIO%FSATXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%WSURFXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_QSEVAXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_QINSURXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_DWATERXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_PRCPXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_ECANXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_ETRANXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_EDIRXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_GLAFLWXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SFCRUNOFF) < horizontal_loop_extent .or. &
        size(NoahmpIO%UDRUNOFF) < horizontal_loop_extent .or. &
        size(NoahmpIO%QTDRAIN) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACSNOW) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACSNOM) < horizontal_loop_extent .or. &
        size(NoahmpIO%RECHXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%DEEPRECHXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRSIVOL) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRMIVOL) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRFIVOL) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRELOSS) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal scalar water arrays are too small'
       return
    end if

    if (lbound(NoahmpIO%SNICEXY, 1) > -vertical_dimension_of_surface_snow+1 .or. &
        ubound(NoahmpIO%SNICEXY, 1) < 0 .or. &
        size(NoahmpIO%SNICEXY, 2) < horizontal_loop_extent .or. &
        lbound(NoahmpIO%SNLIQXY, 1) > -vertical_dimension_of_surface_snow+1 .or. &
        ubound(NoahmpIO%SNLIQXY, 1) < 0 .or. &
        size(NoahmpIO%SNLIQXY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%SH2O, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%SH2O, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMOIS, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%SMOIS, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMOISEQ, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%SMOISEQ, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%ACC_ETRANIXY, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%ACC_ETRANIXY, 2) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal layered water arrays are too small'
       return
    end if

	if (.not. allocated(NoahmpIO%PGSXY) .or. .not. allocated(NoahmpIO%LFMASSXY) .or. &
		.not. allocated(NoahmpIO%RTMASSXY) .or. .not. allocated(NoahmpIO%STMASSXY) .or. &
		.not. allocated(NoahmpIO%WOODXY) .or. .not. allocated(NoahmpIO%STBLCPXY) .or. &
		.not. allocated(NoahmpIO%FASTCPXY) .or. .not. allocated(NoahmpIO%GRAINXY) .or. &
		.not. allocated(NoahmpIO%GDDXY) .or. .not. allocated(NoahmpIO%PLANTING) .or. &
		.not. allocated(NoahmpIO%HARVEST) .or. .not. allocated(NoahmpIO%SEASON_GDD) .or. &
		.not. allocated(NoahmpIO%NEEXY) .or. .not. allocated(NoahmpIO%GPPXY) .or. &
		.not. allocated(NoahmpIO%NPPXY) .or. .not. allocated(NoahmpIO%PSNXY)) then
	   errflg = 1
	   errmsg = 'NoahmpModular_run: internal biochemistry and crop arrays are not allocated'
	   return
	end if

    if (size(NoahmpIO%PGSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%LFMASSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%RTMASSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%STMASSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%WOODXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%STBLCPXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%FASTCPXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%GRAINXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%GDDXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%PLANTING) < horizontal_loop_extent .or. &
        size(NoahmpIO%HARVEST) < horizontal_loop_extent .or. &
        size(NoahmpIO%SEASON_GDD) < horizontal_loop_extent .or. &
        size(NoahmpIO%NEEXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%GPPXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%NPPXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%PSNXY) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal biochemistry and crop arrays are too small'
       return
    end if

	if (.not. allocated(NoahmpIO%quartz_3D) .or. .not. allocated(NoahmpIO%BEXP_3D) .or. &
		.not. allocated(NoahmpIO%SMCDRY_3D) .or. .not. allocated(NoahmpIO%SMCWLT_3D) .or. &
		.not. allocated(NoahmpIO%SMCREF_3D) .or. .not. allocated(NoahmpIO%SMCMAX_3D) .or. &
		.not. allocated(NoahmpIO%DKSAT_3D) .or. .not. allocated(NoahmpIO%DWSAT_3D) .or. &
		.not. allocated(NoahmpIO%PSISAT_3D) .or. .not. allocated(NoahmpIO%REFDK_2D) .or. &
		.not. allocated(NoahmpIO%REFKDT_2D) .or. .not. allocated(NoahmpIO%BVIC_2D) .or. &
		.not. allocated(NoahmpIO%AXAJ_2D) .or. .not. allocated(NoahmpIO%BXAJ_2D) .or. &
		.not. allocated(NoahmpIO%XXAJ_2D) .or. .not. allocated(NoahmpIO%BDVIC_2D) .or. &
		.not. allocated(NoahmpIO%GDVIC_2D) .or. .not. allocated(NoahmpIO%BBVIC_2D) .or. &
		.not. allocated(NoahmpIO%IRR_FRAC_2D) .or. .not. allocated(NoahmpIO%IRR_HAR_2D) .or. &
		.not. allocated(NoahmpIO%IRR_LAI_2D) .or. .not. allocated(NoahmpIO%IRR_MAD_2D) .or. &
		.not. allocated(NoahmpIO%FILOSS_2D) .or. .not. allocated(NoahmpIO%SPRIR_RATE_2D) .or. &
		.not. allocated(NoahmpIO%MICIR_RATE_2D) .or. .not. allocated(NoahmpIO%FIRTFAC_2D) .or. &
		.not. allocated(NoahmpIO%IR_RAIN_2D) .or. .not. allocated(NoahmpIO%KLAT_FAC) .or. &
		.not. allocated(NoahmpIO%TDSMC_FAC) .or. .not. allocated(NoahmpIO%TD_DC) .or. &
		.not. allocated(NoahmpIO%TD_DCOEF) .or. .not. allocated(NoahmpIO%TD_DDRAIN) .or. &
		.not. allocated(NoahmpIO%TD_RADI) .or. .not. allocated(NoahmpIO%TD_SPAC) .or. &
		.not. allocated(NoahmpIO%FSATMX) .or. .not. allocated(NoahmpIO%WCAP)) then
	   errflg = 1
	   errmsg = 'NoahmpModular_run: internal spatial-parameter arrays are not allocated'
	   return
	end if

    if (size(NoahmpIO%quartz_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%quartz_3D, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%BEXP_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%BEXP_3D, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMCDRY_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%SMCDRY_3D, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMCWLT_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%SMCWLT_3D, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMCREF_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%SMCREF_3D, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMCMAX_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%SMCMAX_3D, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%DKSAT_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%DKSAT_3D, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%DWSAT_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%DWSAT_3D, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%PSISAT_3D, 1) < vertical_dimension_of_soil .or. &
        size(NoahmpIO%PSISAT_3D, 2) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal layered spatial-parameter arrays are too small'
       return
    end if

    if (size(NoahmpIO%REFDK_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%REFKDT_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%BVIC_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%AXAJ_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%BXAJ_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%XXAJ_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%BDVIC_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%GDVIC_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%BBVIC_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRR_FRAC_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRR_HAR_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRR_LAI_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRR_MAD_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%FILOSS_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%SPRIR_RATE_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%MICIR_RATE_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%FIRTFAC_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%IR_RAIN_2D) < horizontal_loop_extent .or. &
        size(NoahmpIO%KLAT_FAC) < horizontal_loop_extent .or. &
        size(NoahmpIO%TDSMC_FAC) < horizontal_loop_extent .or. &
        size(NoahmpIO%TD_DC) < horizontal_loop_extent .or. &
        size(NoahmpIO%TD_DCOEF) < horizontal_loop_extent .or. &
        size(NoahmpIO%TD_DDRAIN) < horizontal_loop_extent .or. &
        size(NoahmpIO%TD_RADI) < horizontal_loop_extent .or. &
        size(NoahmpIO%TD_SPAC) < horizontal_loop_extent .or. &
        size(NoahmpIO%FSATMX) < horizontal_loop_extent .or. &
        size(NoahmpIO%WCAP) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal scalar spatial-parameter arrays are too small'
       return
    end if

	if (.not. allocated(NoahmpIO%FORCTLSM) .or. .not. allocated(NoahmpIO%FORCQLSM) .or. &
		.not. allocated(NoahmpIO%FORCPLSM) .or. .not. allocated(NoahmpIO%FORCWLSM) .or. &
		.not. allocated(NoahmpIO%TSK) .or. .not. allocated(NoahmpIO%EMISS) .or. &
		.not. allocated(NoahmpIO%Z0) .or. .not. allocated(NoahmpIO%T2MVXY) .or. &
		.not. allocated(NoahmpIO%T2MBXY) .or. .not. allocated(NoahmpIO%Q2MVXY) .or. &
		.not. allocated(NoahmpIO%Q2MBXY) .or. .not. allocated(NoahmpIO%FVEGXY) .or. &
		.not. allocated(NoahmpIO%RSSUNXY) .or. .not. allocated(NoahmpIO%RSSHAXY) .or. &
		.not. allocated(NoahmpIO%BGAPXY) .or. .not. allocated(NoahmpIO%WGAPXY) .or. &
		.not. allocated(NoahmpIO%TGVXY) .or. .not. allocated(NoahmpIO%TGBXY) .or. &
		.not. allocated(NoahmpIO%CHVXY) .or. .not. allocated(NoahmpIO%CHBXY) .or. &
		.not. allocated(NoahmpIO%CHLEAFXY) .or. .not. allocated(NoahmpIO%CHUCXY) .or. &
		.not. allocated(NoahmpIO%CHV2XY) .or. .not. allocated(NoahmpIO%CHB2XY) .or. &
		.not. allocated(NoahmpIO%ALBEDO) .or. .not. allocated(NoahmpIO%RS) .or. &
		.not. allocated(NoahmpIO%SOILENERGY) .or. .not. allocated(NoahmpIO%SNOWENERGY) .or. &
		.not. allocated(NoahmpIO%LH) .or. .not. allocated(NoahmpIO%HFX) .or. &
		.not. allocated(NoahmpIO%GRDFLX) .or. .not. allocated(NoahmpIO%FSAXY) .or. &
		.not. allocated(NoahmpIO%FIRAXY) .or. .not. allocated(NoahmpIO%APARXY) .or. &
		.not. allocated(NoahmpIO%SAVXY) .or. .not. allocated(NoahmpIO%SAGXY) .or. &
		.not. allocated(NoahmpIO%IRCXY) .or. .not. allocated(NoahmpIO%IRGXY) .or. &
		.not. allocated(NoahmpIO%SHCXY) .or. .not. allocated(NoahmpIO%SHGXY) .or. &
		.not. allocated(NoahmpIO%EVGXY) .or. .not. allocated(NoahmpIO%GHVXY) .or. &
		.not. allocated(NoahmpIO%IRBXY) .or. .not. allocated(NoahmpIO%SHBXY) .or. &
		.not. allocated(NoahmpIO%EVBXY) .or. .not. allocated(NoahmpIO%GHBXY) .or. &
		.not. allocated(NoahmpIO%TRXY) .or. .not. allocated(NoahmpIO%EVCXY) .or. &
		.not. allocated(NoahmpIO%CANHSXY) .or. .not. allocated(NoahmpIO%PAHXY) .or. &
		.not. allocated(NoahmpIO%PAHGXY) .or. .not. allocated(NoahmpIO%PAHVXY) .or. &
		.not. allocated(NoahmpIO%PAHBXY) .or. .not. allocated(NoahmpIO%EFLXBXY) .or. &
		.not. allocated(NoahmpIO%QFX) .or. .not. allocated(NoahmpIO%SMSTAV) .or. &
		.not. allocated(NoahmpIO%SMSTOT) .or. .not. allocated(NoahmpIO%SNOWC) .or. &
		.not. allocated(NoahmpIO%CANWAT) .or. .not. allocated(NoahmpIO%QSNOWXY) .or. &
		.not. allocated(NoahmpIO%QRAINXY) .or. .not. allocated(NoahmpIO%RUNSFXY) .or. &
		.not. allocated(NoahmpIO%RUNSBXY) .or. .not. allocated(NoahmpIO%ECANXY) .or. &
		.not. allocated(NoahmpIO%EDIRXY) .or. .not. allocated(NoahmpIO%ETRANXY) .or. &
		.not. allocated(NoahmpIO%QINTSXY) .or. .not. allocated(NoahmpIO%QINTRXY) .or. &
		.not. allocated(NoahmpIO%QDRIPSXY) .or. .not. allocated(NoahmpIO%QDRIPRXY) .or. &
		.not. allocated(NoahmpIO%QTHROSXY) .or. .not. allocated(NoahmpIO%QTHRORXY) .or. &
		.not. allocated(NoahmpIO%QSNSUBXY) .or. .not. allocated(NoahmpIO%QSNFROXY) .or. &
		.not. allocated(NoahmpIO%QSUBCXY) .or. .not. allocated(NoahmpIO%QFROCXY) .or. &
		.not. allocated(NoahmpIO%QEVACXY) .or. .not. allocated(NoahmpIO%QDEWCXY) .or. &
		.not. allocated(NoahmpIO%QFRZCXY) .or. .not. allocated(NoahmpIO%QMELTCXY) .or. &
		.not. allocated(NoahmpIO%QSNBOTXY) .or. .not. allocated(NoahmpIO%QMELTXY) .or. &
		.not. allocated(NoahmpIO%PONDINGXY) .or. .not. allocated(NoahmpIO%FPICEXY) .or. &
		.not. allocated(NoahmpIO%RAINLSM) .or. .not. allocated(NoahmpIO%SNOWLSM) .or. &
		.not. allocated(NoahmpIO%ALBSFCDIRXY) .or. .not. allocated(NoahmpIO%ALBSFCDIFXY) .or. &
		.not. allocated(NoahmpIO%ALBSNOWDIRXY) .or. .not. allocated(NoahmpIO%ALBSNOWDIFXY)) then
	   errflg = 1
	   errmsg = 'NoahmpModular_run: internal diagnostic output arrays are not allocated'
	   return
	end if

    if ( &
        size(NoahmpIO%FORCTLSM) < horizontal_loop_extent .or. &
        size(NoahmpIO%FORCQLSM) < horizontal_loop_extent .or. &
        size(NoahmpIO%FORCPLSM) < horizontal_loop_extent .or. &
        size(NoahmpIO%FORCWLSM) < horizontal_loop_extent .or. &
        size(NoahmpIO%TSK) < horizontal_loop_extent .or. &
        size(NoahmpIO%EMISS) < horizontal_loop_extent .or. &
        size(NoahmpIO%Z0) < horizontal_loop_extent .or. &
        size(NoahmpIO%T2MVXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%T2MBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%Q2MVXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%Q2MBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%FVEGXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%RSSUNXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%RSSHAXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%BGAPXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%WGAPXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%TGVXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%TGBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CHVXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CHBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CHLEAFXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CHUCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CHV2XY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CHB2XY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ALBEDO) < horizontal_loop_extent .or. &
        size(NoahmpIO%RS) < horizontal_loop_extent .or. &
        size(NoahmpIO%SOILENERGY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SNOWENERGY) < horizontal_loop_extent .or. &
        size(NoahmpIO%LH) < horizontal_loop_extent .or. &
        size(NoahmpIO%HFX) < horizontal_loop_extent .or. &
        size(NoahmpIO%GRDFLX) < horizontal_loop_extent .or. &
        size(NoahmpIO%FSAXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%FIRAXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%APARXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SAVXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SAGXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRGXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SHCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SHGXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%EVGXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%GHVXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%IRBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%SHBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%EVBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%GHBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%TRXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%EVCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%CANHSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%PAHXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%PAHGXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%PAHVXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%PAHBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%EFLXBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QFX) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMSTAV) < horizontal_loop_extent .or. &
        size(NoahmpIO%SMSTOT) < horizontal_loop_extent .or. &
        size(NoahmpIO%SNOWC) < horizontal_loop_extent .or. &
        size(NoahmpIO%CANWAT) < horizontal_loop_extent .or. &
        size(NoahmpIO%QSNOWXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QRAINXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%RUNSFXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%RUNSBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ECANXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%EDIRXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%ETRANXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QINTSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QINTRXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QDRIPSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QDRIPRXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QTHROSXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QTHRORXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QSNSUBXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QSNFROXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QSUBCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QFROCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QEVACXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QDEWCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QFRZCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QMELTCXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QSNBOTXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%QMELTXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%PONDINGXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%FPICEXY) < horizontal_loop_extent .or. &
        size(NoahmpIO%RAINLSM) < horizontal_loop_extent .or. &
        size(NoahmpIO%SNOWLSM) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal scalar diagnostic arrays are too small'
       return
    end if

    if (size(NoahmpIO%ALBSFCDIRXY, 1) < number_of_shortwave_radiation_bands .or. &
        size(NoahmpIO%ALBSFCDIRXY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%ALBSFCDIFXY, 1) < number_of_shortwave_radiation_bands .or. &
        size(NoahmpIO%ALBSFCDIFXY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%ALBSNOWDIRXY, 1) < number_of_shortwave_radiation_bands .or. &
        size(NoahmpIO%ALBSNOWDIRXY, 2) < horizontal_loop_extent .or. &
        size(NoahmpIO%ALBSNOWDIFXY, 1) < number_of_shortwave_radiation_bands .or. &
        size(NoahmpIO%ALBSNOWDIFXY, 2) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: internal banded diagnostic arrays are too small'
       return
    end if

    NoahmpIO%DTBL = real(timestep_for_physics, kind=kind_phys)

    NoahmpIO%calculate_soil = .true.
    NoahmpIO%soil_update_steps = 1
    NoahmpIO%JULIAN = forecast_julian_day
    NoahmpIO%YEARLEN = number_of_days_in_current_year
    NoahmpIO%DX(1:number_of_active_columns) = &
         characteristic_grid_lengthscale(active_source_column(1:number_of_active_columns))
    NoahmpIO%DY(1:number_of_active_columns) = NoahmpIO%DX(1:number_of_active_columns)
    NoahmpIO%NSOIL = vertical_dimension_of_soil
    NoahmpIO%NSNOW = vertical_dimension_of_surface_snow
    NoahmpIO%NUMRAD = number_of_shortwave_radiation_bands
    NoahmpIO%ITIMESTEP = timestep_index
    NoahmpIO%THSFC_LOC = flag_for_reference_pressure_theta

    NoahmpIO%SLOPETYP(1:number_of_active_columns) = &
         surface_slope_classification(active_source_column(1:number_of_active_columns))

    NoahmpIO%ZSOIL(1:vertical_dimension_of_soil) = &
         depth_of_soil_layer_interfaces(1:vertical_dimension_of_soil)
    NoahmpIO%IVGTYP(1:number_of_active_columns) = vegetation_category(active_source_column(1:number_of_active_columns))
    NoahmpIO%ISLTYP(1:number_of_active_columns) = soil_category(active_source_column(1:number_of_active_columns))
    NoahmpIO%SOILCOL(1:number_of_active_columns) = &
         soil_color_classification(active_source_column(1:number_of_active_columns))

    if (NoahmpIO%IOPT_CROP > 0 .or. NoahmpIO%IOPT_IRR == 2) then
       NoahmpIO%CROPCAT(1:number_of_active_columns) = &
            merge(NoahmpIO%DEFAULT_CROP_TABLE, 0, &
                  NoahmpIO%IVGTYP(1:number_of_active_columns) == NoahmpIO%ISCROP_TABLE)
    else
       NoahmpIO%CROPCAT(1:number_of_active_columns) = 0
    end if
    NoahmpIO%XICE(1:number_of_active_columns) = merge(1.0_kind_phys, 0.0_kind_phys, &
         NoahmpIO%IVGTYP(1:number_of_active_columns) == NoahmpIO%ISICE_TABLE)
    NoahmpIO%ICE(1:number_of_active_columns) = merge(-1, 0, &
         NoahmpIO%IVGTYP(1:number_of_active_columns) == NoahmpIO%ISICE_TABLE)
    NoahmpIO%SOILCL1(1:number_of_active_columns) = NoahmpIO%ISLTYP(1:number_of_active_columns)
    NoahmpIO%SOILCL2(1:number_of_active_columns) = NoahmpIO%ISLTYP(1:number_of_active_columns)
    NoahmpIO%SOILCL3(1:number_of_active_columns) = NoahmpIO%ISLTYP(1:number_of_active_columns)
    NoahmpIO%SOILCL4(1:number_of_active_columns) = NoahmpIO%ISLTYP(1:number_of_active_columns)
    NoahmpIO%source_i(1:number_of_active_columns) = active_source_column(1:number_of_active_columns)
    NoahmpIO%source_j(1:number_of_active_columns) = -9999
    NoahmpIO%ISNOWXY(1:number_of_active_columns) = nint( &
         active_snow_layer_lower_index(active_source_column(1:number_of_active_columns)))
    NoahmpIO%ZSNSOXY(lower_bound_of_vertical_dimension_of_surface_snow:vertical_dimension_of_soil, &
                  1:number_of_active_columns) = &
         transpose(depth_of_snow_soil_layer_interfaces( &
              active_source_column(1:number_of_active_columns), &
              lower_bound_of_vertical_dimension_of_surface_snow:vertical_dimension_of_soil))


    NoahmpIO%LAI(1:number_of_active_columns) = &
         leaf_area_index(active_source_column(1:number_of_active_columns))
    NoahmpIO%XSAIXY(1:number_of_active_columns) = &
         stem_area_index(active_source_column(1:number_of_active_columns))
    NoahmpIO%TGXY(1:number_of_active_columns) = &
         ground_temperature(active_source_column(1:number_of_active_columns))
    NoahmpIO%TVXY(1:number_of_active_columns) = &
         canopy_temperature(active_source_column(1:number_of_active_columns))
    NoahmpIO%TAUSSXY(1:number_of_active_columns) = &
         dimensionless_snow_age(active_source_column(1:number_of_active_columns))
    NoahmpIO%ALBOLDXY(1:number_of_active_columns) = &
         snow_albedo_on_previous_timestep(active_source_column(1:number_of_active_columns))
    NoahmpIO%EAHXY(1:number_of_active_columns) = &
         canopy_air_vapor_pressure(active_source_column(1:number_of_active_columns))
    NoahmpIO%TAHXY(1:number_of_active_columns) = &
         canopy_air_temperature(active_source_column(1:number_of_active_columns))
    NoahmpIO%CHXY(1:number_of_active_columns) = &
         surface_exchange_coefficient_for_heat(active_source_column(1:number_of_active_columns))
    NoahmpIO%CMXY(1:number_of_active_columns) = &
         surface_exchange_coefficient_for_momentum(active_source_column(1:number_of_active_columns))
    NoahmpIO%GVFMAX(1:number_of_active_columns) = &
         maximum_vegetation_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%VEGFRA(1:number_of_active_columns) = &
         vegetation_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%TSLB(1:vertical_dimension_of_soil,1:number_of_active_columns) = &
         transpose(soil_temperature(active_source_column(1:number_of_active_columns),1:vertical_dimension_of_soil))
    NoahmpIO%TSNOXY(lower_bound_of_vertical_dimension_of_surface_snow: &
                 upper_bound_of_vertical_dimension_of_surface_snow, &
                 1:number_of_active_columns) = &
         transpose(snow_temperature(active_source_column(1:number_of_active_columns), &
                          lower_bound_of_vertical_dimension_of_surface_snow: &
                          upper_bound_of_vertical_dimension_of_surface_snow))
    NoahmpIO%ALBSOILDIRXY(1:number_of_shortwave_radiation_bands, &
                              1:number_of_active_columns) = &
         transpose(direct_soil_albedo(active_source_column(1:number_of_active_columns), &
                                      1:number_of_shortwave_radiation_bands))
    NoahmpIO%ALBSOILDIFXY(1:number_of_shortwave_radiation_bands, &
                              1:number_of_active_columns) = &
         transpose(diffuse_soil_albedo(active_source_column(1:number_of_active_columns), &
                                       1:number_of_shortwave_radiation_bands))

    NoahmpIO%CANLIQXY(1:number_of_active_columns) = &
         canopy_liquid_water(active_source_column(1:number_of_active_columns))
    NoahmpIO%CANICEXY(1:number_of_active_columns) = &
         canopy_ice_water(active_source_column(1:number_of_active_columns))
    NoahmpIO%FWETXY(1:number_of_active_columns) = &
         canopy_wet_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%SNOW(1:number_of_active_columns) = &
         snow_water_equivalent(active_source_column(1:number_of_active_columns))
    NoahmpIO%SNEQVOXY(1:number_of_active_columns) = &
         snow_water_equivalent_on_previous_timestep(active_source_column(1:number_of_active_columns))
    NoahmpIO%SNOWH(1:number_of_active_columns) = 0.001_kind_phys * &
         snow_depth(active_source_column(1:number_of_active_columns))
    NoahmpIO%FIFRACT(1:number_of_active_columns) = &
         flood_irrigation_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRWATFI(1:number_of_active_columns) = &
         flood_irrigation_water_amount(active_source_column(1:number_of_active_columns))
    NoahmpIO%MIFRACT(1:number_of_active_columns) = &
         micro_irrigation_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRWATMI(1:number_of_active_columns) = &
         micro_irrigation_water_amount(active_source_column(1:number_of_active_columns))
    NoahmpIO%SIFRACT(1:number_of_active_columns) = &
         sprinkler_irrigation_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRWATSI(1:number_of_active_columns) = &
         sprinkler_irrigation_water_amount(active_source_column(1:number_of_active_columns))
    NoahmpIO%ZWTXY(1:number_of_active_columns) = &
         water_table_depth(active_source_column(1:number_of_active_columns))
    NoahmpIO%SMCWTDXY(1:number_of_active_columns) = &
         soil_moisture_below_soil_column(active_source_column(1:number_of_active_columns))
    NoahmpIO%TD_FRACTION(1:number_of_active_columns) = &
         tile_drainage_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%WAXY(1:number_of_active_columns) = &
         aquifer_water_storage(active_source_column(1:number_of_active_columns))
    NoahmpIO%WTXY(1:number_of_active_columns) = &
         soil_aquifer_water_storage(active_source_column(1:number_of_active_columns))
    NoahmpIO%WSLAKEXY(1:number_of_active_columns) = &
         lake_water_storage(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRFRACT(1:number_of_active_columns) = &
         grid_irrigation_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRNUMSI(1:number_of_active_columns) = &
         sprinkler_irrigation_event_count(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRNUMMI(1:number_of_active_columns) = &
         micro_irrigation_event_count(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRNUMFI(1:number_of_active_columns) = &
         flood_irrigation_event_count(active_source_column(1:number_of_active_columns))
    NoahmpIO%FSATXY(1:number_of_active_columns) = &
         soil_saturated_fraction(active_source_column(1:number_of_active_columns))
    NoahmpIO%WSURFXY(1:number_of_active_columns) = &
         wetland_water_storage(active_source_column(1:number_of_active_columns))
    NoahmpIO%SNICEXY(lower_bound_of_vertical_dimension_of_surface_snow: &
                  upper_bound_of_vertical_dimension_of_surface_snow, &
                  1:number_of_active_columns) = &
         transpose(snow_layer_ice(active_source_column(1:number_of_active_columns), &
                        lower_bound_of_vertical_dimension_of_surface_snow: &
                        upper_bound_of_vertical_dimension_of_surface_snow))
    NoahmpIO%SNLIQXY(lower_bound_of_vertical_dimension_of_surface_snow: &
                  upper_bound_of_vertical_dimension_of_surface_snow, &
                  1:number_of_active_columns) = &
         transpose(snow_layer_liquid_water(active_source_column(1:number_of_active_columns), &
                                 lower_bound_of_vertical_dimension_of_surface_snow: &
                                 upper_bound_of_vertical_dimension_of_surface_snow))
    NoahmpIO%SH2O(1:vertical_dimension_of_soil,1:number_of_active_columns) = &
         transpose(soil_liquid_water(active_source_column(1:number_of_active_columns),1:vertical_dimension_of_soil))
    NoahmpIO%SMOIS(1:vertical_dimension_of_soil,1:number_of_active_columns) = &
         transpose(soil_moisture(active_source_column(1:number_of_active_columns),1:vertical_dimension_of_soil))
    NoahmpIO%SMOISEQ(1:vertical_dimension_of_soil,1:number_of_active_columns) = &
         transpose(equilibrium_soil_moisture(active_source_column(1:number_of_active_columns), &
                                   1:vertical_dimension_of_soil))
    NoahmpIO%ACC_SSOILXY(1:number_of_active_columns) = &
         accumulated_ground_heat_flux(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRRSPLH(1:number_of_active_columns) = &
         sprinkler_heat_accumulation(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_QSEVAXY(1:number_of_active_columns) = &
         accumulated_soil_surface_evaporation(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_QINSURXY(1:number_of_active_columns) = &
         accumulated_soil_surface_inflow(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_DWATERXY(1:number_of_active_columns) = &
         accumulated_surface_water_change(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_PRCPXY(1:number_of_active_columns) = &
         accumulated_precipitation(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_ECANXY(1:number_of_active_columns) = &
         accumulated_canopy_evaporation(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_ETRANXY(1:number_of_active_columns) = &
         accumulated_transpiration(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_EDIRXY(1:number_of_active_columns) = &
         accumulated_ground_evaporation(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_GLAFLWXY(1:number_of_active_columns) = &
         accumulated_glacier_excess_flow(active_source_column(1:number_of_active_columns))
    NoahmpIO%SFCRUNOFF(1:number_of_active_columns) = &
         accumulated_surface_runoff(active_source_column(1:number_of_active_columns))
    NoahmpIO%UDRUNOFF(1:number_of_active_columns) = &
         accumulated_subsurface_runoff(active_source_column(1:number_of_active_columns))
    NoahmpIO%QTDRAIN(1:number_of_active_columns) = &
         accumulated_tile_drainage(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACSNOW(1:number_of_active_columns) = &
         accumulated_snowfall(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACSNOM(1:number_of_active_columns) = &
         accumulated_snowmelt(active_source_column(1:number_of_active_columns))
    NoahmpIO%RECHXY(1:number_of_active_columns) = &
         accumulated_shallow_groundwater_recharge(active_source_column(1:number_of_active_columns))
    NoahmpIO%DEEPRECHXY(1:number_of_active_columns) = &
         accumulated_deep_groundwater_recharge(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRSIVOL(1:number_of_active_columns) = &
         accumulated_sprinkler_irrigation(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRMIVOL(1:number_of_active_columns) = &
         accumulated_micro_irrigation(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRFIVOL(1:number_of_active_columns) = &
         accumulated_flood_irrigation(active_source_column(1:number_of_active_columns))
    NoahmpIO%IRELOSS(1:number_of_active_columns) = &
         accumulated_sprinkler_evaporation_loss(active_source_column(1:number_of_active_columns))
    NoahmpIO%ACC_ETRANIXY(1:vertical_dimension_of_soil,1:number_of_active_columns) = &
         transpose(accumulated_soil_layer_transpiration( &
              active_source_column(1:number_of_active_columns),1:vertical_dimension_of_soil))

    ! Preserve the existing per-soil-cycle reset semantics.
    NoahmpIO%ACC_SSOILXY(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%ACC_QINSURXY(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%ACC_QSEVAXY(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%ACC_ETRANIXY(1:vertical_dimension_of_soil,1:number_of_active_columns) = &
         0.0_kind_phys
    NoahmpIO%ACC_DWATERXY(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%ACC_PRCPXY(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%ACC_ECANXY(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%ACC_ETRANXY(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%ACC_EDIRXY(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%ACC_GLAFLWXY(1:number_of_active_columns) = 0.0_kind_phys

    NoahmpIO%PGSXY(1:number_of_active_columns) = &
         plant_growth_stage(active_source_column(1:number_of_active_columns))
    NoahmpIO%GRAINXY(1:number_of_active_columns) = &
         grain_mass(active_source_column(1:number_of_active_columns))
    NoahmpIO%GDDXY(1:number_of_active_columns) = &
         growing_degree_day(active_source_column(1:number_of_active_columns))
    NoahmpIO%LFMASSXY(1:number_of_active_columns) = &
         leaf_mass(active_source_column(1:number_of_active_columns))
    NoahmpIO%RTMASSXY(1:number_of_active_columns) = &
         root_mass(active_source_column(1:number_of_active_columns))
    NoahmpIO%STMASSXY(1:number_of_active_columns) = &
         stem_mass(active_source_column(1:number_of_active_columns))
    NoahmpIO%WOODXY(1:number_of_active_columns) = &
         wood_mass(active_source_column(1:number_of_active_columns))
    NoahmpIO%STBLCPXY(1:number_of_active_columns) = &
         deep_soil_carbon_mass(active_source_column(1:number_of_active_columns))
    NoahmpIO%FASTCPXY(1:number_of_active_columns) = &
         shallow_soil_carbon_mass(active_source_column(1:number_of_active_columns))
    ! Zero-valued spatial overrides retain the canonical crop-table calendar
    ! and growing-degree-day thresholds selected by CROPCAT.
    NoahmpIO%PLANTING(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%HARVEST(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%SEASON_GDD(1:number_of_active_columns) = 0.0_kind_phys

    if (any(NoahmpIO%ISLTYP(1:number_of_active_columns) < 1) .or. &
        any(NoahmpIO%ISLTYP(1:number_of_active_columns) > NoahmpIO%SLCATS_TABLE)) then
       errflg = 1
       errmsg = 'NoahmpModular_run: soil classification is outside the modular table'
       return
    end if
    if (size(snicar_hydrophobic_black_carbon_deposition_flux) < horizontal_loop_extent .or. &
        size(snicar_hydrophilic_black_carbon_deposition_flux) < horizontal_loop_extent .or. &
        size(snicar_hydrophobic_organic_carbon_deposition_flux) < horizontal_loop_extent .or. &
        size(snicar_hydrophilic_organic_carbon_deposition_flux) < horizontal_loop_extent .or. &
        size(snicar_dust_species_1_deposition_flux) < horizontal_loop_extent .or. &
        size(snicar_dust_species_2_deposition_flux) < horizontal_loop_extent .or. &
        size(snicar_dust_species_3_deposition_flux) < horizontal_loop_extent .or. &
        size(snicar_dust_species_4_deposition_flux) < horizontal_loop_extent .or. &
        size(snicar_dust_species_5_deposition_flux) < horizontal_loop_extent) then
       errflg = 1
       errmsg = 'NoahmpModular_run: SNICAR aerosol deposition forcing is too small'
       return
    end if
    do n = 1, number_of_active_columns
       NoahmpIO%quartz_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%QUARTZ_TABLE(NoahmpIO%ISLTYP(n))
       NoahmpIO%BEXP_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%BEXP_TABLE(NoahmpIO%ISLTYP(n))
       NoahmpIO%SMCDRY_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%SMCDRY_TABLE(NoahmpIO%ISLTYP(n))
       NoahmpIO%SMCWLT_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%SMCWLT_TABLE(NoahmpIO%ISLTYP(n))
       NoahmpIO%SMCREF_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%SMCREF_TABLE(NoahmpIO%ISLTYP(n))
       NoahmpIO%SMCMAX_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%SMCMAX_TABLE(NoahmpIO%ISLTYP(n))
       NoahmpIO%DKSAT_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%DKSAT_TABLE(NoahmpIO%ISLTYP(n))
       NoahmpIO%DWSAT_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%DWSAT_TABLE(NoahmpIO%ISLTYP(n))
       NoahmpIO%PSISAT_3D(1:vertical_dimension_of_soil,n) = &
            NoahmpIO%PSISAT_TABLE(NoahmpIO%ISLTYP(n))
    end do
    NoahmpIO%REFDK_2D(1:number_of_active_columns) = NoahmpIO%REFDK_TABLE
    NoahmpIO%REFKDT_2D(1:number_of_active_columns) = NoahmpIO%REFKDT_TABLE
    NoahmpIO%BVIC_2D(1:number_of_active_columns) = &
         NoahmpIO%BVIC_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%AXAJ_2D(1:number_of_active_columns) = &
         NoahmpIO%AXAJ_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%BXAJ_2D(1:number_of_active_columns) = &
         NoahmpIO%BXAJ_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%XXAJ_2D(1:number_of_active_columns) = &
         NoahmpIO%XXAJ_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%BDVIC_2D(1:number_of_active_columns) = &
         NoahmpIO%BDVIC_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%GDVIC_2D(1:number_of_active_columns) = &
         NoahmpIO%GDVIC_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%BBVIC_2D(1:number_of_active_columns) = &
         NoahmpIO%BBVIC_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%IRR_FRAC_2D(1:number_of_active_columns) = NoahmpIO%IRR_FRAC_TABLE
    NoahmpIO%IRR_HAR_2D(1:number_of_active_columns) = &
         real(NoahmpIO%IRR_HAR_TABLE, kind=kind_phys)
    NoahmpIO%IRR_LAI_2D(1:number_of_active_columns) = NoahmpIO%IRR_LAI_TABLE
    NoahmpIO%IRR_MAD_2D(1:number_of_active_columns) = NoahmpIO%IRR_MAD_TABLE
    NoahmpIO%FILOSS_2D(1:number_of_active_columns) = NoahmpIO%FILOSS_TABLE
    NoahmpIO%SPRIR_RATE_2D(1:number_of_active_columns) = NoahmpIO%SPRIR_RATE_TABLE
    NoahmpIO%MICIR_RATE_2D(1:number_of_active_columns) = NoahmpIO%MICIR_RATE_TABLE
    NoahmpIO%FIRTFAC_2D(1:number_of_active_columns) = NoahmpIO%FIRTFAC_TABLE
    NoahmpIO%IR_RAIN_2D(1:number_of_active_columns) = NoahmpIO%IR_RAIN_TABLE
    NoahmpIO%KLAT_FAC(1:number_of_active_columns) = &
         NoahmpIO%KLAT_FAC_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%TDSMC_FAC(1:number_of_active_columns) = &
         NoahmpIO%TDSMC_FAC_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%TD_DC(1:number_of_active_columns) = &
         NoahmpIO%TD_DC_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%TD_DCOEF(1:number_of_active_columns) = &
         NoahmpIO%TD_DCOEF_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%TD_DDRAIN(1:number_of_active_columns) = &
         NoahmpIO%TD_DDRAIN_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%TD_RADI(1:number_of_active_columns) = &
         NoahmpIO%TD_RADI_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%TD_SPAC(1:number_of_active_columns) = &
         NoahmpIO%TD_SPAC_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    NoahmpIO%FSATMX(1:number_of_active_columns) = NoahmpIO%FSATMX_TABLE
    NoahmpIO%WCAP(1:number_of_active_columns) = NoahmpIO%WCAP_TABLE

    NoahmpIO%XLAT(1:number_of_active_columns) = latitude(active_source_column(1:number_of_active_columns))
    NoahmpIO%COSZEN(1:number_of_active_columns) = &
         instantaneous_cosine_of_zenith_angle(active_source_column(1:number_of_active_columns))
    ! Noah-MP does not use this quantity.  Match the legacy CCPP driver
    ! without adding an unused field to the host-model interface.
    NoahmpIO%DZ8W(NoahmpIO%KTS, 1:number_of_active_columns) = -9999.0_kind_phys
    NoahmpIO%FORCZLSM(1:number_of_active_columns) = &
         atmospheric_forcing_height(active_source_column(1:number_of_active_columns))

    NoahmpIO%SWDOWN(1:number_of_active_columns) = &
         surface_downwelling_shortwave_flux(active_source_column(1:number_of_active_columns))
    NoahmpIO%GLW(1:number_of_active_columns) = &
         surface_downwelling_longwave_flux(active_source_column(1:number_of_active_columns))
    if (NoahmpIO%IOPT_ALB == 3) then
       NoahmpIO%DepBChydrophoXY(1:number_of_active_columns) = &
         snicar_hydrophobic_black_carbon_deposition_flux( &
              active_source_column(1:number_of_active_columns))
       NoahmpIO%DepBChydrophiXY(1:number_of_active_columns) = &
         snicar_hydrophilic_black_carbon_deposition_flux( &
              active_source_column(1:number_of_active_columns))
       NoahmpIO%DepOChydrophoXY(1:number_of_active_columns) = &
         snicar_hydrophobic_organic_carbon_deposition_flux( &
              active_source_column(1:number_of_active_columns))
       NoahmpIO%DepOChydrophiXY(1:number_of_active_columns) = &
         snicar_hydrophilic_organic_carbon_deposition_flux( &
              active_source_column(1:number_of_active_columns))
       NoahmpIO%DepDust1XY(1:number_of_active_columns) = &
         snicar_dust_species_1_deposition_flux( &
              active_source_column(1:number_of_active_columns))
       NoahmpIO%DepDust2XY(1:number_of_active_columns) = &
         snicar_dust_species_2_deposition_flux( &
              active_source_column(1:number_of_active_columns))
       NoahmpIO%DepDust3XY(1:number_of_active_columns) = &
         snicar_dust_species_3_deposition_flux( &
              active_source_column(1:number_of_active_columns))
       NoahmpIO%DepDust4XY(1:number_of_active_columns) = &
         snicar_dust_species_4_deposition_flux( &
              active_source_column(1:number_of_active_columns))
       NoahmpIO%DepDust5XY(1:number_of_active_columns) = &
         snicar_dust_species_5_deposition_flux( &
              active_source_column(1:number_of_active_columns))
    endif
    shortwave_component_total(1:number_of_active_columns) = &
         direct_nir_shortwave_flux_at_surface( &
              active_source_column(1:number_of_active_columns)) + &
         diffuse_nir_shortwave_flux_at_surface( &
              active_source_column(1:number_of_active_columns)) + &
         direct_vis_shortwave_flux_at_surface( &
              active_source_column(1:number_of_active_columns)) + &
         diffuse_vis_shortwave_flux_at_surface( &
              active_source_column(1:number_of_active_columns))
    ! At night these fractions do not affect the zero shortwave flux. Retain
    ! the modular driver defaults so every packed value remains initialized.
    NoahmpIO%RadSwDirFrac(1:number_of_active_columns) = 0.7_kind_phys
    NoahmpIO%RadSwVisFrac(1:number_of_active_columns) = 0.5_kind_phys
	
    where (shortwave_component_total(1:number_of_active_columns) > tiny(1.0_kind_phys))
       NoahmpIO%RadSwDirFrac(1:number_of_active_columns) = &
            (direct_nir_shortwave_flux_at_surface( &
                 active_source_column(1:number_of_active_columns)) + &
             direct_vis_shortwave_flux_at_surface( &
                 active_source_column(1:number_of_active_columns))) / &
            shortwave_component_total(1:number_of_active_columns)
       NoahmpIO%RadSwVisFrac(1:number_of_active_columns) = &
            (direct_vis_shortwave_flux_at_surface( &
                 active_source_column(1:number_of_active_columns)) + &
             diffuse_vis_shortwave_flux_at_surface( &
                 active_source_column(1:number_of_active_columns))) / &
            shortwave_component_total(1:number_of_active_columns)
    end where
	
    NoahmpIO%PS(1:number_of_active_columns) = &
         surface_air_pressure(active_source_column(1:number_of_active_columns))
    NoahmpIO%PRSL1(1:number_of_active_columns) = &
         air_pressure_at_surface_adjacent_layer(active_source_column(1:number_of_active_columns))
    NoahmpIO%PBLH(1:number_of_active_columns) = &
         atmosphere_boundary_layer_thickness(active_source_column(1:number_of_active_columns))
    NoahmpIO%TMN(1:number_of_active_columns) = &
         deep_soil_temperature(active_source_column(1:number_of_active_columns))
    NoahmpIO%T_PHY(NoahmpIO%KTS, 1:number_of_active_columns) = &
         air_temperature_at_lowest_model_layer(active_source_column(1:number_of_active_columns))
    NoahmpIO%QV_CURR(NoahmpIO%KTS, 1:number_of_active_columns) = &
         specific_humidity_at_lowest_model_layer(active_source_column(1:number_of_active_columns))
    NoahmpIO%U_PHY(NoahmpIO%KTS, 1:number_of_active_columns) = &
         eastward_wind_at_lowest_model_layer(active_source_column(1:number_of_active_columns))
    NoahmpIO%V_PHY(NoahmpIO%KTS, 1:number_of_active_columns) = &
         northward_wind_at_lowest_model_layer(active_source_column(1:number_of_active_columns))
    NoahmpIO%MP_RAINC(1:number_of_active_columns) = &
         convective_precipitation_rate_on_previous_timestep(active_source_column(1:number_of_active_columns))
    NoahmpIO%MP_RAINNC(1:number_of_active_columns) = &
         explicit_precipitation_rate_on_previous_timestep(active_source_column(1:number_of_active_columns))
    NoahmpIO%MP_SHCV(1:number_of_active_columns) = 0.0_kind_phys
    NoahmpIO%RAINBL(1:number_of_active_columns) = 1000.0_kind_phys * &
         precipitation_amount_on_dynamics_timestep( &
              active_source_column(1:number_of_active_columns)) / timestep_for_physics
    NoahmpIO%MP_SNOW(1:number_of_active_columns) = &
         snowfall_rate_on_previous_timestep(active_source_column(1:number_of_active_columns))
    NoahmpIO%MP_GRAUP(1:number_of_active_columns) = &
         graupel_precipitation_rate_on_previous_timestep(active_source_column(1:number_of_active_columns))
    NoahmpIO%MP_HAIL(1:number_of_active_columns) = &
         ice_precipitation_rate_on_previous_timestep(active_source_column(1:number_of_active_columns))
    NoahmpIO%SR(1:number_of_active_columns) = &
         precipitation_type(active_source_column(1:number_of_active_columns))

    call NoahmpDriverMain(NoahmpIO, &
                        packed_apply_urban_irrigation(1:number_of_active_columns), &
                        errmsg, errflg)

    if (errflg /= 0) return

    air_density(1:number_of_active_columns) = &
         NoahmpIO%FORCPLSM(1:number_of_active_columns) / &
         (ConstGasDryAir * NoahmpIO%FORCTLSM(1:number_of_active_columns) * &
          (1.0_kind_phys + 0.61_kind_phys * NoahmpIO%FORCQLSM(1:number_of_active_columns)))
    runoff_accumulation_time(1:number_of_active_columns) = NoahmpIO%DTBL
    if (NoahmpIO%calculate_soil) then
       where (NoahmpIO%ICE(1:number_of_active_columns) == 0)
          runoff_accumulation_time(1:number_of_active_columns) = &
               NoahmpIO%DTBL * real(NoahmpIO%soil_update_steps, kind=kind_phys)
       end where
    end if


    active_snow_layer_lower_index(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ISNOWXY(1:number_of_active_columns)
    depth_of_snow_soil_layer_interfaces( &
         active_source_column(1:number_of_active_columns), &
         lower_bound_of_vertical_dimension_of_surface_snow:vertical_dimension_of_soil) = &
         transpose(NoahmpIO%ZSNSOXY( &
              lower_bound_of_vertical_dimension_of_surface_snow:vertical_dimension_of_soil, &
              1:number_of_active_columns))

    leaf_area_index(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%LAI(1:number_of_active_columns)
    stem_area_index(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%XSAIXY(1:number_of_active_columns)
    specific_humidity_at_surface(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%QSFC(1:number_of_active_columns)
    ground_temperature(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%TGXY(1:number_of_active_columns)
    canopy_temperature(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%TVXY(1:number_of_active_columns)
    dimensionless_snow_age(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%TAUSSXY(1:number_of_active_columns)
    snow_albedo_on_previous_timestep(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ALBOLDXY(1:number_of_active_columns)
    canopy_air_vapor_pressure(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%EAHXY(1:number_of_active_columns)
    canopy_air_temperature(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%TAHXY(1:number_of_active_columns)
    surface_exchange_coefficient_for_heat(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%CHXY(1:number_of_active_columns)
    surface_exchange_coefficient_for_momentum(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%CMXY(1:number_of_active_columns)
    maximum_vegetation_fraction(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%GVFMAX(1:number_of_active_columns)
    vegetation_fraction(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%VEGFRA(1:number_of_active_columns)
    soil_temperature(active_source_column(1:number_of_active_columns),1:vertical_dimension_of_soil) = &
         transpose(NoahmpIO%TSLB(1:vertical_dimension_of_soil,1:number_of_active_columns))
    snow_temperature(active_source_column(1:number_of_active_columns), &
                     lower_bound_of_vertical_dimension_of_surface_snow: &
                     upper_bound_of_vertical_dimension_of_surface_snow) = &
         transpose(NoahmpIO%TSNOXY(lower_bound_of_vertical_dimension_of_surface_snow: &
                      upper_bound_of_vertical_dimension_of_surface_snow, &
                      1:number_of_active_columns))
					  
    direct_soil_albedo(active_source_column(1:number_of_active_columns), &
                             1:number_of_shortwave_radiation_bands) = &
         transpose(NoahmpIO%ALBSOILDIRXY(1:number_of_shortwave_radiation_bands, &
                                         1:number_of_active_columns))
    diffuse_soil_albedo(active_source_column(1:number_of_active_columns), &
                        1:number_of_shortwave_radiation_bands) = &
         transpose(NoahmpIO%ALBSOILDIFXY(1:number_of_shortwave_radiation_bands, &
                                         1:number_of_active_columns))

    canopy_liquid_water(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%CANLIQXY(1:number_of_active_columns)
    canopy_ice_water(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%CANICEXY(1:number_of_active_columns)
    canopy_wet_fraction(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%FWETXY(1:number_of_active_columns)
    snow_water_equivalent(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%SNOW(1:number_of_active_columns)
    snow_water_equivalent_on_previous_timestep(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%SNEQVOXY(1:number_of_active_columns)
    snow_depth(active_source_column(1:number_of_active_columns)) = &
         1000.0_kind_phys * NoahmpIO%SNOWH(1:number_of_active_columns)
    flood_irrigation_water_amount(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRWATFI(1:number_of_active_columns)
    micro_irrigation_water_amount(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRWATMI(1:number_of_active_columns)
    sprinkler_irrigation_water_amount(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRWATSI(1:number_of_active_columns)
    water_table_depth(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ZWTXY(1:number_of_active_columns)
    soil_moisture_below_soil_column(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%SMCWTDXY(1:number_of_active_columns)
    aquifer_water_storage(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%WAXY(1:number_of_active_columns)
    soil_aquifer_water_storage(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%WTXY(1:number_of_active_columns)
    lake_water_storage(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%WSLAKEXY(1:number_of_active_columns)
    sprinkler_irrigation_event_count(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRNUMSI(1:number_of_active_columns)
    micro_irrigation_event_count(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRNUMMI(1:number_of_active_columns)
    flood_irrigation_event_count(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRNUMFI(1:number_of_active_columns)
    soil_saturated_fraction(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%FSATXY(1:number_of_active_columns)
    wetland_water_storage(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%WSURFXY(1:number_of_active_columns)
    snow_layer_ice(active_source_column(1:number_of_active_columns), &
                   lower_bound_of_vertical_dimension_of_surface_snow: &
                   upper_bound_of_vertical_dimension_of_surface_snow) = &
         transpose(NoahmpIO%SNICEXY(lower_bound_of_vertical_dimension_of_surface_snow: &
                       upper_bound_of_vertical_dimension_of_surface_snow, &
                       1:number_of_active_columns))
    snow_layer_liquid_water(active_source_column(1:number_of_active_columns), &
                            lower_bound_of_vertical_dimension_of_surface_snow: &
                            upper_bound_of_vertical_dimension_of_surface_snow) = &
         transpose(NoahmpIO%SNLIQXY(lower_bound_of_vertical_dimension_of_surface_snow: &
                       upper_bound_of_vertical_dimension_of_surface_snow, &
                       1:number_of_active_columns))
    soil_liquid_water(active_source_column(1:number_of_active_columns),1:vertical_dimension_of_soil) = &
         transpose(NoahmpIO%SH2O(1:vertical_dimension_of_soil,1:number_of_active_columns))
    soil_moisture(active_source_column(1:number_of_active_columns),1:vertical_dimension_of_soil) = &
         transpose(NoahmpIO%SMOIS(1:vertical_dimension_of_soil,1:number_of_active_columns))
    accumulated_ground_heat_flux(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_SSOILXY(1:number_of_active_columns)
    sprinkler_heat_accumulation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRRSPLH(1:number_of_active_columns)
    accumulated_soil_surface_evaporation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_QSEVAXY(1:number_of_active_columns)
    accumulated_soil_surface_inflow(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_QINSURXY(1:number_of_active_columns)
    accumulated_surface_water_change(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_DWATERXY(1:number_of_active_columns)
    accumulated_precipitation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_PRCPXY(1:number_of_active_columns)
    accumulated_canopy_evaporation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_ECANXY(1:number_of_active_columns)
    accumulated_transpiration(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_ETRANXY(1:number_of_active_columns)
    accumulated_ground_evaporation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_EDIRXY(1:number_of_active_columns)
    accumulated_glacier_excess_flow(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACC_GLAFLWXY(1:number_of_active_columns)
    accumulated_surface_runoff(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%SFCRUNOFF(1:number_of_active_columns)
    accumulated_subsurface_runoff(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%UDRUNOFF(1:number_of_active_columns)
    accumulated_tile_drainage(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%QTDRAIN(1:number_of_active_columns)
    accumulated_snowfall(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACSNOW(1:number_of_active_columns)
    accumulated_snowmelt(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ACSNOM(1:number_of_active_columns)
    accumulated_shallow_groundwater_recharge(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%RECHXY(1:number_of_active_columns)
    accumulated_deep_groundwater_recharge(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%DEEPRECHXY(1:number_of_active_columns)
    accumulated_sprinkler_irrigation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRSIVOL(1:number_of_active_columns)
    accumulated_micro_irrigation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRMIVOL(1:number_of_active_columns)
    accumulated_flood_irrigation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRFIVOL(1:number_of_active_columns)
    accumulated_sprinkler_evaporation_loss(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%IRELOSS(1:number_of_active_columns)
    accumulated_soil_layer_transpiration( &
         active_source_column(1:number_of_active_columns),1:vertical_dimension_of_soil) = &
         transpose(NoahmpIO%ACC_ETRANIXY(1:vertical_dimension_of_soil, &
                                          1:number_of_active_columns))

    plant_growth_stage(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%PGSXY(1:number_of_active_columns)
    grain_mass(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%GRAINXY(1:number_of_active_columns)
    growing_degree_day(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%GDDXY(1:number_of_active_columns)
	
    leaf_mass(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%LFMASSXY(1:number_of_active_columns)
    root_mass(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%RTMASSXY(1:number_of_active_columns)
    stem_mass(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%STMASSXY(1:number_of_active_columns)
    wood_mass(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%WOODXY(1:number_of_active_columns)
    deep_soil_carbon_mass(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%STBLCPXY(1:number_of_active_columns)
    shallow_soil_carbon_mass(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%FASTCPXY(1:number_of_active_columns)

    surface_radiative_temperature(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%TSK(1:number_of_active_columns)
    surface_emissivity(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%EMISS(1:number_of_active_columns)
    surface_roughness_length(active_source_column(1:number_of_active_columns)) = &
         100.0_kind_phys * NoahmpIO%Z0(1:number_of_active_columns)
    temperature_at_2m_from_noahmp(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%T2MVXY(1:number_of_active_columns) * NoahmpIO%VEGFRA(1:number_of_active_columns) + &
         NoahmpIO%T2MBXY(1:number_of_active_columns) * &
         (1.0_kind_phys - NoahmpIO%VEGFRA(1:number_of_active_columns))
    specific_humidity_at_2m_from_noahmp(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%Q2MVXY(1:number_of_active_columns) / &
         (1.0_kind_phys + NoahmpIO%Q2MVXY(1:number_of_active_columns)) * &
         NoahmpIO%VEGFRA(1:number_of_active_columns) + &
         NoahmpIO%Q2MBXY(1:number_of_active_columns) / &
         (1.0_kind_phys + NoahmpIO%Q2MBXY(1:number_of_active_columns)) * &
         (1.0_kind_phys - NoahmpIO%VEGFRA(1:number_of_active_columns))
    surface_drag_wind_speed_for_momentum(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%CMXY(1:number_of_active_columns) * NoahmpIO%FORCWLSM(1:number_of_active_columns)
    surface_drag_mass_flux_for_heat_and_moisture( &
         active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%CHXY(1:number_of_active_columns) * air_density(1:number_of_active_columns)
		 
    where (NoahmpIO%ALBEDO(1:number_of_active_columns) > 0.0_kind_phys)
       surface_albedo(active_source_column(1:number_of_active_columns)) = &
            NoahmpIO%ALBEDO(1:number_of_active_columns)
    end where
	
    canopy_resistance(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%RS(1:number_of_active_columns)
    latent_heat_flux(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%LH(1:number_of_active_columns) / &
         (air_density(1:number_of_active_columns) * ConstLatHeatEvap)
    sensible_heat_flux(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%HFX(1:number_of_active_columns) / &
         (air_density(1:number_of_active_columns) * ConstHeatCapacAir)
    ground_heat_flux(active_source_column(1:number_of_active_columns)) = &
         -NoahmpIO%GRDFLX(1:number_of_active_columns)
    soil_upward_latent_heat_flux(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%EVGXY(1:number_of_active_columns) * NoahmpIO%VEGFRA(1:number_of_active_columns) + &
         NoahmpIO%EVBXY(1:number_of_active_columns) * &
         (1.0_kind_phys - NoahmpIO%VEGFRA(1:number_of_active_columns))
    transpiration_latent_heat_flux(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%TRXY(1:number_of_active_columns)
    canopy_upward_latent_heat_flux(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%EVCXY(1:number_of_active_columns)
    precipitation_advected_heat_flux(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%PAHXY(1:number_of_active_columns)
    wilting_soil_moisture(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%SMCWLT_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    reference_soil_moisture(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%SMCREF_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    normalized_soil_wetness(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%SMOIS(1,1:number_of_active_columns) / &
         NoahmpIO%SMCMAX_TABLE(NoahmpIO%ISLTYP(1:number_of_active_columns))
    soil_moisture_content(active_source_column(1:number_of_active_columns)) = &
         -NoahmpIO%ZSOIL(1) * NoahmpIO%SMOIS(1,1:number_of_active_columns) * 1000.0_kind_phys
		 
    do soil_layer = 2, vertical_dimension_of_soil
       soil_moisture_content(active_source_column(1:number_of_active_columns)) = &
            soil_moisture_content(active_source_column(1:number_of_active_columns)) + &
            (NoahmpIO%ZSOIL(soil_layer-1) - NoahmpIO%ZSOIL(soil_layer)) * &
            NoahmpIO%SMOIS(soil_layer,1:number_of_active_columns) * 1000.0_kind_phys
    end do
	
    snow_cover_fraction(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%SNOWC(1:number_of_active_columns)
    total_canopy_water(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%CANWAT(1:number_of_active_columns)
    ground_snowfall_rate(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%QSNOWXY(1:number_of_active_columns)
    surface_runoff(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%RUNSFXY(1:number_of_active_columns)
    surface_runoff_flux(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%RUNSFXY(1:number_of_active_columns) / &
         runoff_accumulation_time(1:number_of_active_columns)
    subsurface_runoff(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%RUNSBXY(1:number_of_active_columns) / &
         runoff_accumulation_time(1:number_of_active_columns)
    net_canopy_evaporation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ECANXY(1:number_of_active_columns)
    net_ground_evaporation(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%EDIRXY(1:number_of_active_columns)
    transpiration(active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%ETRANXY(1:number_of_active_columns)
    snow_deposition_sublimation_upward_latent_heat_flux( &
         active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%QSNSUBXY(1:number_of_active_columns) * ConstLatHeatSublim
    snow_freezing_rain_upward_latent_heat_flux( &
         active_source_column(1:number_of_active_columns)) = &
         NoahmpIO%QSNBOTXY(1:number_of_active_columns) * ConstLatHeatFusion
		 
    where (NoahmpIO%ALBEDO(1:number_of_active_columns) > 0.0_kind_phys)
       direct_visible_surface_albedo(active_source_column(1:number_of_active_columns)) = &
            NoahmpIO%ALBSFCDIRXY(1,1:number_of_active_columns)
       direct_nir_surface_albedo(active_source_column(1:number_of_active_columns)) = &
            NoahmpIO%ALBSFCDIRXY(2,1:number_of_active_columns)
       diffuse_visible_surface_albedo(active_source_column(1:number_of_active_columns)) = &
            NoahmpIO%ALBSFCDIFXY(1,1:number_of_active_columns)
       diffuse_nir_surface_albedo(active_source_column(1:number_of_active_columns)) = &
            NoahmpIO%ALBSFCDIFXY(2,1:number_of_active_columns)
    end where

    end associate

  end subroutine noahmp_run


  !> \section arg_table_noahmp_final Argument Table
  !! \htmlinclude noahmp_final.html
  !!
  subroutine noahmp_final(instance_number, errmsg, errflg)

    integer, intent(in)           :: instance_number
    character(len=*), intent(out) :: errmsg
    integer, intent(out)          :: errflg

    integer :: context_index

    errmsg = ''
    errflg = 0

    if (.not. allocated(NoahmpContexts)) return
    if (instance_number < 1 .or. instance_number > size(NoahmpContexts)) then
       errflg = 1
       errmsg = 'noahmp_final: instance_number is outside the allocated context range'
       return
    end if

    if (allocated(NoahmpContexts(instance_number)%io)) then
       deallocate(NoahmpContexts(instance_number)%io)
    end if

    do context_index = 1, size(NoahmpContexts)
       if (allocated(NoahmpContexts(context_index)%io)) return
    end do
    deallocate(NoahmpContexts)

  end subroutine noahmp_final

end module noahmp
