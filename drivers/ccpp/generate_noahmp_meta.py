#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from ccpp_sources import sources as ccpp_sources


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "noahmp.F90"
OUTPUT = HERE / "noahmp.meta"


@dataclass(frozen=True)
class Argument:
    local_name: str
    type_name: str
    kind: str | None
    intent: str
    rank: int
    optional: bool


def preprocess_source(text: str, enabled_macros: set[str]) -> str:
    """Return the source seen with the requested interface macros enabled.

    Keeping this small preprocessor here makes the metadata generator consume
    the same public argument lists as the CCPP Fortran compilation.  It is not
    intended to replace cpp generally; the public interface currently uses
    only #ifdef, #ifndef, #else, and #endif.
    """
    frames: list[tuple[bool, bool]] = []
    active = True
    output: list[str] = []

    for line_number, line in enumerate(text.splitlines(keepends=True), start=1):
        directive = re.match(r"^\s*#\s*(ifdef|ifndef|else|endif)\b\s*(\w+)?", line)
        if directive is None:
            if active:
                output.append(line)
            continue

        keyword, macro = directive.groups()
        if keyword in {"ifdef", "ifndef"}:
            if macro is None:
                raise RuntimeError(f"Missing macro after #{keyword} on line {line_number}")
            condition = macro in enabled_macros
            if keyword == "ifndef":
                condition = not condition
            frames.append((active, condition))
            active = active and condition
        elif keyword == "else":
            if not frames:
                raise RuntimeError(f"Unmatched #else on line {line_number}")
            parent_active, condition = frames[-1]
            active = parent_active and not condition
        else:
            if not frames:
                raise RuntimeError(f"Unmatched #endif on line {line_number}")
            parent_active, _ = frames.pop()
            active = parent_active

    if frames:
        raise RuntimeError("Unterminated preprocessor conditional in noahmp.F90")

    return "".join(output)


def preprocess_ccpp_source(text: str) -> str:
    """Return the CCPP interface, with CCPP on and WRF_HYDRO off."""
    return preprocess_source(text, {"CCPP"})


def split_top_level(text: str) -> list[str]:
    fields: list[str] = []
    start = 0
    level = 0
    for index, character in enumerate(text):
        if character == "(":
            level += 1
        elif character == ")":
            level -= 1
        elif character == "," and level == 0:
            fields.append(text[start:index].strip())
            start = index + 1
    fields.append(text[start:].strip())
    return [field for field in fields if field]


def parse_scheme(source: str, scheme_name: str) -> list[Argument]:
    flattened = re.sub(r"&\s*\r?\n\s*", "", source)
    match = re.search(
        rf"(?ims)^\s*subroutine\s+{scheme_name}\s*\((.*?)\)\s*"
        rf"(.*?)^\s*end\s+subroutine\s+{scheme_name}\b",
        flattened,
    )
    if not match:
        raise RuntimeError(f"Unable to find {scheme_name} in {SOURCE}")

    ordered_names = split_top_level(match.group(1))
    body = match.group(2)
    declarations: dict[str, Argument] = {}
    declaration_pattern = re.compile(
        r"(?im)^\s*"
        r"(integer|logical|real\s*\(\s*kind\s*=\s*kind_phys\s*\)|"
        r"character\s*\(\s*len\s*=\s*\*\s*\))"
        r"\s*,\s*intent\s*\(\s*(inout|in|out)\s*\)"
        r"\s*(?:,\s*(optional))?\s*::\s*([^\r\n]+)"
    )

    for declaration in declaration_pattern.finditer(body):
        fortran_type = declaration.group(1).lower().replace(" ", "")
        intent = declaration.group(2).lower()
        if fortran_type.startswith("real"):
            type_name = "real"
            kind = "kind_phys"
        elif fortran_type.startswith("character"):
            type_name = "character"
            kind = "len=*"
        else:
            type_name = fortran_type
            kind = None

        optional = declaration.group(3) is not None
        for variable in split_top_level(declaration.group(4)):
            local_name = variable.split("(", maxsplit=1)[0].strip().lower()
            rank = variable.count(":")
            declarations[local_name] = Argument(
                local_name=local_name,
                type_name=type_name,
                kind=kind,
                intent=intent,
                rank=rank,
                optional=optional,
            )

    missing = [name for name in ordered_names if name.lower() not in declarations]
    if missing:
        raise RuntimeError(
            f"Missing declarations in {scheme_name}: {', '.join(missing)}"
        )

    arguments = [declarations[name.lower()] for name in ordered_names]
    if len({argument.local_name for argument in arguments}) != len(arguments):
        raise RuntimeError(f"Duplicate dummy argument in {scheme_name}")
    return arguments


# Configuration quantities whose run-phase names represent different fields.
INIT_STANDARD_NAME_OVERRIDES = {
    "dynamic_vegetation_option":
        "control_for_land_surface_scheme_dynamic_vegetation",
    "rain_snow_partition_option":
        "control_for_land_surface_scheme_precipitation_type_partition",
    "soil_water_transpiration_option":
        "control_for_land_surface_scheme_soil_moisture_factor_stomatal_resistance",
    "ground_resistance_evaporation_option":
        "control_for_land_surface_scheme_ground_evaporation_resistance",
    "stomata_resistance_option":
        "control_for_land_surface_scheme_canopy_stomatal_resistance",
    "snow_albedo_option":
        "control_for_land_surface_scheme_surface_snow_albedo",
    "canopy_radiation_transfer_option":
        "control_for_land_surface_scheme_radiative_transfer",
    "snow_soil_temperature_time_option":
        "control_for_land_surface_scheme_soil_and_snow_temperature_time_integration",
    "snow_thermal_conductivity_option":
        "control_for_land_surface_scheme_snow_thermal_conductivity",
    "soil_temperature_bottom_option":
        "control_for_land_surface_scheme_lower_boundary_soil_temperature",
    "soil_supercooled_water_option":
        "control_for_land_surface_scheme_supercooled_liquid_water",
    "frozen_soil_permeability_option":
        "control_for_land_surface_scheme_frozen_soil_permeability",
    "dynamic_vic_infiltration_option":
        "control_for_land_surface_scheme_dynamic_vic_infiltration",
    "tile_drainage_option":
        "control_for_land_surface_scheme_tile_drainage",
    "irrigation_option":
        "control_for_land_surface_scheme_irrigation",
    "irrigation_method_option":
        "control_for_land_surface_scheme_irrigation_method",
    "crop_model_option": "control_for_land_surface_scheme_crop_model",
    "soil_property_option":
        "control_for_land_surface_scheme_soil_parameter_treatment",
    "pedotransfer_option":
        "control_for_land_surface_scheme_pedotransfer_function",
    "snow_compaction_option":
        "control_for_land_surface_scheme_snow_compaction",
    "wetland_model_option":
        "control_for_land_surface_scheme_wetland",
    "snow_cover_fraction_option":
        "control_for_land_surface_scheme_snow_cover_fraction",
    "glacier_treatment_option":
        "control_for_land_surface_scheme_glacier_treatment",
    "urban_physics_option": "control_for_urban_physics",
    "surface_runoff_option":
        "control_for_land_surface_scheme_surface_runoff",
    "subsurface_runoff_option":
        "control_for_land_surface_scheme_subsurface_runoff",
    "surface_drag_option":
        "control_for_land_surface_scheme_surface_layer_drag_coefficient",
    "surface_stability_function_option":
        "control_for_surface_layer_scheme_stability_function",
    "surface_thermal_roughness_option":
        "control_for_thermal_roughness_lengths_over_land",
    "snicar_snow_shape_option":
        "control_for_land_surface_scheme_snicar_snow_shape",
    "snicar_rt_solver_option":
        "control_for_land_surface_scheme_snicar_radiative_transfer_solver",
    "snicar_band_number_option":
        "control_for_land_surface_scheme_snicar_spectral_bands",
    "snicar_solar_spectrum_option":
        "control_for_land_surface_scheme_snicar_solar_spectrum",
    "snicar_snow_optics_option":
        "control_for_land_surface_scheme_snicar_snow_optics",
    "snicar_dust_optics_option":
        "control_for_land_surface_scheme_snicar_dust_optics",
    "snicar_snow_bc_internal_mixing":
        "flag_for_snicar_internal_mixing_of_black_carbon_in_snow",
    "snicar_snow_dust_internal_mixing":
        "flag_for_snicar_internal_mixing_of_dust_in_snow",
    "snicar_use_aerosol": "flag_for_snicar_aerosol_effects",
    "snicar_use_organic_carbon": "flag_for_snicar_organic_carbon",
    "snicar_read_aerosol_table":
        "flag_for_snicar_aerosol_deposition_from_table",
    "water_category": "index_of_water_vegetation_category",
    "barren_category": "index_of_barren_vegetation_category",
    "ice_category": "index_of_ice_vegetation_category",
    "crop_category": "index_of_crop_vegetation_category",
    "evergreen_broadleaf_forest_category":
        "index_of_evergreen_broadleaf_forest_vegetation_category",
    "urban_category": "index_of_urban_vegetation_category",
    "natural_vegetation_category": "index_of_natural_vegetation_category",
    "urban_category_begin": "index_of_first_urban_vegetation_category",
    "default_crop_category": "index_of_default_crop_vegetation_category",
    "runoff_slope_category": "index_of_runoff_slope_category",
}

# A few descriptive local names would otherwise collide with registered CCPP
# quantities that have different physical units.  Keep the modular Noah-MP
# boundary explicit rather than silently advertising an incompatible field.
RUN_STANDARD_NAME_OVERRIDES = {
    "timestep_index": "index_of_timestep",
    "depth_of_soil_layer_interfaces": "depth_of_soil_layers",
    "active_snow_layer_lower_index": "number_of_snow_layers",
    "depth_of_snow_soil_layer_interfaces":
        "depth_from_snow_surface_at_bottom_interface",
    "vegetation_category": "vegetation_type_classification",
    "soil_category": "soil_type_classification",
    "dimensionless_snow_age": "dimensionless_age_of_surface_snow",
    "snow_albedo_on_previous_timestep":
        "surface_albedo_assuming_deep_snow_on_previous_timestep",
    "canopy_air_vapor_pressure": "air_vapor_pressure_in_canopy",
    "canopy_air_temperature": "air_temperature_in_canopy",
    "maximum_vegetation_fraction": "max_vegetation_area_fraction",
    "vegetation_fraction": "vegetation_area_fraction",
    "snow_temperature": "temperature_in_surface_snow",
    "canopy_liquid_water": "canopy_intercepted_liquid_water",
    "canopy_ice_water": "canopy_intercepted_ice_mass",
    "canopy_wet_fraction": "wet_canopy_area_fraction",
    "snow_water_equivalent":
        "water_equivalent_accumulated_snow_depth_over_land",
    "snow_water_equivalent_on_previous_timestep":
        "lwe_thickness_of_snowfall_amount_on_previous_timestep",
    "snow_layer_ice": "lwe_thickness_of_ice_in_surface_snow",
    "snow_layer_liquid_water":
        "lwe_thickness_of_liquid_water_in_surface_snow",
    "soil_liquid_water": "volume_fraction_of_unfrozen_water_in_soil",
    "soil_moisture": "volume_fraction_of_condensed_water_in_soil",
    "equilibrium_soil_moisture": "volumetric_equilibrium_soil_moisture",
    "atmospheric_forcing_height": "height_above_ground_at_lowest_model_layer",
    "air_temperature_at_lowest_model_layer":
        "air_temperature_at_surface_adjacent_layer",
    "specific_humidity_at_lowest_model_layer":
        "specific_humidity_at_surface_adjacent_layer",
    "eastward_wind_at_lowest_model_layer":
        "x_wind_at_surface_adjacent_layer",
    "northward_wind_at_lowest_model_layer":
        "y_wind_at_surface_adjacent_layer",
    "surface_exchange_coefficient_for_heat":
        "surface_drag_coefficient_for_heat_and_moisture_for_noahmp",
    "surface_exchange_coefficient_for_momentum":
        "surface_drag_coefficient_for_momentum_for_noahmp",
    "soil_moisture_below_soil_column":
        "volumetric_soil_moisture_between_soil_bottom_and_water_table",
    "aquifer_water_storage": "water_storage_in_aquifer",
    "soil_aquifer_water_storage":
        "water_storage_in_aquifer_and_saturated_soil",
    "lake_water_storage": "water_storage_in_lake",
    "leaf_mass": "leaf_mass_content",
    "root_mass": "fine_root_mass_content",
    "stem_mass": "stem_mass_content",
    "wood_mass": "wood_mass_content",
    "deep_soil_carbon_mass": "slow_soil_pool_mass_content_of_carbon",
    "shallow_soil_carbon_mass": "fast_soil_pool_mass_content_of_carbon",
    "ground_snowfall_rate": "lwe_snowfall_rate",
    "snow_depth": "surface_snow_thickness_water_equivalent_over_land",
    "surface_radiative_temperature": "surface_skin_temperature_over_land",
    "surface_emissivity": "surface_longwave_emissivity_over_land",
    "surface_roughness_length": "surface_roughness_length_over_land",
    "specific_humidity_at_surface": "surface_specific_humidity_over_land",
    "ground_heat_flux": "upward_heat_flux_in_soil_over_land",
    "subsurface_runoff": "subsurface_runoff_flux",
    "latent_heat_flux":
        "kinematic_surface_upward_latent_heat_flux_over_land",
    "sensible_heat_flux":
        "kinematic_surface_upward_sensible_heat_flux_over_land",
    "surface_drag_wind_speed_for_momentum":
        "surface_drag_wind_speed_for_momentum_in_air_over_land",
    "surface_drag_mass_flux_for_heat_and_moisture":
        "surface_drag_mass_flux_for_heat_and_moisture_in_air_over_land",
    "precipitation_advected_heat_flux": "total_precipitation_advected_heat",
    "net_canopy_evaporation": "evaporation_of_intercepted_water",
    "net_ground_evaporation": "soil_surface_evaporation_rate",
    "transpiration": "transpiration_rate",
    "wilting_soil_moisture":
        "volume_fraction_of_condensed_water_in_soil_at_wilting_point",
    "reference_soil_moisture":
        "threshold_volume_fraction_of_condensed_water_in_soil",
    "surface_albedo":
        "surface_albedo_for_diffused_shortwave_on_radiation_timestep",
    "canopy_resistance": "aerodynamic_resistance_in_canopy",
    "snow_cover_fraction": "surface_snow_area_fraction_over_land",
    "total_canopy_water": "canopy_water_amount",
    "transpiration_latent_heat_flux": "transpiration_flux",
    "direct_visible_surface_albedo":
        "surface_albedo_direct_visible_over_land",
    "direct_nir_surface_albedo": "surface_albedo_direct_NIR_over_land",
    "diffuse_visible_surface_albedo":
        "surface_albedo_diffuse_visible_over_land",
    "diffuse_nir_surface_albedo": "surface_albedo_diffuse_NIR_over_land",
    # These SCM standard names exceed Fortran's 63-character identifier
    # limit, so the scheme uses shorter local argument names.
    "direct_nir_shortwave_flux_at_surface":
        "surface_downwelling_direct_nir_shortwave_flux_on_radiation_timestep",
    "diffuse_nir_shortwave_flux_at_surface":
        "surface_downwelling_diffuse_nir_shortwave_flux_on_radiation_timestep",
    "direct_vis_shortwave_flux_at_surface":
        "surface_downwelling_direct_uv_and_vis_shortwave_flux_on_radiation_timestep",
    "diffuse_vis_shortwave_flux_at_surface":
        "surface_downwelling_diffuse_uv_and_vis_shortwave_flux_on_radiation_timestep",
    "precipitation_amount_on_dynamics_timestep":
        "nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep",
    "total_precipitation_rate_at_surface":
        "total_precipitation_rate_at_surface_for_noahmp",
}


def standard_name(argument: Argument, scheme_name: str) -> str:
    if argument.local_name == "errmsg":
        return "ccpp_error_message"
    if argument.local_name == "errflg":
        return "ccpp_error_code"
    if scheme_name == "noahmp_init":
        return INIT_STANDARD_NAME_OVERRIDES.get(
            argument.local_name, argument.local_name
        )
    if scheme_name == "noahmp_run":
        return RUN_STANDARD_NAME_OVERRIDES.get(
            argument.local_name, argument.local_name
        )
    return argument.local_name


def long_name(argument: Argument) -> str:
    if argument.local_name == "errmsg":
        return "error message for error handling in CCPP"
    if argument.local_name == "errflg":
        return "error code for error handling in CCPP"
    return argument.local_name.replace("_", " ")


SOIL_LAYERED = {
    "soil_temperature",
    "soil_liquid_water",
    "soil_moisture",
    "equilibrium_soil_moisture",
    "accumulated_soil_layer_transpiration",
    "soil_quartz_content",
    "soil_exponent_b",
    "dry_soil_moisture",
    "wilting_soil_moisture",
    "reference_soil_moisture",
    "saturated_soil_moisture",
    "saturated_soil_hydraulic_conductivity",
    "saturated_soil_water_diffusivity",
    "saturated_soil_matric_potential",
}

SNOW_LAYERED = {
    "snow_temperature",
    "snow_layer_ice",
    "snow_layer_liquid_water",
}

RADIATION_BANDED = {
    "direct_soil_albedo",
    "diffuse_soil_albedo",
    "direct_surface_albedo",
    "diffuse_surface_albedo",
    "direct_snow_albedo",
    "diffuse_snow_albedo",
}


def dimensions(argument: Argument, scheme_name: str) -> str:
    if argument.rank == 0:
        return "()"
    if scheme_name != "noahmp_run":
        raise RuntimeError(
            f"Unexpected array {argument.local_name} outside the run phase"
        )
    if argument.local_name == "depth_of_soil_layer_interfaces":
        return "(vertical_dimension_of_soil_internal_to_land_surface_scheme)"
    if argument.local_name == "soil_category_at_layer":
        return "(horizontal_loop_extent,4)"
    if argument.local_name == "equilibrium_soil_moisture":
        return (
            "(horizontal_loop_extent,"
            "vertical_dimension_of_soil_internal_to_land_surface_scheme)"
        )
    if argument.local_name == "depth_of_snow_soil_layer_interfaces":
        return (
            "(horizontal_loop_extent,"
            "lower_bound_of_vertical_dimension_of_surface_snow:"
            "vertical_dimension_of_soil_internal_to_land_surface_scheme)"
        )
    if argument.local_name in SOIL_LAYERED and argument.rank == 2:
        return "(horizontal_loop_extent,vertical_dimension_of_soil)"
    if argument.local_name in SNOW_LAYERED:
        return (
            "(horizontal_loop_extent,"
            "lower_bound_of_vertical_dimension_of_surface_snow:"
            "upper_bound_of_vertical_dimension_of_surface_snow)"
        )
    if argument.local_name in RADIATION_BANDED:
        # Noah-MP and both legacy and modular table readers fix this extent at
        # two bands. A literal keeps the CCPP host interface free of a
        # redundant dimension variable.
        return "(horizontal_loop_extent,2)"
    if argument.rank == 1:
        return "(horizontal_loop_extent)"
    raise RuntimeError(f"No dimension mapping for {argument.local_name}")


REAL_UNITS = {
    "timestep_for_physics": "s",
    "leaf_area_index": "none",
    "stem_area_index": "none",
    "instantaneous_cosine_of_zenith_angle": "none",
    "active_snow_layer_lower_index": "count",
    "dimensionless_snow_age": "none",
    "snow_albedo_on_previous_timestep": "frac",
    "canopy_wet_fraction": "none",
    "maximum_vegetation_fraction": "frac",
    "vegetation_fraction": "frac",
    # Time and coordinates.
    "forecast_julian_day": "days",
    "precipitation_amount_on_dynamics_timestep": "m",
    "direct_nir_shortwave_flux_at_surface": "W m-2",
    "diffuse_nir_shortwave_flux_at_surface": "W m-2",
    "direct_vis_shortwave_flux_at_surface": "W m-2",
    "diffuse_vis_shortwave_flux_at_surface": "W m-2",
    "snicar_hydrophobic_black_carbon_deposition_flux": "kg m-2 s-1",
    "snicar_hydrophilic_black_carbon_deposition_flux": "kg m-2 s-1",
    "snicar_hydrophobic_organic_carbon_deposition_flux": "kg m-2 s-1",
    "snicar_hydrophilic_organic_carbon_deposition_flux": "kg m-2 s-1",
    "snicar_dust_species_1_deposition_flux": "kg m-2 s-1",
    "snicar_dust_species_2_deposition_flux": "kg m-2 s-1",
    "snicar_dust_species_3_deposition_flux": "kg m-2 s-1",
    "snicar_dust_species_4_deposition_flux": "kg m-2 s-1",
    "snicar_dust_species_5_deposition_flux": "kg m-2 s-1",
    "characteristic_grid_lengthscale": "m",
    "x_direction_grid_spacing": "m",
    "y_direction_grid_spacing": "m",
    "soil_temperature_bottom_depth": "m",
    "depth_of_soil_layer_interfaces": "m",
    "depth_of_snow_soil_layer_interfaces": "m",
    "atmospheric_forcing_height": "m",
    "snow_depth": "mm",
    "water_table_depth": "m",
    "mean_capillary_drive": "m",
    "tile_drain_depth": "m",
    "tile_radius": "m",
    "tile_spacing": "m",
    "maximum_wetland_storage": "m",
    "surface_roughness_length": "cm",
    "atmosphere_bottom_layer_thickness": "m",
    "atmosphere_boundary_layer_thickness": "m",
    "latitude": "radian",
    "planting_day": "d",
    "harvest_day": "d",
    "irrigation_stop_days_before_harvest": "d",
    "growing_degree_day": "K d",
    "seasonal_growing_degree_days": "K d",

    # Atmospheric and surface thermodynamic quantities.
    "specific_humidity_at_surface": "kg kg-1",
    "canopy_air_vapor_pressure": "Pa",
    "output_specific_humidity": "kg kg-1",
    "vegetated_water_vapor_mixing_ratio_at_2m": "kg kg-1",
    "bare_ground_water_vapor_mixing_ratio_at_2m": "kg kg-1",
    "specific_humidity_at_lowest_model_layer": "kg kg-1",
    "output_air_pressure": "Pa",
    "surface_air_pressure": "Pa",
    "air_pressure_at_surface_adjacent_layer": "Pa",
    "output_wind_speed": "m s-1",
    "eastward_wind_at_lowest_model_layer": "m s-1",
    "northward_wind_at_lowest_model_layer": "m s-1",

    # Water and snow state.
    "canopy_liquid_water": "mm",
    "canopy_ice_water": "mm",
    "snow_water_equivalent": "mm",
    "snow_water_equivalent_on_previous_timestep": "mm",
    "flood_irrigation_water_amount": "m",
    "micro_irrigation_water_amount": "m",
    "sprinkler_irrigation_water_amount": "m",
    "soil_moisture_below_soil_column": "m3 m-3",
    "aquifer_water_storage": "mm",
    "soil_aquifer_water_storage": "mm",
    "lake_water_storage": "mm",
    "wetland_water_storage": "mm",
    "snow_layer_ice": "mm",
    "snow_layer_liquid_water": "mm",
    "accumulated_soil_surface_evaporation": "m s-1",
    "accumulated_soil_surface_inflow": "m s-1",
    "accumulated_surface_water_change": "mm",
    "accumulated_precipitation": "mm",
    "accumulated_canopy_evaporation": "mm",
    "accumulated_transpiration": "mm",
    "accumulated_ground_evaporation": "mm",
    "accumulated_soil_layer_transpiration": "m s-1",
    "accumulated_glacier_excess_flow": "mm",
    "accumulated_surface_runoff": "m",
    "accumulated_subsurface_runoff": "m",
    "accumulated_tile_drainage": "mm",
    "accumulated_snowfall": "mm",
    "accumulated_snowmelt": "mm",
    "accumulated_shallow_groundwater_recharge": "mm",
    "accumulated_deep_groundwater_recharge": "mm",
    "accumulated_sprinkler_irrigation": "mm",
    "accumulated_micro_irrigation": "mm",
    "accumulated_flood_irrigation": "mm",
    "accumulated_sprinkler_evaporation_loss": "m",

    # Vegetation and carbon state/fluxes.
    "leaf_mass": "g m-2",
    "root_mass": "g m-2",
    "stem_mass": "g m-2",
    "wood_mass": "g m-2",
    "deep_soil_carbon_mass": "g m-2",
    "shallow_soil_carbon_mass": "g m-2",
    "grain_mass": "g m-2",
    "net_ecosystem_exchange": "g m-2 s-1",
    "gross_primary_productivity": "g m-2 s-1",
    "net_primary_productivity": "g m-2 s-1",
    "total_photosynthesis": "umol m-2 s-1",

    # Soil hydraulic parameters.
    "soil_liquid_water": "frac",
    "soil_moisture": "frac",
    "equilibrium_soil_moisture": "m3 m-3",
    "dry_soil_moisture": "m3 m-3",
    "wilting_soil_moisture": "frac",
    "reference_soil_moisture": "frac",
    "saturated_soil_moisture": "m3 m-3",
    "saturated_soil_hydraulic_conductivity": "m s-1",
    "saturated_soil_water_diffusivity": "m2 s-1",
    "saturated_soil_matric_potential": "m",
    "reference_soil_conductivity": "m s-1",
    "sprinkler_irrigation_rate": "mm h-1",
    "micro_irrigation_rate": "mm h-1",
    "irrigation_stop_precipitation_threshold": "mm",
    "simple_tile_drainage_coefficient": "mm d-1",
    "hooghoudt_tile_drainage_coefficient": "m d-1",

    # Diagnostic energy and exchange quantities.
    "surface_exchange_coefficient_for_heat": "none",
    "surface_exchange_coefficient_for_momentum": "none",
    "heat_exchange_coefficient_above_canopy": "m s-1",
    "heat_exchange_coefficient_bare_ground": "m s-1",
    "leaf_heat_exchange_coefficient": "m s-1",
    "under_canopy_heat_exchange_coefficient": "m s-1",
    "vegetated_heat_exchange_coefficient_at_2m": "m s-1",
    "bare_ground_heat_exchange_coefficient_at_2m": "m s-1",
    "sunlit_stomatal_resistance": "s m-1",
    "shaded_stomatal_resistance": "s m-1",
    "canopy_resistance": "s m-1",
    "surface_emissivity": "frac",
    "surface_albedo": "frac",
    "snow_cover_fraction": "frac",
    "direct_visible_surface_albedo": "frac",
    "direct_nir_surface_albedo": "frac",
    "diffuse_visible_surface_albedo": "frac",
    "diffuse_nir_surface_albedo": "frac",
    "soil_energy_storage": "kJ m-2",
    "snow_energy_storage": "kJ m-2",
    "total_evapotranspiration": "kg m-2 s-1",
    "total_soil_moisture": "mm",
    "soil_moisture_content": "kg m-2",
    "normalized_soil_wetness": "frac",
    "total_canopy_water": "kg m-2",
    "surface_runoff": "kg m-2",
    "temperature_at_2m_from_noahmp": "K",
    "specific_humidity_at_2m_from_noahmp": "kg kg-1",
    "surface_drag_wind_speed_for_momentum": "m s-1",
    "surface_drag_mass_flux_for_heat_and_moisture": "kg m-2 s-1",
    "latent_heat_flux": "kg kg-1 m s-1",
    "sensible_heat_flux": "K m s-1",
    "subsurface_runoff": "kg m-2 s-1",
    "surface_runoff_flux": "kg m-2 s-1",
    "net_canopy_evaporation": "kg m-2 s-1",
    "net_ground_evaporation": "kg m-2 s-1",
    "transpiration": "kg m-2 s-1",
    "surface_ponding": "mm",
}


TEMPERATURE_FIELDS = {
    "ground_temperature",
    "canopy_temperature",
    "canopy_air_temperature",
    "soil_temperature",
    "snow_temperature",
    "output_air_temperature",
    "surface_radiative_temperature",
    "vegetated_air_temperature_at_2m",
    "bare_ground_air_temperature_at_2m",
    "temperature_at_2m_from_noahmp",
    "vegetated_ground_temperature",
    "bare_ground_temperature",
    "deep_soil_temperature",
    "air_temperature_at_lowest_model_layer",
}

HEAT_FLUX_FIELDS = {
    "accumulated_ground_heat_flux",
    "sprinkler_heat_accumulation",
    "latent_heat_flux",
    "sensible_heat_flux",
    "ground_heat_flux",
    "absorbed_shortwave_radiation",
    "net_longwave_radiation",
    "absorbed_photosynthetic_radiation",
    "vegetation_absorbed_shortwave_radiation",
    "ground_absorbed_shortwave_radiation",
    "canopy_net_longwave_radiation",
    "vegetated_ground_net_longwave_radiation",
    "canopy_sensible_heat_flux",
    "vegetated_ground_sensible_heat_flux",
    "vegetated_ground_latent_heat_flux",
    "vegetated_ground_heat_flux",
    "bare_ground_net_longwave_radiation",
    "bare_ground_sensible_heat_flux",
    "bare_ground_latent_heat_flux",
    "bare_ground_heat_flux",
    "transpiration_latent_heat_flux",
    "canopy_evaporation_latent_heat_flux",
    "canopy_heat_storage_change",
    "precipitation_advected_heat_flux",
    "vegetated_ground_precipitation_advected_heat_flux",
    "canopy_precipitation_advected_heat_flux",
    "bare_ground_precipitation_advected_heat_flux",
    "bottom_soil_heat_flux",
    "surface_downwelling_shortwave_flux",
    "surface_downwelling_longwave_flux",
    "soil_upward_latent_heat_flux",
    "canopy_upward_latent_heat_flux",
    "snow_deposition_sublimation_upward_latent_heat_flux",
    "snow_freezing_rain_upward_latent_heat_flux",
}

WATER_FLUX_FIELDS = {
    "ground_snowfall_rate",
    "ground_rainfall_rate",
    "net_canopy_evaporation",
    "net_ground_evaporation",
    "transpiration",
    "canopy_snow_interception",
    "canopy_rain_interception",
    "canopy_snow_drip",
    "canopy_rain_drip",
    "snow_throughfall",
    "rain_throughfall",
    "snow_sublimation",
    "snow_frost",
    "canopy_ice_sublimation",
    "canopy_ice_frost",
    "canopy_liquid_evaporation",
    "canopy_liquid_dew",
    "canopy_liquid_freeze",
    "canopy_ice_melt",
    "snow_bottom_outflow",
    "ground_snowmelt",
    "rainfall_at_reference_height",
    "snowfall_at_reference_height",
    "convective_precipitation_rate_on_previous_timestep",
    "explicit_precipitation_rate_on_previous_timestep",
    "shallow_convective_precipitation_rate_on_previous_timestep",
    "total_precipitation_rate_at_surface",
    "snowfall_rate_on_previous_timestep",
    "graupel_precipitation_rate_on_previous_timestep",
    "ice_precipitation_rate_on_previous_timestep",
    "hail_precipitation_rate_on_previous_timestep",
}


def units(argument: Argument) -> str:
    name = argument.local_name
    if name == "errmsg":
        return "none"
    if name == "errflg":
        return "1"
    if argument.type_name == "character":
        return "none"
    if argument.type_name == "logical":
        return "flag"
    if argument.type_name == "integer":
        if name == "timestep_for_physics":
            return "s"
        if name == "number_of_days_in_current_year":
            return "days"
        if name in {
            "lower_bound_of_vertical_dimension_of_surface_snow",
            "upper_bound_of_vertical_dimension_of_surface_snow",
        }:
            return "count"
        if (
            name.startswith("number_of_")
            or name in {
                "horizontal_dimension",
                "horizontal_loop_extent",
                "vertical_dimension_of_soil",
                "vertical_dimension_of_surface_snow",
                "sprinkler_irrigation_event_count",
                "micro_irrigation_event_count",
                "flood_irrigation_event_count",
            }
        ):
            return "count"
        if name == "surface_thermal_roughness_option":
            return "1"
        return "index"
    if name in REAL_UNITS:
        return REAL_UNITS[name]
    if name in TEMPERATURE_FIELDS:
        return "K"
    if name in HEAT_FLUX_FIELDS:
        return "W m-2"
    if name in WATER_FLUX_FIELDS:
        return "mm s-1"
    if name == "soil_category_at_layer":
        return "index"
    if name == "precipitation_type":
        return "flag"
    # Remaining real quantities are fractions, albedos, indices expressed as
    # reals, or empirical dimensionless parameters.
    return "1"


def render_table(scheme_name: str, arguments: list[Argument]) -> list[str]:
    lines = [
        "########################################################################",
        "[ccpp-arg-table]",
        f"  name = {scheme_name}",
        "  type = scheme",
    ]
    for argument in arguments:
        lines.extend(
            [
                f"[{argument.local_name}]",
                f"  standard_name = {standard_name(argument, scheme_name)}",
                f"  long_name = {long_name(argument)}",
                f"  units = {units(argument)}",
                f"  dimensions = {dimensions(argument, scheme_name)}",
                f"  type = {argument.type_name}",
            ]
        )
        if argument.kind is not None:
            lines.append(f"  kind = {argument.kind}")
        lines.append(f"  intent = {argument.intent}")
        if argument.optional:
            lines.append("  optional = True")
    return lines


def validate(arguments: dict[str, list[Argument]], metadata: str) -> None:
    table_pattern = re.compile(
        r"(?ms)^\[ccpp-arg-table\]\s*\n\s*name\s*=\s*(\w+)\s*\n"
        r"\s*type\s*=\s*scheme\s*\n(.*?)(?=^#{8,}|\Z)"
    )
    metadata_tables: dict[str, list[str]] = {}
    for table in table_pattern.finditer(metadata):
        metadata_tables[table.group(1)] = re.findall(
            r"(?m)^\[([^\]]+)\]\s*$", table.group(2)
        )

    for scheme_name, source_arguments in arguments.items():
        source_names = [argument.local_name for argument in source_arguments]
        metadata_names = metadata_tables.get(scheme_name)
        if metadata_names is None:
            raise RuntimeError(f"Missing metadata table for {scheme_name}")
        if metadata_names != source_names:
            raise RuntimeError(
                f"Argument order mismatch for {scheme_name}:\n"
                f"source={source_names}\nmetadata={metadata_names}"
            )
        for argument in source_arguments:
            if argument.type_name == "real" and argument.kind != "kind_phys":
                raise RuntimeError(
                    f"Real argument {argument.local_name} does not use kind_phys"
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check", action="store_true",
        help="verify the existing metadata matches the interface and source dependencies without writing it",
    )
    options = parser.parse_args()
    raw_source = SOURCE.read_text(encoding="utf-8")
    # CCPP prebuild reads noahmp.meta only through these source annotations.
    # A matching metadata file is insufficient if an annotation was removed.
    for scheme_name in ("noahmp_init", "noahmp_run", "noahmp_final"):
        anchor = (
            rf"(?m)^\s*!>\s*\\section\s+arg_table_{scheme_name}\b[^\n]*\n"
            rf"\s*!!\s*\\htmlinclude\s+{scheme_name}\.html\s*$"
        )
        if not re.search(anchor, raw_source):
            raise RuntimeError(f"Missing CCPP metadata annotation for {scheme_name}")
    source = preprocess_ccpp_source(raw_source)
    schemes = {
        name: parse_scheme(source, name)
        for name in ("noahmp_init", "noahmp_run", "noahmp_final")
    }
    expected_counts = {
        "noahmp_init": 45,
        "noahmp_run": 140,
        "noahmp_final": 2,
    }
    actual_counts = {name: len(args) for name, args in schemes.items()}
    if actual_counts != expected_counts:
        raise RuntimeError(
            f"Unexpected lifecycle argument counts: {actual_counts}; "
            f"expected {expected_counts}"
        )

    lines = [
        "[ccpp-table-properties]",
        "  name = noahmp",
        "  type = scheme",
        "",
        "# Base modular Noah-MP interface. WRF_HYDRO is intentionally disabled.",
        "# Scheme and host metadata must use identical standard names and units.",
    ]
    # Only sources reachable from this scheme are registered. Machine is a
    # host dependency, so the HRLDAS stand-in must never enter this closure.
    import os
    dependencies = [os.path.relpath(path, HERE).replace(os.sep, "/")
                    for path in ccpp_sources() if path != SOURCE]
    lines[3:3] = [f"  dependencies = {path}" for path in dependencies]
    for scheme_name, scheme_arguments in schemes.items():
        lines.extend(render_table(scheme_name, scheme_arguments))
    metadata = "\n".join(lines) + "\n"
    validate(schemes, metadata)
    if options.check:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != metadata:
            raise RuntimeError(
                f"{OUTPUT} is missing or stale; run generate_noahmp_meta.py without --check to regenerate it"
            )
    else:
        OUTPUT.write_text(metadata, encoding="utf-8", newline="\n")

    counts = ", ".join(
        f"{name}={len(args)}" for name, args in schemes.items()
    )
    action = "Verified" if options.check else "Wrote"
    print(f"{action} {OUTPUT} ({counts})")


if __name__ == "__main__":
    main()
