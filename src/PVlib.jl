module PVlib

using HTTP: get
using JSON: parse
using Dates: Date, DateTime, Second, @dateformat_str, UTC, dayofyear, value
using TimeZones: ZonedDateTime, TimeZone, Millisecond, month, day, Minute
using CSV: File
using Plots: plot
using Statistics: mean
using LinearAlgebra
using RecipesBase
using FLOWMath: abs_smooth, ksmax, ksmin
import ChainRulesCore: rrule, NoTangent
using ChainRulesCore: ignore_derivatives

const DEFAULT_ABS_DELTA = 1.0e-6
const DEFAULT_KS_HARDNESS = 50.0

pv_smooth_abs(x; delta::Real = DEFAULT_ABS_DELTA) = abs_smooth(x, delta * one(x))
pv_smooth_max(a, b; hardness::Real = DEFAULT_KS_HARDNESS) = ksmax([a, b], hardness)
pv_smooth_min(a, b; hardness::Real = DEFAULT_KS_HARDNESS) = ksmin([a, b], hardness)
pv_smooth_clamp(x, lo, hi; hardness::Real = DEFAULT_KS_HARDNESS) =
    pv_smooth_min(pv_smooth_max(x, lo; hardness = hardness), hi; hardness = hardness)
pv_smooth_step(x; hardness::Real = DEFAULT_KS_HARDNESS) = (one(x) + tanh(hardness * x)) / 2

datapath(parts...) = joinpath(pkgdir(@__MODULE__), "data", parts...)

include("weather.jl")
include("solar_position_algorithm.jl")
include("irradiance.jl")
include("temperature.jl")
include("solar_module.jl")
include("inverter.jl")
include("visualization.jl")
include("floating_utilities.jl")

export WeatherSample,
    TotalIrradiance,
    SolarPosition,
    ModuleTemperature,
    CellTemperature,
    SolarModule,
    SolarInverter,
    EffectiveIrradiance,
    DCComponents,
    ACPower,
    Panel,
    BoxObstacle
export get_meteorological_data_pvgis, get_meteorological_data_nsrdb, get_absolute_airmass
export get_solar_position
export get_projection, haydavies, get_extra_radiation, get_poa_ground_diffuse
export get_angle_of_incidence, get_poa_direct, get_total_irradiance
export sapm_module_temperature, sapm_cell_temperature
export read_solar_module, sapm_effective_irradiance, sapm_dc_components
export read_solar_inverter, sandia_ac_power
export panel_tilt_azimuth,
    panel_tilt_azimuth_3dof,
    get_ocean_surface_albedo,
    rolling_average_sapm_cell_temperature,
    get_shaded_fraction,
    get_power_norm,
    sapm_dc_components_shaded

end
