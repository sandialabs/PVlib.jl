# PVlib

## Overview

PVlib is a Julia package that provides functions and data structures for simulating the performance of photovoltaic (PV) energy systems and accomplishing related solar-resource and PV modeling tasks. PVlib is developed based on the widely-used pvlib python project ([pvlib/pvlib-python](https://github.com/pvlib/pvlib-python)) and intended to encompass a basic PV modeling approach while allowing for automatic differentiability within the Julia ecosystem.

The source code for PVlib is hosted on GitHub. Please see the [Installation](@ref) instructions in this documentation for help getting started.

For examples of how to use PVlib, see the quick start guide. This documentation assumes general familiarity with Julia, including common packages such as DataFrames and Dates/TimeZones. If you are new to Julia, the official Julia documentation and community tutorials are a good place to start.

PVlib follows a variable naming convention to promote consistency throughout the library (e.g., ghi, dni, dhi, poa_global, apparent_zenith).

This package also include functions for interfacing with floating structures in the `floating_utilities.jl` file.

## Installation

PVlib can be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run

```
pkg> add PVlib
```

## Quick start

## Example: PV power for Albuquerque, NM (Jan 1, 2023)

Load dependencies:

```@example pvlib_quick
using PVlib
using Dates
using TimeZones
using Plots
using JSON
```

Set site and time parameters:

```@example pvlib_quick
latitude = 35.1
longitude = -106.6
altitude = 1500.0 # m
start_date = Date(2023, 1, 1)
end_date = Date(2023, 1, 1)
tz = TimeZone("America/Denver")
```

Fetch weather data. Usually this would be done by accessing the PVGIS or NSRDB data sources. However, for documentation purposes, an example dataset has been saved and is loaded here:

```@example pvlib_quick
#weather_data = get_meteorological_data_pvgis(
#    latitude, longitude, start_date, end_date, tz, false
#);

weather_data_txt = read("../../data/weather_data.json", String)
weather_rows = JSON.parse(weather_data_txt)
weather_data = [
    WeatherSample(
        ZonedDateTime(r["time"]),
        isnothing(r["ghi"]) ? missing : r["ghi"],
        isnothing(r["dni"]) ? missing : r["dni"],
        isnothing(r["dhi"]) ? missing : r["dhi"],
        isnothing(r["temp_air"]) ? missing : r["temp_air"],
        isnothing(r["temp_dew"]) ? missing : r["temp_dew"],
        isnothing(r["relative_humidity"]) ? missing : r["relative_humidity"],
        isnothing(r["pressure"]) ? missing : r["pressure"],
        isnothing(r["wind_speed"]) ? missing : r["wind_speed"],
        isnothing(r["wind_direction"]) ? missing : r["wind_direction"],
        isnothing(r["albedo"]) ? missing : r["albedo"],
    )
    for r in weather_rows
]
```

Calculate solar positions:

```@example pvlib_quick
solar_pos = get_solar_position(latitude, longitude, altitude, weather_data);
```

Calculate plane-of-array irradiance for a south-facing array tilted at an angle equal to the latitude:

```@example pvlib_quick
surface_tilt = latitude
surface_azimuth = 180.0
total_irradiance = get_total_irradiance(
    surface_tilt, surface_azimuth, weather_data, solar_pos, 0.25
);
```

Load the desired PV module and inverter:

```@example pvlib_quick
module_filename = "sam-library-sandia-modules-2015-6-30.csv"
module_name = "Canadian Solar CS5P-220M [ 2009]"
pv_module = read_solar_module(module_name, module_filename);

inverter_filename = "sam-library-cec-inverters-2019-03-05.csv"
inverter_name = "ABB: MICRO-0.25-I-OUTD-US-208 [208V]"
pv_inverter = read_solar_inverter(inverter_name, inverter_filename);
```

Calculate the cell temperature and effective irradiance:

```@example pvlib_quick
cell_temp = sapm_cell_temperature(total_irradiance, weather_data);

effective_irradiance = sapm_effective_irradiance(
    total_irradiance, pv_module, solar_pos, surface_tilt, surface_azimuth, altitude
);
```

Calculate the DC power components and AC power:

```@example pvlib_quick
dc_components = sapm_dc_components(pv_module, effective_irradiance, cell_temp);
ac_power = sandia_ac_power(pv_inverter, dc_components);
```

Plot DC and AC power:

```@example pvlib_quick
plot(
    getfield.(dc_components, :time), 
    getfield.(dc_components, :p_mp),
    label="DC",
    xlabel="Time",
    ylabel="Power (W)",
    legend=:topright,
    xrotation=-30,
    bottom_margin=8Plots.mm,
)
plot!(getfield.(ac_power, :time), getfield.(ac_power, :ac_power); label="AC")
```

## Floating PV

This package was specifically built to enable modeling of floating solar modules in Julia. Functions to enable floating solar modeling include conversion from 6 degree-of-freedom motions to solar module tilt and azimuth and allowing for smaller timesteps to account for ocean wave timescale simulations.

## AD and Integration Notes

The current PV power path in `PVlib` is explicit and differentiable with respect to the numeric inputs used by `SIRENOpt`, including irradiance, surface tilt, surface azimuth, and array-area scaling applied outside the package. `FLOWMath`-based smooth approximations are used in the irradiance, SAPM, and inverter path so `ForwardDiff` can propagate gradients through the end-to-end AC power calculation.

`SIRENOpt` currently uses `PVlib` at the solar-array boundary: weather and solar position are mapped to photovoltaic DC or AC power output, and the higher-level ontology keeps generator and converter losses explicit in the system model. No implicit solve is used in this path.

## Reference

Refer to [pvlib/pvlib-python](https://github.com/pvlib/pvlib-python) for more details 
on the theory behind the photovoltaic modeling.

See the [API reference](@ref API) for all exported functions.
