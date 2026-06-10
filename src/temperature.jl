"""
    ModuleTemperature{T}

Module back-surface temperature for a single timestamp.

Typically computed from POA irradiance and weather inputs using
`sapm_module_temperature`.

# Fields:
* `time`: Timestamp
* `module_temperature`: Module temperature (°C)
"""
struct ModuleTemperature{T}
    time::Any
    module_temperature::T
end

function rrule(::Type{ModuleTemperature}, time, module_temperature)
    y = ModuleTemperature{typeof(module_temperature)}(time, module_temperature)

    function pullback(ȳ)
        return NoTangent(), NoTangent(), getproperty(ȳ, :module_temperature)
    end

    return y, pullback
end

"""
    CellTemperature{T}

Cell temperature for a single timestamp.

Typically computed from module temperature and POA irradiance using
`sapm_cell_temperature`.

# Fields:
* `time`: Timestamp
* `cell_temperature`: Cell temperature (°C)
"""
struct CellTemperature{T}
    time::Any
    cell_temperature::T
end

function rrule(::Type{CellTemperature}, time, cell_temperature)
    y = CellTemperature{typeof(cell_temperature)}(time, cell_temperature)

    function pullback(ȳ)
        return NoTangent(), NoTangent(), getproperty(ȳ, :cell_temperature)
    end

    return y, pullback
end

"""
    sapm_module_temperature(
        total_irradiance::TotalIrradiance,
        weather_data::WeatherSample;
        a::Real=-3.47,
        b::Real=-0.0594
    ) -> ModuleTemperature

    sapm_module_temperature(
        total_irradiance::AbstractVector{<:TotalIrradiance},
        weather_data::AbstractVector{<:WeatherSample};
        a::Real=-3.47,
        b::Real=-0.0594
    ) -> Vector{<:ModuleTemperature}

    sapm_module_temperature(
        total_irradiance::AbstractVector{<:TotalIrradiance},
        weather_data::WeatherSample;
        a::Real=-3.47,
        b::Real=-0.0594
    ) -> Vector{<:ModuleTemperature}

Compute module back-surface temperature using the Sandia Module Temperature Model.

# Arguments:
* `total_irradiance`: POA irradiance input(s) containing at least `poa_global` (W/m^2)
* `weather_data`: Weather input(s) containing at least `wind_speed` (m/s) and `temp_air` (°C)
* `a`, `b`: Empirical coefficients for the temperature model

# Returns:
* Module temperature(s) (°C), as `ModuleTemperature` or `Vector{<:ModuleTemperature}`

```jldoctest
julia> using TimeZones

julia> total_irradiance = TotalIrradiance(ZonedDateTime(2023, 1, 1, tz"America/Denver"),
           1000.0, 0.0, 0.0, 0.0, 0.0);

julia> weather_data = WeatherSample{Float64}(
            ZonedDateTime(2023, 1, 1, tz"America/Denver"),
            1.0, 1.0, 1.0,               # ghi, dni, dhi (unused here)
            25.0,                        # temp_air (C)
            1.0,                         # relative_humidity (unused here)
            20.0,                        # temp_dewpoint (unused here)
            101325.0,                    # pressure (unused here)
            1.0,                         # wind_speed (m/s)
            1.0,                         # wind_direction (unused here)
            0.1,                         # albedo (unused here)
        );

julia> module_temp = sapm_module_temperature(total_irradiance, weather_data);

julia> isapprox(module_temp.module_temperature, 54.322504092500964; atol=1e-6)
true
```
"""
function sapm_module_temperature(
    total_irradiance::TotalIrradiance,
    weather_data::WeatherSample;
    a::Real = -3.47,
    b::Real = -0.0594,
)

    module_temperature =
        total_irradiance.poa_global * exp(a + b * weather_data.wind_speed) +
        weather_data.temp_air

    return ModuleTemperature(weather_data.time, module_temperature)
end

function sapm_module_temperature(
    total_irradiance::AbstractVector{<:TotalIrradiance},
    weather_data::AbstractVector{<:WeatherSample};
    a::Real = -3.47,
    b::Real = -0.0594,
)

    n = length(total_irradiance)
    length(weather_data) == n || throw(
        DimensionMismatch(
            "weather_data length $(length(weather_data)) != total_irradiance length $n",
        ),
    )

    return map(1:n) do i
        sapm_module_temperature(total_irradiance[i], weather_data[i]; a = a, b = b)
    end
end

function sapm_module_temperature(
    total_irradiance::AbstractVector{<:TotalIrradiance},
    weather_data::WeatherSample;
    a::Real = -3.47,
    b::Real = -0.0594,
)

    n = length(total_irradiance)

    return map(1:n) do ii
        sapm_module_temperature(total_irradiance[ii], weather_data; a = a, b = b)
    end
end

"""
    sapm_cell_temperature(
        total_irradiance::TotalIrradiance,
        weather_data::WeatherSample;
        a::Real=-3.47,
        b::Real=-0.0594,
        deltaT::Real=3.0,
        irrad_ref::Real=1000.0
    ) -> CellTemperature

    sapm_cell_temperature(
        total_irradiance::AbstractVector{<:TotalIrradiance},
        weather_data::AbstractVector{<:WeatherSample};
        a::Real=-3.47,
        b::Real=-0.0594,
        deltaT::Real=3.0,
        irrad_ref::Real=1000.0
    ) -> Vector{<:CellTemperature}

    sapm_cell_temperature(
        total_irradiance::AbstractVector{<:TotalIrradiance},
        weather_data::WeatherSample;
        a::Real=-3.47,
        b::Real=-0.0594,
        deltaT::Real=3.0,
        irrad_ref::Real=1000.0
    ) -> Vector{<:CellTemperature}

Compute SAPM cell temperature from SAPM module temperature and POA irradiance.

# Arguments:
* `total_irradiance`: POA irradiance input(s) containing at least `poa_global` (W/m^2)
* `weather_data`: Weather input(s) containing at least `wind_speed` (m/s) and `temp_air` (°C)
* `a`, `b`: Empirical coefficients for the module temperature model
* `deltaT`: Cell-to-module temperature difference at reference irradiance (°C)
* `irrad_ref`: Reference irradiance used for scaling `deltaT` (W/m^2)

# Returns:
* Cell temperature(s) (°C), as `CellTemperature` or `Vector{<:CellTemperature}`

```jldoctest
julia> using TimeZones

julia> total_irradiance = TotalIrradiance(ZonedDateTime(2023, 1, 1, tz"America/Denver"),
           1000.0, 0.0, 0.0, 0.0, 0.0);

julia> weather_data = WeatherSample{Float64}(
            ZonedDateTime(2023, 1, 1, tz"America/Denver"),
            1.0, 1.0, 1.0,               # ghi, dni, dhi (unused here)
            25.0,                        # temp_air (C)
            1.0,                         # relative_humidity (unused here)
            20.0,                        # temp_dewpoint (unused here)
            101325.0,                    # pressure (unused here)
            1.0,                         # wind_speed (m/s)
            1.0,                         # wind_direction (unused here)
            0.1,                         # albedo (unused here)
        );

julia> cell_temp = sapm_cell_temperature(total_irradiance, weather_data);

julia> isapprox(cell_temp.cell_temperature, 57.322504092500964; atol=1e-6)
true
```
"""
function sapm_cell_temperature(
    total_irradiance::TotalIrradiance,
    weather_data::WeatherSample;
    a::Real = -3.47,
    b::Real = -0.0594,
    deltaT::Real = 3.0,
    irrad_ref::Real = 1000.0,
)

    module_temperature =
        sapm_module_temperature(total_irradiance, weather_data; a = a, b = b)
    cell_temperature =
        module_temperature.module_temperature +
        (total_irradiance.poa_global / irrad_ref) * deltaT

    return CellTemperature(total_irradiance.time, cell_temperature)
end


function sapm_cell_temperature(
    total_irradiance::AbstractVector{<:TotalIrradiance},
    weather_data::AbstractVector{<:WeatherSample};
    a::Real = -3.47,
    b::Real = -0.0594,
    deltaT::Real = 3.0,
    irrad_ref::Real = 1000.0,
)

    [ti.time for ti in total_irradiance] == [wd.time for wd in weather_data] || throw(
        ArgumentError(
            "Time mismatch: weather=$([wd.time for wd in weather_data]) irradiance=$([ti.time for ti in total_irradiance])",
        ),
    )

    n = length(total_irradiance)
    length(weather_data) == n || throw(
        DimensionMismatch(
            "weather_data length $(length(weather_data)) != total_irradiance length $n",
        ),
    )

    return map(1:n) do i
        sapm_cell_temperature(
            total_irradiance[i],
            weather_data[i];
            a = a,
            b = b,
            deltaT = deltaT,
            irrad_ref = irrad_ref,
        )
    end
end

function sapm_cell_temperature(
    total_irradiance::AbstractVector{<:TotalIrradiance},
    weather_data::WeatherSample;
    a::Real = -3.47,
    b::Real = -0.0594,
    deltaT::Real = 3.0,
    irrad_ref::Real = 1000.0,
)

    n = length(total_irradiance)

    return map(1:n) do ii
        sapm_cell_temperature(
            total_irradiance[ii],
            weather_data;
            a = a,
            b = b,
            deltaT = deltaT,
            irrad_ref = irrad_ref,
        )
    end
end

# ModuleTemperature pretty printer (shared implementation)
function _show_moduletemp_pretty(io::IO, mt::ModuleTemperature)
    println(io, "Module Temperature:")
    println(io, "──────────────────")
    println(io, "time     │ ", mt.time)
    println(io, "module temperature │ ", mt.module_temperature)
end

function Base.show(io::IO, mt::ModuleTemperature)
    _show_moduletemp_pretty(io, mt)
end

function Base.show(io::IO, ::MIME"text/plain", mt::ModuleTemperature)
    _show_moduletemp_pretty(io, mt)
end

# ModuleTemperature vector pretty printer (shared implementation)
function _show_moduletemp_vec_pretty(io::IO, v::Vector{<:ModuleTemperature})
    header = ["time", "module_temperature"]
    _show_table(
        io,
        "ModuleTemperature",
        header,
        i -> begin
            x = v[i]
            [string(x.time), string(x.module_temperature)]
        end,
        length(v),
    )
end

function Base.show(io::IO, v::Vector{<:ModuleTemperature})
    _show_moduletemp_vec_pretty(io, v)
end

function Base.show(io::IO, ::MIME"text/plain", v::Vector{<:ModuleTemperature})
    _show_moduletemp_vec_pretty(io, v)
end

# CellTemperature pretty printer (shared implementation)
function _show_celltemp_pretty(io::IO, ct::CellTemperature)
    println(io, "Cell Temperature:")
    println(io, "──────────────")
    println(io, "time   │ ", ct.time)
    print(io, "cell temperature │ ", ct.cell_temperature)
end

function Base.show(io::IO, ct::CellTemperature)
    _show_celltemp_pretty(io, ct)
end

function Base.show(io::IO, ::MIME"text/plain", ct::CellTemperature)
    _show_celltemp_pretty(io, ct)
end

# CellTemperature vector pretty printer (shared implementation)
function _show_celltemp_vec_pretty(io::IO, v::Vector{<:CellTemperature})
    header = ["time", "cell_temperature"]
    _show_table(
        io,
        "CellTemperature",
        header,
        i -> begin
            x = v[i]
            [string(x.time), string(x.cell_temperature)]
        end,
        length(v),
    )
end

function Base.show(io::IO, v::Vector{<:CellTemperature})
    _show_celltemp_vec_pretty(io, v)
end

function Base.show(io::IO, ::MIME"text/plain", v::Vector{<:CellTemperature})
    _show_celltemp_vec_pretty(io, v)
end
