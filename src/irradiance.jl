"""
    TotalIrradiance{T}

Total plane-of-array (POA) irradiance components for a single timestamp.

Typically computed from weather inputs and solar position using
`get_total_irradiance`.

# Fields:
* `time`: Timestamp
* `poa_global`: Total POA irradiance (W/m^2)
* `poa_direct`: POA direct component (W/m^2)
* `poa_diffuse`: POA diffuse component (W/m^2)
* `poa_sky_diffuse`: POA sky diffuse component (W/m^2)
* `poa_ground_diffuse`: POA ground-reflected diffuse component (W/m^2)
"""
struct TotalIrradiance{T}
    time::Any
    poa_global::T
    poa_direct::T
    poa_diffuse::T
    poa_sky_diffuse::T
    poa_ground_diffuse::T
end

function rrule(
    ::Type{TotalIrradiance},
    time,
    poa_global,
    poa_direct,
    poa_diffuse,
    poa_sky_diffuse,
    poa_ground_diffuse,
)
    y = TotalIrradiance{
        promote_type(
            typeof(poa_global),
            typeof(poa_direct),
            typeof(poa_diffuse),
            typeof(poa_sky_diffuse),
            typeof(poa_ground_diffuse),
        ),
    }(
        time,
        poa_global,
        poa_direct,
        poa_diffuse,
        poa_sky_diffuse,
        poa_ground_diffuse,
    )

    function pullback(ȳ)
        return (
            NoTangent(),
            NoTangent(),
            getproperty(ȳ, :poa_global),
            getproperty(ȳ, :poa_direct),
            getproperty(ȳ, :poa_diffuse),
            getproperty(ȳ, :poa_sky_diffuse),
            getproperty(ȳ, :poa_ground_diffuse),
        )
    end

    return y, pullback
end

const solar_constant = 1366.1  # W/m^2

"""
    get_projection(
        surface_tilt::Real,
        surface_azimuth::Real,
        solar_zenith::Real,
        solar_azimuth::Real
    ) -> Real

Compute the cosine projection (dot product) between a tilted surface normal and
the solar vector.

The returned value is clamped to `[-1, 1]`.

# Arguments:
* `surface_tilt`: Surface tilt angle (degrees)
* `surface_azimuth`: Surface azimuth angle (degrees)
* `solar_zenith`: Solar zenith angle (degrees)
* `solar_azimuth`: Solar azimuth angle (degrees)

# Returns:
* Scalar projection based on the angle of incidence (unitless)
"""
function get_projection(
    surface_tilt::Real,
    surface_azimuth::Real,
    solar_zenith::Real,
    solar_azimuth::Real,
)

    projection =
        cosd(surface_tilt) * cosd(solar_zenith) +
        sind(surface_tilt) * sind(solar_zenith) * cosd(solar_azimuth - surface_azimuth)
    projection = pv_smooth_clamp(projection, -one(projection), one(projection))

    return projection
end

"""
    haydavies(
        surface_tilt::Real,
        surface_azimuth::Real,
        dhi::Real,
        dni::Real,
        dni_extra::Real,
        solar_zenith::Union{Nothing, Real}=nothing,
        solar_azimuth::Union{Nothing, Real}=nothing,
    ) -> Real

    haydavies(
        surface_tilt::Real,
        surface_azimuth::Real,
        dhi::AbstractVector{<:Real},
        dni::AbstractVector{<:Real},
        dni_extra::AbstractVector{<:Real},
        solar_zenith::Union{Nothing, AbstractVector{<:Real}}=nothing,
        solar_azimuth::Union{Nothing, AbstractVector{<:Real}}=nothing,
    ) -> AbstractVector{<:Real}

    haydavies(
        surface_tilt::AbstractVector{<:Real},
        surface_azimuth::Real,
        dhi::AbstractVector{<:Real},
        dni::AbstractVector{<:Real},
        dni_extra::AbstractVector{<:Real},
        solar_zenith::Union{Nothing, AbstractVector{<:Real}}=nothing,
        solar_azimuth::Union{Nothing, AbstractVector{<:Real}}=nothing,
    ) -> AbstractVector{<:Real}

    haydavies(
        surface_tilt::Real,
        surface_azimuth::AbstractVector{<:Real},
        dhi::AbstractVector{<:Real},
        dni::AbstractVector{<:Real},
        dni_extra::AbstractVector{<:Real},
        solar_zenith::Union{Nothing, AbstractVector{<:Real}}=nothing,
        solar_azimuth::Union{Nothing, AbstractVector{<:Real}}=nothing,
    ) -> AbstractVector{<:Real}

    haydavies(
        surface_tilt::AbstractVector{<:Real},
        surface_azimuth::AbstractVector{<:Real},
        dhi::AbstractVector{<:Real},
        dni::AbstractVector{<:Real},
        dni_extra::AbstractVector{<:Real},
        solar_zenith::Union{Nothing, AbstractVector{<:Real}}=nothing,
        solar_azimuth::Union{Nothing, AbstractVector{<:Real}}=nothing,
    ) -> AbstractVector{<:Real}

Compute plane-of-array (POA) sky diffuse irradiance using the Hay & Davies model.

# Arguments:
* `surface_tilt`: Surface tilt angle(s) (degrees)
* `surface_azimuth`: Surface azimuth angle(s) (degrees)
* `dhi`: Diffuse horizontal irradiance (W/m^2)
* `dni`: Direct normal irradiance (W/m^2)
* `dni_extra`: Extraterrestrial direct normal irradiance (W/m^2)
* `solar_zenith`: Solar zenith angle(s) (degrees)
* `solar_azimuth`: Solar azimuth angle(s) (degrees)

# Returns:
* POA sky diffuse irradiance value(s) (W/m^2)
"""
function haydavies(
    surface_tilt::Real,
    surface_azimuth::Real,
    dhi::Real,
    dni::Real,
    dni_extra::Real,
    solar_zenith::Union{Nothing,Real} = nothing,
    solar_azimuth::Union{Nothing,Real} = nothing,
)

    # Calculate the ratio of tilted and horizontal beam irradiance
    projection = get_projection(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth)
    cos_tt = projection

    cos_tt = pv_smooth_max(cos_tt, zero(cos_tt))
    cos_solar_zenith = cosd(solar_zenith)
    Rb = cos_tt / pv_smooth_max(cos_solar_zenith, oftype(cos_solar_zenith, 0.01745))

    # Anisotropy Index
    AI = dni ./ dni_extra

    # Terms for the second term of equation 7
    term1 = 1 .- AI
    term2 = 0.5 * (1 + cosd(surface_tilt))

    poa_isotropic = pv_smooth_max(dhi * term1 * term2, zero(dhi))
    poa_circumsolar = pv_smooth_max(dhi * AI * Rb, zero(dhi))

    sky_diffuse = poa_isotropic .+ poa_circumsolar

    return sky_diffuse
end

function haydavies(
    surface_tilt::Real,
    surface_azimuth::Real,
    dhi::AbstractVector{<:Real},
    dni::AbstractVector{<:Real},
    dni_extra::AbstractVector{<:Real},
    solar_zenith::Union{Nothing,AbstractVector{<:Real}} = nothing,
    solar_azimuth::Union{Nothing,AbstractVector{<:Real}} = nothing,
)

    sky_diffuse = similar(dhi)

    for ii in eachindex(dhi)
        # Call the scalar version for each timestep
        sky_diffuse[ii] = haydavies(
            surface_tilt,
            surface_azimuth,
            dhi[ii],
            dni[ii],
            dni_extra[ii],
            solar_zenith[ii],
            solar_azimuth[ii],
        )
    end

    return sky_diffuse
end

function haydavies(
    surface_tilt::AbstractVector{<:Real},
    surface_azimuth::Real,
    dhi::AbstractVector{<:Real},
    dni::AbstractVector{<:Real},
    dni_extra::AbstractVector{<:Real},
    solar_zenith::Union{Nothing,AbstractVector{<:Real}} = nothing,
    solar_azimuth::Union{Nothing,AbstractVector{<:Real}} = nothing,
)

    sky_diffuse = similar(dhi)

    for ii in eachindex(dhi)
        # Call the scalar version for each timestep
        sky_diffuse[ii] = haydavies(
            surface_tilt[ii],
            surface_azimuth,
            dhi[ii],
            dni[ii],
            dni_extra[ii],
            solar_zenith[ii],
            solar_azimuth[ii],
        )
    end

    return sky_diffuse
end

function haydavies(
    surface_tilt::Real,
    surface_azimuth::AbstractVector{<:Real},
    dhi::AbstractVector{<:Real},
    dni::AbstractVector{<:Real},
    dni_extra::AbstractVector{<:Real},
    solar_zenith::Union{Nothing,AbstractVector{<:Real}} = nothing,
    solar_azimuth::Union{Nothing,AbstractVector{<:Real}} = nothing,
)

    sky_diffuse = similar(dhi)

    for ii in eachindex(dhi)
        # Call the scalar version for each timestep
        sky_diffuse[ii] = haydavies(
            surface_tilt,
            surface_azimuth[ii],
            dhi[ii],
            dni[ii],
            dni_extra[ii],
            solar_zenith[ii],
            solar_azimuth[ii],
        )
    end

    return sky_diffuse
end

function haydavies(
    surface_tilt::AbstractVector{<:Real},
    surface_azimuth::AbstractVector{<:Real},
    dhi::AbstractVector{<:Real},
    dni::AbstractVector{<:Real},
    dni_extra::AbstractVector{<:Real},
    solar_zenith::Union{Nothing,AbstractVector{<:Real}} = nothing,
    solar_azimuth::Union{Nothing,AbstractVector{<:Real}} = nothing,
)

    sky_diffuse = similar(dhi)

    for ii in eachindex(dhi)
        # Call the scalar version for each timestep
        sky_diffuse[ii] = haydavies(
            surface_tilt[ii],
            surface_azimuth[ii],
            dhi[ii],
            dni[ii],
            dni_extra[ii],
            solar_zenith[ii],
            solar_azimuth[ii],
        )
    end

    return sky_diffuse
end

"""
    get_extra_radiation(time::ZonedDateTime) -> Real
    get_extra_radiation(time::Vector{<:ZonedDateTime}) -> Vector{<:Real}

Compute extraterrestrial direct normal irradiance (DNI) using a day-of-year
correction to the solar constant.

# Arguments:
* `time`: Timestamp(s) as `ZonedDateTime`

# Returns:
* Extraterrestrial irradiance value(s) (W/m^2)
"""
function get_extra_radiation(time::ZonedDateTime)
    # Convert to day of year
    doy = dayofyear(time)

    B = (2 * π / 365) * (doy - 1)
    RoverR0sqrd = (
        1.00011 +
        0.034221 * cos(B) +
        0.00128 * sin(B) +
        0.000719 * cos(2 * B) +
        7.7e-05 * sin(2 * B)
    )

    Ea = solar_constant * RoverR0sqrd

    return Ea
end

function get_extra_radiation(time::Vector{<:ZonedDateTime})

    Ea = [get_extra_radiation(t) for t in time]

    return Ea
end

"""
    get_poa_ground_diffuse(
        ghi::Real,
        albedo::Real,
        surface_tilt::Real
    ) -> Real

    get_poa_ground_diffuse(
        ghi::Vector{<:Real},
        albedo::Real,
        surface_tilt::Real
    ) -> Vector{<:Real}

    get_poa_ground_diffuse(
        ghi::Vector{<:Real},
        albedo::Real,
        surface_tilt::AbstractVector{<:Real}
    ) -> Vector{<:Real}

Compute plane-of-array (POA) ground-reflected diffuse irradiance.

# Arguments:
* `ghi`: Global horizontal irradiance (W/m^2)
* `albedo`: Ground albedo (unitless, typically in `[0, 1]`)
* `surface_tilt`: Surface tilt angle(s) (degrees)

# Returns:
* POA ground diffuse irradiance value(s) (W/m^2)
"""
function get_poa_ground_diffuse(ghi::Real, albedo::Real, surface_tilt::Real)

    poa_ground_diffuse = ghi * albedo * (1 - cosd(surface_tilt)) * 0.5

    return poa_ground_diffuse
end

function get_poa_ground_diffuse(ghi::Vector{<:Real}, albedo::Real, surface_tilt::Real)

    poa_ground_diffuse =
        [get_poa_ground_diffuse(ghi[i], albedo, surface_tilt) for i in eachindex(ghi)]

    return poa_ground_diffuse
end

function get_poa_ground_diffuse(
    ghi::Vector{<:Real},
    albedo::Real,
    surface_tilt::AbstractVector{<:Real},
)

    poa_ground_diffuse =
        [get_poa_ground_diffuse(ghi[i], albedo, surface_tilt[i]) for i in eachindex(ghi)]

    return poa_ground_diffuse
end

"""
    get_angle_of_incidence(
        surface_tilt::Real,
        surface_azimuth::Real,
        solar_zenith::Real,
        solar_azimuth::Real
    ) -> Real

    get_angle_of_incidence(
        surface_tilt::Real,
        surface_azimuth::Real,
        solar_zenith::Vector{<:Real},
        solar_azimuth::Vector{<:Real}
    ) -> Vector{<:Real}

    get_angle_of_incidence(
        surface_tilt::Vector{<:Real},
        surface_azimuth::Real,
        solar_zenith::Vector{<:Real},
        solar_azimuth::Vector{<:Real}
    ) -> Vector{<:Real}

    get_angle_of_incidence(
        surface_tilt::Real,
        surface_azimuth::Vector{<:Real},
        solar_zenith::Vector{<:Real},
        solar_azimuth::Vector{<:Real}
    ) -> Vector{<:Real}

    get_angle_of_incidence(
        surface_tilt::Vector{<:Real},
        surface_azimuth::Vector{<:Real},
        solar_zenith::Vector{<:Real},
        solar_azimuth::Vector{<:Real}
    ) -> Vector{<:Real}

Compute the angle of incidence (AOI) between incoming solar rays and a tilted
surface.

# Arguments:
* `surface_tilt`: Surface tilt angle(s) (degrees)
* `surface_azimuth`: Surface azimuth angle(s) (degrees)
* `solar_zenith`: Solar zenith angle(s) (degrees)
* `solar_azimuth`: Solar azimuth angle(s) (degrees)

# Returns:
* Angle of incidence value(s) (degrees)
"""
function get_angle_of_incidence(
    surface_tilt::Real,
    surface_azimuth::Real,
    solar_zenith::Real,
    solar_azimuth::Real,
)

    projection = get_projection(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth)
    angle_of_incidence = acosd(projection)

    return angle_of_incidence
end

function get_angle_of_incidence(
    surface_tilt::Real,
    surface_azimuth::Real,
    solar_zenith::Vector{<:Real},
    solar_azimuth::Vector{<:Real},
)

    angle_of_incidence = [
        get_angle_of_incidence(
            surface_tilt,
            surface_azimuth,
            solar_zenith[i],
            solar_azimuth[i],
        ) for i in eachindex(solar_zenith)
    ]

    return angle_of_incidence
end

function get_angle_of_incidence(
    surface_tilt::Vector{<:Real},
    surface_azimuth::Real,
    solar_zenith::Vector{<:Real},
    solar_azimuth::Vector{<:Real},
)

    angle_of_incidence = [
        get_angle_of_incidence(
            surface_tilt[i],
            surface_azimuth,
            solar_zenith[i],
            solar_azimuth[i],
        ) for i in eachindex(solar_zenith)
    ]

    return angle_of_incidence
end

function get_angle_of_incidence(
    surface_tilt::Real,
    surface_azimuth::Vector{<:Real},
    solar_zenith::Vector{<:Real},
    solar_azimuth::Vector{<:Real},
)

    angle_of_incidence = [
        get_angle_of_incidence(
            surface_tilt,
            surface_azimuth[i],
            solar_zenith[i],
            solar_azimuth[i],
        ) for i in eachindex(solar_zenith)
    ]

    return angle_of_incidence
end

function get_angle_of_incidence(
    surface_tilt::Vector{<:Real},
    surface_azimuth::Vector{<:Real},
    solar_zenith::Vector{<:Real},
    solar_azimuth::Vector{<:Real},
)

    angle_of_incidence = [
        get_angle_of_incidence(
            surface_tilt[i],
            surface_azimuth[i],
            solar_zenith[i],
            solar_azimuth[i],
        ) for i in eachindex(solar_zenith)
    ]

    return angle_of_incidence
end

"""
    get_poa_direct(dni::Real, angle_of_incidence::Real) -> Real
    get_poa_direct(dni::Vector{<:Real}, angle_of_incidence::Vector{<:Real}) -> Vector{<:Real}

Compute plane-of-array (POA) direct irradiance from direct normal irradiance (DNI)
and angle of incidence (AOI).

The output is floored at `0` to avoid negative irradiance when the sun is behind
the plane.

# Arguments:
* `dni`: Direct normal irradiance (W/m^2)
* `angle_of_incidence`: Angle of incidence (degrees)

# Returns:
* POA direct irradiance value(s) (W/m^2)
"""
function get_poa_direct(dni::Real, angle_of_incidence::Real)

    poa_direct = dni * cosd(angle_of_incidence)
    poa_direct = pv_smooth_max(poa_direct, zero(poa_direct))

    return poa_direct
end

function get_poa_direct(dni::Vector{<:Real}, angle_of_incidence::Vector{<:Real})

    poa_direct = [get_poa_direct(dni[i], angle_of_incidence[i]) for i in eachindex(dni)]

    return poa_direct
end

"""
    get_total_irradiance(
        surface_tilt::Real,
        surface_azimuth::Real,
        weather_data::WeatherSample,
        solar_position::SolarPosition,
        albedo::Real,
        time::ZonedDateTime
    ) -> TotalIrradiance

    get_total_irradiance(
        surface_tilt::Real,
        surface_azimuth::Real,
        weather_data::WeatherSample,
        solar_position::SolarPosition,
        albedo::Real
    ) -> TotalIrradiance

    get_total_irradiance(
        surface_tilt::Real,
        surface_azimuth::Real,
        weather_data::Vector{<:WeatherSample},
        solar_position::Vector{<:SolarPosition},
        albedo::Real
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::AbstractVector{<:Real},
        surface_azimuth::Real,
        weather_data::Vector{<:WeatherSample},
        solar_position::Vector{<:SolarPosition},
        albedo::Real
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::Real,
        surface_azimuth::AbstractVector{<:Real},
        weather_data::Vector{<:WeatherSample},
        solar_position::Vector{<:SolarPosition},
        albedo::Real
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::AbstractVector{<:Real},
        surface_azimuth::AbstractVector{<:Real},
        weather_data::Vector{<:WeatherSample},
        solar_position::Vector{<:SolarPosition},
        albedo::Real
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::Real,
        surface_azimuth::Real,
        weather_data::Vector{<:WeatherSample},
        solar_position::Vector{<:SolarPosition},
        albedo::AbstractVector{<:Real}
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::AbstractVector{<:Real},
        surface_azimuth::Real,
        weather_data::Vector{<:WeatherSample},
        solar_position::Vector{<:SolarPosition},
        albedo::AbstractVector{<:Real}
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::Real,
        surface_azimuth::AbstractVector{<:Real},
        weather_data::Vector{<:WeatherSample},
        solar_position::Vector{<:SolarPosition},
        albedo::AbstractVector{<:Real}
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::AbstractVector{<:Real},
        surface_azimuth::AbstractVector{<:Real},
        weather_data::Vector{<:WeatherSample},
        solar_position::Vector{<:SolarPosition},
        albedo::AbstractVector{<:Real}
    ) -> Vector{<:TotalIrradiance}

        get_total_irradiance(
        surface_tilt::AbstractVector{<:Real},
        surface_azimuth::Real,
        weather_data::WeatherSample,
        solar_position::SolarPosition,
        albedo::Real,
        sim_time::AbstractVector{<:Real}
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::Real,
        surface_azimuth::AbstractVector{<:Real},
        weather_data::WeatherSample,
        solar_position::SolarPosition,
        albedo::Real,
        sim_time::AbstractVector{<:AReal}
    ) -> Vector{<:TotalIrradiance}

    get_total_irradiance(
        surface_tilt::AbstractVector{<:Real},
        surface_azimuth::AbstractVector{<:Real},
        weather_data::WeatherSample,
        solar_position::SolarPosition,
        albedo::Real,
        sim_time::AbstractVector{<:Real}
    ) -> Vector{<:TotalIrradiance}

Compute total plane-of-array (POA) irradiance components using weather inputs,
solar position, the Hay & Davies sky diffuse model, and a ground-reflection term.

Returns a `TotalIrradiance` (or vector of `TotalIrradiance`) with fields:
`time`, `poa_global`, `poa_direct`, `poa_diffuse`, `poa_sky_diffuse`,
`poa_ground_diffuse`.

# Arguments:
* `surface_tilt`: Surface tilt angle(s) (degrees)
* `surface_azimuth`: Surface azimuth angle(s) (degrees)
* `weather_data`: Weather input(s) containing at least `time`, `ghi`, `dhi`, `dni`
* `solar_position`: Solar position input(s) containing at least `apparent_zenith`, `azimuth`
* `albedo`: Ground albedo (unitless)
* `time`: Optional explicit timestamp for the returned `TotalIrradiance`. If omitted, `weather_data.time` is used.
* `sim_time`: Simulation time offsets in seconds from `weather_data.time`. Required for overloads where a single `WeatherSample` and `SolarPosition` are reused across multiple simulated times.

# Returns:
* Total POA irradiance components (W/m^2)

```jldoctest
julia> using TimeZones

julia> weather_data = WeatherSample(
           ZonedDateTime(2020, 6, 1, 12, 0, 0, tz"UTC"),
           800.0,       # ghi
           100.0,       # dhi
           900.0,       # dni
           25.0,        # temp_air
           50.0,        # relative_humidity
           20.0,        # temp_dewpoint
           101325.0,    # pressure
           5.0,         # wind_speed
           180.0,       # wind_direction
           0.1,         # albedo
       );

julia> solar_pos = SolarPosition(
           ZonedDateTime(2020, 6, 1, 12, 0, 0, tz"UTC"),
           30.0,   # apparent_zenith
           30.0,   # zenith
           60.0,   # apparent_elevation
           60.0,   # elevation
           180.0,  # azimuth
           1.0,    # equation_of_time
       );

julia> out = get_total_irradiance(30.0, 180.0, weather_data, solar_pos, 0.2);

julia> fieldnames(typeof(out)) == (:time,:poa_global,:poa_direct,:poa_diffuse,:poa_sky_diffuse,:poa_ground_diffuse)
true

julia> out.poa_global ≈ out.poa_direct + out.poa_diffuse
true

julia> out.poa_diffuse ≈ out.poa_sky_diffuse + out.poa_ground_diffuse
true

julia> isapprox(out.poa_global, 962.9920322637975; atol=1e-6)
true
```
"""
function get_total_irradiance(
    surface_tilt::Real,
    surface_azimuth::Real,
    weather_data::WeatherSample,
    solar_position::SolarPosition,
    albedo::Real,
    time::ZonedDateTime,
)
    dni_extra = get_extra_radiation(time)
    angle_of_incidence = get_angle_of_incidence(
        surface_tilt,
        surface_azimuth,
        solar_position.apparent_zenith,
        solar_position.azimuth,
    )

    poa_sky_diffuse = haydavies(
        surface_tilt,
        surface_azimuth,
        weather_data.dhi,
        weather_data.dni,
        dni_extra,
        solar_position.apparent_zenith,
        solar_position.azimuth,
    )
    poa_ground_diffuse = get_poa_ground_diffuse(weather_data.ghi, albedo, surface_tilt)
    poa_diffuse = poa_sky_diffuse + poa_ground_diffuse
    poa_direct = get_poa_direct(weather_data.dni, angle_of_incidence)
    poa_global = poa_direct + poa_diffuse

    return TotalIrradiance(
        time,
        poa_global,
        poa_direct,
        poa_diffuse,
        poa_sky_diffuse,
        poa_ground_diffuse,
    )
end

function get_total_irradiance(
    surface_tilt::Real,
    surface_azimuth::Real,
    weather_data::WeatherSample,
    solar_position::SolarPosition,
    albedo::Real,
)

    return get_total_irradiance(
        surface_tilt,
        surface_azimuth,
        weather_data,
        solar_position,
        albedo,
        weather_data.time,
    )
end

function get_total_irradiance(
    surface_tilt::Real,
    surface_azimuth::Real,
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
    albedo::Real,
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt,
            surface_azimuth,
            weather_data[ii],
            solar_position[ii],
            albedo,
            weather_data[ii].time,
        )
    end
end

function get_total_irradiance(
    surface_tilt::AbstractVector{<:Real},
    surface_azimuth::Real,
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
    albedo::Real,
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    length(surface_tilt) == n || throw(
        DimensionMismatch(
            "surface_tilt has length $(length(surface_tilt)) but weather_data has $n rows",
        ),
    )

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt[ii],
            surface_azimuth,
            weather_data[ii],
            solar_position[ii],
            albedo,
            weather_data[ii].time,
        )
    end
end

function get_total_irradiance(
    surface_tilt::Real,
    surface_azimuth::AbstractVector{<:Real},
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
    albedo::Real,
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    length(surface_azimuth) == n || throw(
        DimensionMismatch(
            "surface_azimuth has length $(length(surface_azimuth)) but weather_data has $n rows",
        ),
    )

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt,
            surface_azimuth[ii],
            weather_data[ii],
            solar_position[ii],
            albedo,
            weather_data[ii].time,
        )
    end
end

function get_total_irradiance(
    surface_tilt::AbstractVector{<:Real},
    surface_azimuth::AbstractVector{<:Real},
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
    albedo::Real,
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    length(surface_tilt) == n || throw(
        DimensionMismatch(
            "surface_tilt has length $(length(surface_tilt)) but weather_data has $n rows",
        ),
    )

    length(surface_azimuth) == n || throw(
        DimensionMismatch(
            "surface_azimuth has length $(length(surface_azimuth)) but weather_data has $n rows",
        ),
    )

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt[ii],
            surface_azimuth[ii],
            weather_data[ii],
            solar_position[ii],
            albedo,
            weather_data[ii].time,
        )
    end
end

function get_total_irradiance(
    surface_tilt::Real,
    surface_azimuth::Real,
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
    albedo::AbstractVector{<:Real},
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    length(albedo) == n || throw(
        DimensionMismatch(
            "albedo has length $(length(albedo)) but weather_data has $n rows",
        ),
    )

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt,
            surface_azimuth,
            weather_data[ii],
            solar_position[ii],
            albedo[ii],
            weather_data[ii].time,
        )
    end
end

function get_total_irradiance(
    surface_tilt::AbstractVector{<:Real},
    surface_azimuth::Real,
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
    albedo::AbstractVector{<:Real},
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    length(surface_tilt) == n || throw(
        DimensionMismatch(
            "surface_tilt has length $(length(surface_tilt)) but weather_data has $n rows",
        ),
    )

    length(albedo) == n || throw(
        DimensionMismatch(
            "albedo has length $(length(albedo)) but weather_data has $n rows",
        ),
    )

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt[ii],
            surface_azimuth,
            weather_data[ii],
            solar_position[ii],
            albedo[ii],
            weather_data[ii].time,
        )
    end
end

function get_total_irradiance(
    surface_tilt::Real,
    surface_azimuth::AbstractVector{<:Real},
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
    albedo::AbstractVector{<:Real},
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    length(surface_azimuth) == n || throw(
        DimensionMismatch(
            "surface_azimuth has length $(length(surface_azimuth)) but weather_data has $n rows",
        ),
    )

    length(albedo) == n || throw(
        DimensionMismatch(
            "albedo has length $(length(albedo)) but weather_data has $n rows",
        ),
    )

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt,
            surface_azimuth[ii],
            weather_data[ii],
            solar_position[ii],
            albedo[ii],
            weather_data[ii].time,
        )
    end
end

function get_total_irradiance(
    surface_tilt::AbstractVector{<:Real},
    surface_azimuth::AbstractVector{<:Real},
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
    albedo::AbstractVector{<:Real},
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    length(surface_tilt) == n || throw(
        DimensionMismatch(
            "surface_tilt has length $(length(surface_tilt)) but weather_data has $n rows",
        ),
    )

    length(surface_azimuth) == n || throw(
        DimensionMismatch(
            "surface_azimuth has length $(length(surface_azimuth)) but weather_data has $n rows",
        ),
    )

    length(albedo) == n || throw(
        DimensionMismatch(
            "albedo has length $(length(albedo)) but weather_data has $n rows",
        ),
    )

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt[ii],
            surface_azimuth[ii],
            weather_data[ii],
            solar_position[ii],
            albedo[ii],
            weather_data[ii].time,
        )
    end
end

function get_total_irradiance(
    surface_tilt::AbstractVector{<:Real},
    surface_azimuth::Real,
    weather_data::WeatherSample,
    solar_position::SolarPosition,
    albedo::Real,
    sim_time::AbstractVector{<:Real},
)

    n = length(surface_tilt)

    length(sim_time) == n ||
        throw(DimensionMismatch("sim_time has length $(length(sim_time)) but expected $n"))

    time_zdt = ignore_derivatives() do
        _sim_to_zoneddt(weather_data.time, sim_time)
    end

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt[ii],
            surface_azimuth,
            weather_data,
            solar_position,
            albedo,
            time_zdt[ii],
        )
    end
end

function get_total_irradiance(
    surface_tilt::Real,
    surface_azimuth::AbstractVector{<:Real},
    weather_data::WeatherSample,
    solar_position::SolarPosition,
    albedo::Real,
    sim_time::AbstractVector{<:Real},
)

    n = length(surface_azimuth)

    length(sim_time) == n ||
        throw(DimensionMismatch("sim_time has length $(length(sim_time)) but expected $n"))

    time_zdt = ignore_derivatives() do
        _sim_to_zoneddt(weather_data.time, sim_time)
    end

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt,
            surface_azimuth[ii],
            weather_data,
            solar_position,
            albedo,
            time_zdt[ii],
        )
    end
end

function get_total_irradiance(
    surface_tilt::AbstractVector{<:Real},
    surface_azimuth::AbstractVector{<:Real},
    weather_data::WeatherSample,
    solar_position::SolarPosition,
    albedo::Real,
    sim_time::AbstractVector{<:Real},
)

    n = length(surface_tilt)

    length(surface_azimuth) == n || throw(
        DimensionMismatch(
            "surface_azimuth has length $(length(surface_azimuth)) but expected $n",
        ),
    )
    length(sim_time) == n ||
        throw(DimensionMismatch("sim_time has length $(length(sim_time)) but expected $n"))

    time_zdt = ignore_derivatives() do
        _sim_to_zoneddt(weather_data.time, sim_time)
    end

    return map(1:n) do ii
        get_total_irradiance(
            surface_tilt[ii],
            surface_azimuth[ii],
            weather_data,
            solar_position,
            albedo,
            time_zdt[ii],
        )
    end
end

# TotalIrradiance pretty printer (shared implementation)
function _show_total_irradiance_pretty(io::IO, ti::TotalIrradiance)
    println(io, "TotalIrradiance")
    println(io, "────────────")
    println(io, "time               │ ", ti.time)
    println(io, "poa_global         │ ", ti.poa_global)
    println(io, "poa_direct         │ ", ti.poa_direct)
    println(io, "poa_diffuse        │ ", ti.poa_diffuse)
    println(io, "poa_sky_diffuse    │ ", ti.poa_sky_diffuse)
    println(io, "poa_ground_diffuse │ ", ti.poa_ground_diffuse)
end

function Base.show(io::IO, ti::TotalIrradiance)
    _show_total_irradiance_pretty(io, ti)
end

function Base.show(io::IO, ::MIME"text/plain", ti::TotalIrradiance)
    _show_total_irradiance_pretty(io, ti)
end

# TotalIrradiance vector pretty printer (shared implementation)
function show_total_irradiance_vec_pretty(io::IO, v::Vector{<:TotalIrradiance})
    header = [
        "time",
        "poa_global",
        "poa_direct",
        "poa_diffuse",
        "poa_sky_diffuse",
        "poa_ground_diffuse",
    ]

    _show_table(
        io,
        "TotalIrradiance",
        header,
        i -> begin
            x = v[i]
            [
                string(x.time),
                string(x.poa_global),
                string(x.poa_direct),
                string(x.poa_diffuse),
                string(x.poa_sky_diffuse),
                string(x.poa_ground_diffuse),
            ]
        end,
        length(v),
    )
end

function Base.show(io::IO, v::Vector{<:TotalIrradiance})
    show_total_irradiance_vec_pretty(io, v)
end

function Base.show(io::IO, ::MIME"text/plain", v::Vector{<:TotalIrradiance})
    show_total_irradiance_vec_pretty(io, v)
end
