# Function to convert simulation time to zoned date time
function _sim_to_zoneddt(t0::ZonedDateTime, sim_time::AbstractVector{<:Real})
    # If any timestep is < 1 ms, rounding to milliseconds can collapse distinct times.
    if length(sim_time) ≥ 2
        dtmin = minimum(abs.(diff(sim_time)))
        if dtmin > 0 && dtmin < 1e-3
            @warn "sim_time contains timesteps smaller than 1 ms; conversion rounds to milliseconds and may lose time resolution." dtmin
        end
    end

    zdt = t0 .+ Millisecond.(round.(Int, 1000 .* sim_time))

    return zdt
end

"""
    panel_tilt_azimuth(
        motion::AbstractMatrix{<:Real},
        install_tilt_deg::Real = 0.0,
        install_azimuth_deg::Real = 0.0,
        yaw_offset_deg::Real = 0.0
    ) -> (tilt, azimuth)

Compute the panel tilt and azimuth time series from a 6-DOF motion array.

This function interprets columns 4–6 of `motion` as roll, pitch, and yaw,
respectively, and applies `yaw_offset_deg` as an offset to the yaw channel. The
panel normal in the world frame is computed by rotating the installed panel
normal by the body motion. Tilt is computed from the vertical component of the
normal, and azimuth is computed from the horizontal components.

Azimuth is returned in degrees, clockwise from North (ENU), mapped to ([0, 360]).

# Arguments:
* `motion`: Array with roll, pitch, yaw in columns 4–6 (radians), interpreted in ZYX Euler order
* `install_tilt_deg`: Initial panel tilt angle relative to the body frame (degrees)
* `install_azimuth_deg`: Initial panel azimuth angle relative to the body frame (degrees)
* `yaw_offset_deg`: Initial yaw offset added to the yaw channel (degrees)

# Returns:
* `tilt`: Vector of panel tilt angles (degrees), in ([0, 180])
* `azimuth`: Vector of panel azimuth angles (degrees), in ([0, 360)), clockwise from North (ENU)

```jldoctest
julia> motion = 0.2*ones(1, 6);

julia> tilt, azimuth = panel_tilt_azimuth(motion);

julia> isapprox(tilt[1], 16.15129532916803; atol=1e-6)
true

julia> isapprox(azimuth[1], 124.11762388541732; atol=1e-6)
true
```
"""
function panel_tilt_azimuth(
    motion::AbstractMatrix{<:Real},
    install_tilt_deg::Real = 0.0,
    install_azimuth_deg::Real = 0.0,
    yaw_offset_deg::Real = 0.0,
)
    roll = motion[:, 4]
    pitch = motion[:, 5]
    heading_offset_deg = 90 # wave heading offset from ENU North
    yaw = motion[:, 6] .+ deg2rad(yaw_offset_deg)

    # Panel normal in body coordinates from installed tilt/azimuth
    β = deg2rad(install_tilt_deg)
    γ = deg2rad(install_azimuth_deg)

    # ENU convention: x = East, y = North, z = Up
    # azimuth clockwise from North
    n0x = sin(β) * sin(γ)
    n0y = sin(β) * cos(γ)
    n0z = cos(β)

    # Apply body rotation Rz(yaw) * Ry(pitch) * Rx(roll)
    nx =
        (cos.(yaw) .* cos.(pitch)) .* n0x .+
        (cos.(yaw) .* sin.(pitch) .* sin.(roll) .- sin.(yaw) .* cos.(roll)) .* n0y .+
        (cos.(yaw) .* sin.(pitch) .* cos.(roll) .+ sin.(yaw) .* sin.(roll)) .* n0z

    ny =
        (sin.(yaw) .* cos.(pitch)) .* n0x .+
        (sin.(yaw) .* sin.(pitch) .* sin.(roll) .+ cos.(yaw) .* cos.(roll)) .* n0y .+
        (sin.(yaw) .* sin.(pitch) .* cos.(roll) .- cos.(yaw) .* sin.(roll)) .* n0z

    nz =
        (-sin.(pitch)) .* n0x .+ (cos.(pitch) .* sin.(roll)) .* n0y .+
        (cos.(pitch) .* cos.(roll)) .* n0z

    # Unsigned tilt from horizontal magnitude and |nz|
    nh = sqrt.(nx .^ 2 .+ ny .^ 2)
    tilt = atan.(nh, pv_smooth_abs.(nz))

    azimuth_raw = atan.(nx, ny)
    azimuth = mod.(azimuth_raw, 2π)

    return rad2deg.(tilt), rad2deg.(azimuth), nz
end

"""
    get_ocean_surface_albedo(
        global_horizontal_irradiance::Real,
        solar_zenith_angle::Real,
        day_of_year::Real
    ) -> Real

    get_ocean_surface_albedo(
        global_horizontal_irradiance::AbstractVector{<:Real},
        solar_zenith_angle::AbstractVector{<:Real},
        day_of_year::AbstractVector{<:Real}
    ) -> Vector{Float64}

    get_ocean_surface_albedo(
        weather_data,
        solar_position
    ) -> Vector{Float64}

Compute broadband ocean surface albedo using the parameterization of Huang et al. (2024).

The model computes ocean surface albedo as a function of global horizontal irradiance,
solar zenith angle, and day of year. Atmospheric transparency is estimated as the ratio
of downward shortwave irradiance at the surface to top-of-atmosphere horizontal irradiance.

# Arguments:
* `global_horizontal_irradiance`: Global horizontal irradiance, (Q_d), in W/m^2
* `solar_zenith_angle`: Solar zenith angle in degrees
* `day_of_year`: Day of year
* `weather_data`: Weather input(s) containing at least `ghi` (W/m^2) and `time`
* `solar_position`: Solar position input(s) containing at least `zenith` (degrees)

# Returns:
* Ocean surface albedo(s), bounded smoothly to ([0,1])

```jldoctest
julia> albedo = get_ocean_surface_albedo(800.0, 30.0, 180.0);

julia> isapprox(albedo, 0.033189; atol=1e-3)
true
```

# Reference:
* Huang, C. J., Wang, G., Chen, S., Guo, J., & Qiao, F. (2024).
  *An effective parameterization of broadband ocean surface albedo applicable to all skies*.
  Ocean Modelling, 190, 102394.
  https://doi.org/10.1016/j.ocemod.2024.102394
"""
function get_ocean_surface_albedo(
    global_horizontal_irradiance::Real,
    solar_zenith_angle::Real,
    day_of_year::Real,
)
    S0 = 1361 # W/m^2, solar constant
    Qd = global_horizontal_irradiance

    mu = cosd(solar_zenith_angle)
    Qtop = S0 * mu / (1 + 0.0167 * sin(2 * pi * (day_of_year - 93.5) / 365))^2
    ϵ = 1.0e-6
    Qtop_safe = sqrt(Qtop^2 + ϵ^2)
    Beta = Qd / Qtop_safe

    a = 0.79 * exp(-4.3 * mu) - 0.06
    b = -0.096 * exp(-4.3 * mu) + 0.06
    albedo = a * Beta + b

    return pv_smooth_clamp(albedo, 0, 1)
end

function get_ocean_surface_albedo(
    weather_data::WeatherSample,
    solar_position::SolarPosition,
)

    global_horizontal_irradiance = weather_data.ghi
    solar_zenith_angle = solar_position.zenith
    day_of_year = dayofyear(weather_data.time)

    albedo = get_ocean_surface_albedo(
        global_horizontal_irradiance,
        solar_zenith_angle,
        day_of_year,
    )
    return albedo
end

function get_ocean_surface_albedo(
    weather_data::Vector{<:WeatherSample},
    solar_position::Vector{<:SolarPosition},
)

    n = length(weather_data)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but weather_data has $n rows",
        ),
    )

    albedo = [
        get_ocean_surface_albedo(weather_data[ii], solar_position[ii]) for
        ii in eachindex(weather_data)
    ]
    return albedo
end

"""
    rolling_average_sapm_cell_temperature(
        total_irradiance::AbstractVector{<:TotalIrradiance},
        weather_data::WeatherSample;
        time_window::Real=5 * 60.0,
        a::Real=-3.47,
        b::Real=-0.0594,
        deltaT::Real=3.0,
        irrad_ref::Real=1000.0
    ) -> Vector{CellTemperature}

Compute SAPM cell temperature and apply a trailing rolling average over a specified time window.

The model first computes cell temperature using [`sapm_cell_temperature`](@ref), then applies
a trailing rolling mean over the previous `time_window` seconds, including the current sample.
The time step is assumed to be uniform and is inferred from the first two timestamps.

# Arguments:
* `total_irradiance`: Total irradiance input(s) containing at least `time`
* `weather_data`: Weather input(s) required by [`sapm_cell_temperature`](@ref)
* `time_window`: Width of the trailing rolling-average window in seconds
* `a`, `b`: Empirical coefficients for the module temperature model
* `deltaT`: Cell-to-module temperature difference at reference irradiance (°C)
* `irrad_ref`: Reference irradiance used for scaling `deltaT` (W/m^2)

# Returns:
* Vector of `CellTemperature` values with timestamps taken from `total_irradiance`
  and temperatures equal to the trailing rolling average of SAPM cell temperature

```jldoctest
julia> using TimeZones

julia> total_irradiance = [
            TotalIrradiance(
            ZonedDateTime(2023, 1, 1, tz"America/Denver"),
                1000.0, 0.0, 0.0, 0.0, 0.0
            )
            TotalIrradiance(
                ZonedDateTime(2023, 1, 2, tz"America/Denver"),
                1000.0, 0.0, 0.0, 0.0, 0.0
            )
        ];

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

julia> cell_temperature = rolling_average_sapm_cell_temperature(total_irradiance, weather_data);

julia> isapprox(cell_temperature[1].cell_temperature,  57.3225; atol=1e-3)
true
```
"""
function rolling_average_sapm_cell_temperature(
    total_irradiance::AbstractVector{<:TotalIrradiance},
    weather_data::WeatherSample;
    time_window::Real = 5 * 60.0,  # 5 minute rolling average
    a::Real = -3.47,
    b::Real = -0.0594,
    deltaT::Real = 3.0,
    irrad_ref::Real = 1000.0,
)
    cell_temperature = sapm_cell_temperature(
        total_irradiance,
        weather_data;
        a = a,
        b = b,
        deltaT = deltaT,
        irrad_ref = irrad_ref,
    )

    temps = getfield.(cell_temperature, :cell_temperature)
    times = getfield.(cell_temperature, :time)
    dt = value(times[2] - times[1]) / 1000 # uniform time steps

    window_n = max(1, Int(round(time_window / dt)))
    half_window = fld(window_n - 1, 2)

    cell_temperature_rolling_avg = [
        CellTemperature(
            total_irradiance[i].time,
            mean(temps[max(1, i-half_window):min(length(temps), i+half_window)]),
        ) for i in eachindex(temps)
    ]

    return cell_temperature_rolling_avg
end

# Define a simple panel geometry including width and height
struct Panel
    width::Float64
    height::Float64
end

# Define a simple axis-aligned box obstacle with coordinates relative to panel local frame
struct BoxObstacle
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    zmin::Float64
    zmax::Float64
end

# Return a smooth surrogate in [0, 1] for whether a ray intersects the box.
function _ray_intersects_box(
    origin::NTuple{3,<:Real},
    dir::NTuple{3,<:Real},
    box::BoxObstacle;
    height_frac::Real = 0.5,
    xy_hardness::Real = 100.0,
    z_hardness::Real = 100.0,
    eps_z::Real = 1e-6,
)
    zc = box.zmin + height_frac * (box.zmax - box.zmin)

    t = zc / (dir[3] + eps_z)
    xh = origin[1] + t * dir[1]
    yh = origin[2] + t * dir[2]

    sx1 = PVlib.pv_smooth_step(xh - box.xmin; hardness = xy_hardness)
    sx2 = PVlib.pv_smooth_step(box.xmax - xh; hardness = xy_hardness)
    sy1 = PVlib.pv_smooth_step(yh - box.ymin; hardness = xy_hardness)
    sy2 = PVlib.pv_smooth_step(box.ymax - yh; hardness = xy_hardness)

    in_xy = sx1 * sx2 * sy1 * sy2
    forward = PVlib.pv_smooth_step(t; hardness = z_hardness)

    return forward * in_xy
end

# Find the sun vector in world coordinates given the solar zenith and azimuth angles
function _sun_vector_world(zenith_deg::Real, azimuth_deg::Real)
    θ = deg2rad(zenith_deg)
    ϕ = deg2rad(azimuth_deg)

    sx = sin(θ) * sin(ϕ)
    sy = sin(θ) * cos(ϕ)
    sz = cos(θ)

    return (Float64(sx), Float64(sy), Float64(sz))
end

# Transform a vector from world coordinates to panel coordinates given the panel tilt and azimuth angles
function _world_to_panel(s::NTuple{3,<:Real}, panel_tilt_deg::Real, panel_azimuth_deg::Real)
    β = deg2rad(panel_tilt_deg)
    γ = deg2rad(panel_azimuth_deg)

    # panel basis vectors in world coordinates
    x̂ = (cos(γ), -sin(γ), 0.0)
    ŷ = (cos(β) * sin(γ), cos(β) * cos(γ), -sin(β))
    ẑ = (sin(β) * sin(γ), sin(β) * cos(γ), cos(β))

    sx = dot(s, x̂)
    sy = dot(s, ŷ)
    sz = dot(s, ẑ)

    return (sx, sy, sz)
end

# Get the sun vector in panel coordinates given the solar position and panel orientation
function _sun_vector_panel(
    solar_position::SolarPosition,
    panel_tilt_deg::Real,
    panel_azimuth_deg::Real,
)
    sun_vec_world = _sun_vector_world(solar_position.zenith, solar_position.azimuth)
    sun_vec_panel = _world_to_panel(sun_vec_world, panel_tilt_deg, panel_azimuth_deg)
    return sun_vec_panel
end

"""
    get_shaded_fraction(
        solar_position::SolarPosition,
        panel_tilt_deg::Real,
        panel_azimuth_deg::Real,
        panel::Panel,
        obstacle::BoxObstacle;
        nx::Int=40,
        ny::Int=40
    ) -> Float64

    get_shaded_fraction(
        solar_position::SolarPosition,
        panel_tilt_deg::AbstractVector{<:Real},
        panel_azimuth_deg::AbstractVector{<:Real},
        panel::Panel,
        obstacle::BoxObstacle;
        nx::Int=40,
        ny::Int=40
    ) -> Vector{Float64}

Estimate the fraction of panel area shaded by an axis-aligned box obstacle using grid-based ray casting.

The model computes the sun vector in panel coordinates from the solar position and panel
orientation, then samples a uniform `nx x ny` grid of points across the panel surface.
For each point, a ray is cast toward the sun and tested for intersection with the obstacle.
The shaded fraction is estimated as the fraction of sampled points whose rays intersect
the obstacle.

The panel is assumed to lie in its local coordinate system in the plane `z = 0`,
centered at the origin. The obstacle is assumed to be defined in the same panel-local
coordinate system.

# Arguments:
* `solar_position`: Solar position input containing at least `zenith` and `azimuth` in degrees
* `panel_tilt_deg`: Panel tilt angle(s) in degrees
* `panel_azimuth_deg`: Panel azimuth angle(s) in degrees
* `panel`: Panel geometry containing at least `width` and `height`
* `obstacle`: Axis-aligned box obstacle in panel-local coordinates
* `nx`: Number of sample points along the panel width
* `ny`: Number of sample points along the panel height

# Returns:
* Estimated shaded fraction(s) of the panel area, bounded to `[0, 1]`

```jldoctest
julia> using TimeZones

julia> panel = Panel(2.0, 1.0);

julia> obstacle = BoxObstacle(-1.0, 1.0, -0.1, 0.1, 0.3, 0.4);

julia> solar_pos = SolarPosition(
           ZonedDateTime(2020, 6, 1, 12, 0, 0, tz"UTC"),
           30.0,   # apparent_zenith
           30.0,   # zenith
           60.0,   # apparent_elevation
           60.0,   # elevation
           180.0,  # azimuth
           1.0,    # equation_of_time
       );

julia> fshade = get_shaded_fraction(solar_pos, 20.0, 180.0, panel, obstacle);

julia> isapprox(fshade, 0.199933; atol=1e-3)
true
```
"""
function get_shaded_fraction(
    solar_position::SolarPosition,
    panel_tilt_deg::Real,
    panel_azimuth_deg::Real,
    panel::Panel,
    obstacle::BoxObstacle;
    nx::Int = 40,
    ny::Int = 40,
)
    sun_vec = _sun_vector_panel(solar_position, panel_tilt_deg, panel_azimuth_deg)

    shaded_sum = 0
    n_total = nx * ny

    for i = 1:nx, j = 1:ny
        x = -panel.width / 2 + (i - 0.5) * panel.width / nx
        y = -panel.height / 2 + (j - 0.5) * panel.height / ny
        p = (x, y, 0.0)

        shaded = _ray_intersects_box(p, sun_vec, obstacle)
        shaded_sum += shaded
    end

    return shaded_sum / n_total
end

function get_shaded_fraction(
    solar_position::SolarPosition,
    panel_tilt_deg::AbstractVector{<:Real},
    panel_azimuth_deg::AbstractVector{<:Real},
    panel::Panel,
    obstacle::BoxObstacle;
    nx::Int = 40,
    ny::Int = 40,
)
    return [
        get_shaded_fraction(
            solar_position,
            panel_tilt_deg[ii],
            panel_azimuth_deg[ii],
            panel,
            obstacle;
            nx = nx,
            ny = ny,
        ) for ii in eachindex(panel_tilt_deg)
    ]
end

"""
    get_power_norm(
        total_irradiance::TotalIrradiance,
        shaded_fraction::Real,
        n_cells_per_column::Int
    ) -> Real

    get_power_norm(
        total_irradiance::AbstractVector{<:TotalIrradiance},
        shaded_fraction::AbstractVector{<:Real},
        n_cells_per_column::Int
    ) -> Vector{Float64}

Compute a normalized power factor for partially shaded modules using the diffuse fraction
of plane-of-array irradiance and the shaded fraction of the module area.

The model assumes that power decreases approximately linearly with shaded fraction until the
shaded fraction reaches `1 / n_cells_per_column`, after which the normalized power approaches
the diffuse fraction of plane-of-array irradiance. A smooth transition between these two
regimes is applied using `PVlib.pv_smooth_step`.

The diffuse fraction is computed as:
`poa_diffuse / poa_global`

# Arguments:
* `total_irradiance`: Total irradiance input(s) containing at least `poa_global` and `poa_diffuse`
* `shaded_fraction`: Shaded fraction(s) of the module area, typically in `[0, 1]`
* `n_cells_per_column`: Number of cells per column used to define the critical shaded fraction

# Returns:
* Normalized power factor(s), typically in `[0, 1]`, that can be applied to module or array DC power

```jldoctest
julia> using TimeZones

julia> total_irradiance = TotalIrradiance(
           ZonedDateTime(2023, 1, 1, tz"America/Denver"),
           1000.0, 800.0, 200.0, 200.0, 0.0
       );

julia> power_norm = get_power_norm(total_irradiance, 0.0, 9);

julia> isapprox(power_norm, 1.0; atol=1e-4)
true
```
"""
function get_power_norm(
    total_irradiance::TotalIrradiance,
    shaded_fraction::Real,
    n_cells_per_column::Int,
)

    diffuse_fraction = total_irradiance.poa_diffuse / total_irradiance.poa_global
    critical_shaded_fraction = 1 / n_cells_per_column

    p1 = 1 - (1 - diffuse_fraction) * shaded_fraction * n_cells_per_column
    p2 = diffuse_fraction
    smooth_weighting = PVlib.pv_smooth_step(critical_shaded_fraction - shaded_fraction)

    power_norm = smooth_weighting * p1 + (1 - smooth_weighting) * p2
    return power_norm
end

function get_power_norm(
    total_irradiance::AbstractVector{<:TotalIrradiance},
    shaded_fraction::AbstractVector{<:Real},
    n_cells_per_column::Int,
)

    n = length(total_irradiance)

    length(shaded_fraction) == n || throw(
        DimensionMismatch(
            "shaded_fraction has $(length(shaded_fraction)) elements but total_irradiance has $n elements",
        ),
    )

    return [
        get_power_norm(total_irradiance[ii], shaded_fraction[ii], n_cells_per_column) for
        ii = 1:n
    ]
end

"""
    sapm_dc_components_shaded(
        pv_module::SolarModule,
        effective_irradiance::EffectiveIrradiance,
        cell_temperature::CellTemperature,
        power_norm::Real,
        temperature_ref::Real=25.0,
        irradiance_ref::Real=1000.0,
        q::Real=1.602176634e-19,
        kb::Real=1.380649e-23
    ) -> Vector{DCComponents}

    sapm_dc_components_shaded(
        pv_module::SolarModule,
        effective_irradiance::Vector{<:EffectiveIrradiance},
        cell_temperature::Vector{<:CellTemperature},
        power_norm::Vector{<:Real},
        temperature_ref::Real=25.0,
        irradiance_ref::Real=1000.0,
        q::Real=1.602176634e-19,
        kb::Real=1.380649e-23
    ) -> Vector{DCComponents}

Compute SAPM DC output components with an additional shading-based power normalization
applied to the maximum power, `p_mp`.

This function first evaluates `sapm_dc_components` using the provided module,
effective irradiance, and cell temperature inputs. It then scales only the maximum power
component, `p_mp`, by `power_norm`, while leaving the other SAPM DC quantities unchanged.

This is useful when SAPM electrical performance is computed normally, but an externally
calculated shading loss factor should be applied to the module power output.

# Arguments:
* `pv_module`: Sandia module parameters
* `effective_irradiance`: Effective irradiance (typically from `sapm_effective_irradiance`)
* `cell_temperature`: Cell temperature (°C)
* `power_norm`: Scalar or vector of normalized power factors, typically in `[0, 1]`
* `temperature_ref`: Reference cell temperature (°C)
* `irradiance_ref`: Reference irradiance (W/m^2)
* `q`: Elementary charge (C)
* `kb`: Boltzmann constant (J/K)

# Returns:
* A vector of [`DCComponents`](@ref) in which `p_mp` has been multiplied by `power_norm`

```jldoctest
julia> using TimeZones

julia> power_norm = 0.9;

julia> pv_module = read_solar_module("Canadian Solar CS5P-220M [ 2009]");

julia> effective_irradiance = EffectiveIrradiance(ZonedDateTime(2023, 1, 1, tz"America/Denver"), 800);  # effective irradiance (W/m^2)

julia> cell_temperature = CellTemperature(ZonedDateTime(2023, 1, 1, tz"America/Denver"), 25);   # cell temperature (C)

julia> dc_unshaded = sapm_dc_components(
           pv_module,
           effective_irradiance,
           cell_temperature
       );

julia> dc_shaded = sapm_dc_components_shaded(
           pv_module,
           effective_irradiance,
           cell_temperature,
           power_norm
       );

julia> isapprox(dc_shaded.p_mp, dc_unshaded.p_mp * power_norm)
true
```
"""
function sapm_dc_components_shaded(
    pv_module::SolarModule,
    effective_irradiance::EffectiveIrradiance,
    cell_temperature::CellTemperature,
    power_norm::Real,
    temperature_ref::Real = 25.0,
    irradiance_ref::Real = 1000.0,
    q::Real = 1.602176634e-19,
    kb::Real = 1.380649e-23,
)

    dc_components = sapm_dc_components(
        pv_module,
        effective_irradiance,
        cell_temperature,
        temperature_ref,
        irradiance_ref,
        q,
        kb,
    )

    dc_components = DCComponents(
        dc_components.time,
        dc_components.i_sc,
        dc_components.i_mp,
        dc_components.v_oc,
        dc_components.v_mp,
        dc_components.p_mp * power_norm,
        dc_components.i_x,
        dc_components.i_xx,
    )

    return dc_components
end

function sapm_dc_components_shaded(
    pv_module::SolarModule,
    effective_irradiance::Vector{<:EffectiveIrradiance},
    cell_temperature::Vector{<:CellTemperature},
    power_norm::Vector{<:Real},
    temperature_ref::Real = 25.0,
    irradiance_ref::Real = 1000.0,
    q::Real = 1.602176634e-19,
    kb::Real = 1.380649e-23,
)

    n = length(effective_irradiance)

    length(cell_temperature) == n || throw(
        DimensionMismatch(
            "cell_temperature has $(length(cell_temperature)) rows but effective_irradiance has $n rows",
        ),
    )

    length(power_norm) == n || throw(
        DimensionMismatch(
            "power_norm has $(length(power_norm)) rows but effective_irradiance has $n rows",
        ),
    )

    dc_components = sapm_dc_components(
        pv_module,
        effective_irradiance,
        cell_temperature,
        temperature_ref,
        irradiance_ref,
        q,
        kb,
    )

    dc_components = [
        DCComponents(
            dc_components[i].time,
            dc_components[i].i_sc,
            dc_components[i].i_mp,
            dc_components[i].v_oc,
            dc_components[i].v_mp,
            dc_components[i].p_mp * power_norm[i],
            dc_components[i].i_x,
            dc_components[i].i_xx,
        ) for i in eachindex(effective_irradiance)
    ]

    return dc_components
end
