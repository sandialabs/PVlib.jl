# Include the SPA coefficients
include("spa_coefficients.jl")

"""
    SolarPosition{T}

Solar position quantities at a single timestamp.

Typically computed from location and weather inputs using `get_solar_position`.

# Fields:
* `time`: Timestamp associated with this solution
* `apparent_zenith`: Apparent zenith angle (degrees)
* `zenith`: Geometric zenith angle (degrees)
* `apparent_elevation`: Apparent elevation angle (degrees)
* `elevation`: Geometric elevation angle (degrees)
* `azimuth`: Azimuth angle (degrees)
* `equation_of_time`: Equation of time (minutes)
"""
struct SolarPosition{T}
    time::Any
    apparent_zenith::T
    zenith::T
    apparent_elevation::T
    elevation::T
    azimuth::T
    equation_of_time::T

    function SolarPosition(
        time,
        apparent_zenith,
        zenith,
        apparent_elevation,
        elevation,
        azimuth,
        equation_of_time,
    )
        # Promote all numeric outputs to a common type (helps type-stability)
        T = promote_type(
            typeof(apparent_zenith),
            typeof(zenith),
            typeof(apparent_elevation),
            typeof(elevation),
            typeof(azimuth),
            typeof(equation_of_time),
        )

        apz, z, ape, e, az, eot = convert.(
            T,
            (
                apparent_zenith,
                zenith,
                apparent_elevation,
                elevation,
                azimuth,
                equation_of_time,
            ),
        )

        # Example sanity checks (edit/remove as appropriate)
        @assert isfinite(apz) &&
                isfinite(z) &&
                isfinite(ape) &&
                isfinite(e) &&
                isfinite(az) &&
                isfinite(eot) "Non-finite solar position value."
        # optional: enforce typical angle ranges if you want:
        # @assert 0 ≤ z ≤ 180 "zenith out of range"

        return new{T}(time, apz, z, ape, e, az, eot)
    end
end

# Function to calculate the Julian Day from Unix time
function _julian_day(unixtime::Real)
    return unixtime / 86400 + 2440587.5
end

# Function to calculate the Julian Ephemeris Day
function _julian_ephemeris_day(julian_day::Real, delta_t::Real)
    jde = julian_day + delta_t / 86400.0
    return jde
end

# Function to calculate the Julian Century
function _julian_century(julian_day::Real)
    jc = (julian_day - 2451545.0) / 36525.0
    return jc
end

# Function to calculate the Julian Ephemeris Century
function _julian_ephemeris_century(julian_ephemeris_day::Real)
    jce = (julian_ephemeris_day - 2451545.0) / 36525.0
    return jce
end

# Function to calculate the Julian Ephemeris Millennium
function _julian_ephemeris_millennium(julian_ephemeris_century::Real)
    jme = julian_ephemeris_century / 10.0
    return jme
end

# Helper function to sum and multiply cosine values
function _sum_mult_cos_add_mult(arr::AbstractMatrix{<:Real}, x::Real)
    s = 0.0
    for row in eachrow(arr)
        s += row[1] * cos(row[2] + row[3] * x)
    end
    return s
end

# Function to calculate the heliocentric longitude
function _heliocentric_longitude(jme::Real)
    l0 = _sum_mult_cos_add_mult(L0, jme)
    l1 = _sum_mult_cos_add_mult(L1, jme)
    l2 = _sum_mult_cos_add_mult(L2, jme)
    l3 = _sum_mult_cos_add_mult(L3, jme)
    l4 = _sum_mult_cos_add_mult(L4, jme)
    l5 = _sum_mult_cos_add_mult(L5, jme)

    l_rad = (l0 + l1 * jme + l2 * jme^2 + l3 * jme^3 + l4 * jme^4 + l5 * jme^5) / 10^8
    l = rad2deg(l_rad)
    return l % 360
end

function _heliocentric_latitude(jme::Real)
    b0 = _sum_mult_cos_add_mult(B0, jme)
    b1 = _sum_mult_cos_add_mult(B1, jme)

    b_rad = (b0 + b1 * jme) / 10^8
    b = rad2deg(b_rad)
    return b
end

# Function to calculate the heliocentric radius vector
function _heliocentric_radius_vector(jme::Real)
    r0 = _sum_mult_cos_add_mult(R0, jme)
    r1 = _sum_mult_cos_add_mult(R1, jme)
    r2 = _sum_mult_cos_add_mult(R2, jme)
    r3 = _sum_mult_cos_add_mult(R3, jme)
    r4 = _sum_mult_cos_add_mult(R4, jme)

    r = (r0 + r1 * jme + r2 * jme^2 + r3 * jme^3 + r4 * jme^4) / 10^8
    return r
end

# Function to calculate the Geocentric Longitude
function _geocentric_longitude(heliocentric_longitude::Real)
    theta = heliocentric_longitude + 180.0
    return theta % 360.0
end

# Function to calculate the Geocentric Latitude
function _geocentric_latitude(heliocentric_latitude::Real)
    beta = -heliocentric_latitude
    return beta
end

# Function to calculate the Mean Elongation
function _mean_elongation(julian_ephemeris_century::Real)
    x0 = (
        297.85036 + 445267.111480 * julian_ephemeris_century -
        0.0019142 * julian_ephemeris_century^2 + julian_ephemeris_century^3 / 189474
    )
    return x0
end

# Function to calculate the Mean Anomaly of the Sun
function _mean_anomaly_sun(julian_ephemeris_century::Real)
    x1 = (
        357.52772 + 35999.050340 * julian_ephemeris_century -
        0.0001603 * julian_ephemeris_century^2 - julian_ephemeris_century^3 / 300000
    )
    return x1
end

# Function to calculate the Mean Anomaly of the Moon
function _mean_anomaly_moon(julian_ephemeris_century::Real)
    x2 = (
        134.96298 +
        477198.867398 * julian_ephemeris_century +
        0.0086972 * julian_ephemeris_century^2 +
        julian_ephemeris_century^3 / 56250
    )
    return x2
end

# Function to calculate the Moon Argument Latitude
function _moon_argument_latitude(julian_ephemeris_century::Real)
    x3 = (
        93.27191 + 483202.017538 * julian_ephemeris_century -
        0.0036825 * julian_ephemeris_century^2 + julian_ephemeris_century^3 / 327270
    )
    return x3
end

# Function to calculate the Moon Ascending Longitude
function _moon_ascending_longitude(julian_ephemeris_century::Real)
    x4 = (
        125.04452 - 1934.136261 * julian_ephemeris_century +
        0.0020708 * julian_ephemeris_century^2 +
        julian_ephemeris_century^3 / 450000
    )
    return x4
end

# Function to calculate longitude obliquity nutation
function _longitude_obliquity_nutation(
    julian_ephemeris_century::Real,
    x0::Real,
    x1::Real,
    x2::Real,
    x3::Real,
    x4::Real,
)
    delta_psi_sum = 0.0
    delta_eps_sum = 0.0

    for row = 1:size(NUTATION_YTERM_ARRAY, 1)
        a = NUTATION_ABCD_ARRAY[row, 1]
        b = NUTATION_ABCD_ARRAY[row, 2]
        c = NUTATION_ABCD_ARRAY[row, 3]
        d = NUTATION_ABCD_ARRAY[row, 4]

        arg = deg2rad(
            NUTATION_YTERM_ARRAY[row, 1] * x0 +
            NUTATION_YTERM_ARRAY[row, 2] * x1 +
            NUTATION_YTERM_ARRAY[row, 3] * x2 +
            NUTATION_YTERM_ARRAY[row, 4] * x3 +
            NUTATION_YTERM_ARRAY[row, 5] * x4,
        )

        delta_psi_sum += (a + b * julian_ephemeris_century) * sin(arg)
        delta_eps_sum += (c + d * julian_ephemeris_century) * cos(arg)
    end

    delta_psi = delta_psi_sum / 36000000.0
    delta_eps = delta_eps_sum / 36000000.0

    return delta_psi, delta_eps
end

# Function to calculate the Mean Ecliptic Obliquity
function _mean_ecliptic_obliquity(julian_ephemeris_millennium::Real)
    U = julian_ephemeris_millennium / 10.0
    e0 = (
        84381.448 - 4680.93 * U - 1.55 * U^2 + 1999.25 * U^3 - 51.38 * U^4 - 249.67 * U^5 - 39.05 * U^6 +
        7.12 * U^7 +
        27.87 * U^8 +
        5.79 * U^9 +
        2.45 * U^10
    )
    return e0
end

# Function to calculate the True Ecliptic Obliquity
function _true_ecliptic_obliquity(mean_ecliptic_obliquity::Real, obliquity_nutation::Real)
    e0 = mean_ecliptic_obliquity
    deleps = obliquity_nutation
    e = e0 / 3600.0 + deleps
    return e
end

# Function to calculate the Aberration Correction
function _aberration_correction(earth_radius_vector::Real)
    deltau = -20.4898 / (3600 * earth_radius_vector)
    return deltau
end

# Function to calculate the Apparent Sun Longitude
function _apparent_sun_longitude(
    geocentric_longitude::Real,
    longitude_nutation::Real,
    aberration_correction::Real,
)
    lamd = geocentric_longitude + longitude_nutation + aberration_correction
    return lamd
end

# Function to calculate the Mean Sidereal Time
function _mean_sidereal_time(julian_day::Real, julian_century::Real)
    v0 = (
        280.46061837 +
        360.98564736629 * (julian_day - 2451545) +
        0.000387933 * julian_century^2 - julian_century^3 / 38710000
    )
    return v0 % 360.0
end

# Function to calculate the Apparent Sidereal Time
function _apparent_sidereal_time(
    mean_sidereal_time::Real,
    longitude_nutation::Real,
    true_ecliptic_obliquity::Real,
)
    v = mean_sidereal_time + longitude_nutation * cos(deg2rad(true_ecliptic_obliquity))
    return v
end

# Function to calculate the Geocentric Sun Right Ascension
function _geocentric_sun_right_ascension(
    apparent_sun_longitude::Real,
    true_ecliptic_obliquity::Real,
    geocentric_latitude::Real,
)
    true_ecliptic_obliquity_rad = deg2rad(true_ecliptic_obliquity)
    apparent_sun_longitude_rad = deg2rad(apparent_sun_longitude)

    num = (
        sin(apparent_sun_longitude_rad) * cos(true_ecliptic_obliquity_rad) -
        tan(deg2rad(geocentric_latitude)) * sin(true_ecliptic_obliquity_rad)
    )
    alpha = rad2deg(atan(num, cos(apparent_sun_longitude_rad)))
    return alpha % 360
end

# Function to calculate the Geocentric Sun Declination
function _geocentric_sun_declination(
    apparent_sun_longitude::Real,
    true_ecliptic_obliquity::Real,
    geocentric_latitude::Real,
)
    geocentric_latitude_rad = deg2rad(geocentric_latitude)
    true_ecliptic_obliquity_rad = deg2rad(true_ecliptic_obliquity)

    delta = rad2deg(
        asin(
            sin(geocentric_latitude_rad) * cos(true_ecliptic_obliquity_rad) +
            cos(geocentric_latitude_rad) *
            sin(true_ecliptic_obliquity_rad) *
            sin(deg2rad(apparent_sun_longitude)),
        ),
    )
    return delta
end

# Function to calculate the Local Hour Angle
function _local_hour_angle(
    apparent_sidereal_time::Real,
    observer_longitude::Real,
    sun_right_ascension::Real,
)
    """Measured westward from south"""
    H = apparent_sidereal_time + observer_longitude - sun_right_ascension
    return H % 360.0
end

# Function to calculate the Equatorial Horizontal Parallax
function _equatorial_horizontal_parallax(earth_radius_vector::Real)
    xi = 8.794 / (3600 * earth_radius_vector)
    return xi
end

# Function to calculate the U term
function _uterm(observer_latitude::Real)
    u = atan(0.99664719 * tan(deg2rad(observer_latitude)))
    return u
end

# Function to calculate the X term
function _xterm(u::Real, observer_latitude::Real, observer_elevation::Real)
    x = cos(u) + observer_elevation / 6378140 * cos(deg2rad(observer_latitude))
    return x
end

# Function to calculate the Y term
function _yterm(u::Real, observer_latitude::Real, observer_elevation::Real)
    y = 0.99664719 * sin(u) + observer_elevation / 6378140 * sin(deg2rad(observer_latitude))
    return y
end

# Function to calculate the Parallax in Sun Right Ascension
function _parallax_sun_right_ascension(
    xterm::Real,
    equatorial_horizontal_parallax::Real,
    local_hour_angle::Real,
    geocentric_sun_declination::Real,
)
    equatorial_horizontal_parallax_rad = deg2rad(equatorial_horizontal_parallax)
    local_hour_angle_rad = deg2rad(local_hour_angle)

    num = -xterm * sin(equatorial_horizontal_parallax_rad) * sin(local_hour_angle_rad)
    denom =
        cos(deg2rad(geocentric_sun_declination)) -
        xterm * sin(equatorial_horizontal_parallax_rad) * cos(local_hour_angle_rad)
    delta_alpha = rad2deg(atan(num, denom))
    return delta_alpha
end

# Function to calculate the Topocentric Sun Right Ascension
function _topocentric_sun_right_ascension(
    geocentric_sun_right_ascension::Real,
    parallax_sun_right_ascension::Real,
)
    alpha_prime = geocentric_sun_right_ascension + parallax_sun_right_ascension
    return alpha_prime
end

# Function to calculate the Topocentric Sun Declination
function _topocentric_sun_declination(
    geocentric_sun_declination::Real,
    xterm::Real,
    yterm::Real,
    equatorial_horizontal_parallax::Real,
    parallax_sun_right_ascension::Real,
    local_hour_angle::Real,
)
    geocentric_sun_declination_rad = deg2rad(geocentric_sun_declination)
    equatorial_horizontal_parallax_rad = deg2rad(equatorial_horizontal_parallax)

    num = (
        (
            sin(geocentric_sun_declination_rad) -
            yterm * sin(equatorial_horizontal_parallax_rad)
        ) * cos(deg2rad(parallax_sun_right_ascension))
    )
    denom = (
        cos(geocentric_sun_declination_rad) -
        xterm * sin(equatorial_horizontal_parallax_rad) * cos(deg2rad(local_hour_angle))
    )
    delta = rad2deg(atan(num, denom))
    return delta
end

# Function to calculate the Topocentric Local Hour Angle
function _topocentric_local_hour_angle(
    local_hour_angle::Real,
    parallax_sun_right_ascension::Real,
)
    H_prime = local_hour_angle - parallax_sun_right_ascension
    return H_prime
end

# Function to calculate the Topocentric Elevation Angle without Atmosphere
function _topocentric_elevation_angle_without_atmosphere(
    observer_latitude::Real,
    topocentric_sun_declination::Real,
    topocentric_local_hour_angle::Real,
)
    observer_latitude_rad = deg2rad(observer_latitude)
    topocentric_sun_declination_rad = deg2rad(topocentric_sun_declination)

    e0 = rad2deg(
        asin(
            sin(observer_latitude_rad) * sin(topocentric_sun_declination_rad) +
            cos(observer_latitude_rad) *
            cos(topocentric_sun_declination_rad) *
            cos(deg2rad(topocentric_local_hour_angle)),
        ),
    )
    return e0
end

# Function to calculate the Atmospheric Refraction Correction
function _atmospheric_refraction_correction(
    local_pressure::Real,
    local_temp::Real,
    topocentric_elevation_angle_wo_atmosphere::Real,
    atmos_refract::Real,
)
    # Switch sets delta_e when the sun is below the horizon
    switch = topocentric_elevation_angle_wo_atmosphere >= -1.0 * (0.26667 + atmos_refract)
    delta_e = (
        (local_pressure / 1010.0) * (283.0 / (273 + local_temp)) * 1.02 / (
            60 * tan(
                deg2rad(
                    topocentric_elevation_angle_wo_atmosphere +
                    10.3 / (topocentric_elevation_angle_wo_atmosphere + 5.11),
                ),
            )
        ) * switch
    )
    return delta_e
end

# Function to calculate the Topocentric Elevation Angle
function _topocentric_elevation_angle(
    topocentric_elevation_angle_without_atmosphere::Real,
    atmospheric_refraction_correction::Real,
)
    e = topocentric_elevation_angle_without_atmosphere + atmospheric_refraction_correction
    return e
end

# Function to calculate the Topocentric Zenith Angle
function _topocentric_zenith_angle(topocentric_elevation_angle::Real)
    theta = 90 - topocentric_elevation_angle
    return theta
end

# Function to calculate the Topocentric Astronomers Azimuth
function _topocentric_astronomers_azimuth(
    topocentric_local_hour_angle::Real,
    topocentric_sun_declination::Real,
    observer_latitude::Real,
)
    topocentric_local_hour_angle_rad = deg2rad(topocentric_local_hour_angle)
    observer_latitude_rad = deg2rad(observer_latitude)

    num = sin(topocentric_local_hour_angle_rad)
    denom = (
        cos(topocentric_local_hour_angle_rad) * sin(observer_latitude_rad) -
        tan(deg2rad(topocentric_sun_declination)) * cos(observer_latitude_rad)
    )
    gamma = rad2deg(atan(num, denom))
    return gamma % 360
end

# Function to calculate the Topocentric Azimuth Angle
function _topocentric_azimuth_angle(topocentric_astronomers_azimuth::Real)
    phi = topocentric_astronomers_azimuth + 180
    return phi % 360
end

# Function to calculate the Sun Mean Longitude
function _sun_mean_longitude(julian_ephemeris_millennium::Real)
    M = (
        280.4664567 +
        360007.6982779 * julian_ephemeris_millennium +
        0.03032028 * julian_ephemeris_millennium^2 +
        julian_ephemeris_millennium^3 / 49931 - julian_ephemeris_millennium^4 / 15300 -
        julian_ephemeris_millennium^5 / 2000000
    )
    return M
end

# Function to calculate the Equation of Time
function _equation_of_time(
    sun_mean_longitude::Real,
    geocentric_sun_right_ascension::Real,
    longitude_nutation::Real,
    true_ecliptic_obliquity::Real,
)
    E = (
        sun_mean_longitude - 0.0057183 - geocentric_sun_right_ascension +
        longitude_nutation * cos(deg2rad(true_ecliptic_obliquity))
    )
    # Limit between 0 and 360
    E = E % 360.0
    # Convert to minutes
    E *= 4.0
    greater = E > 20.0
    less = E < -20.0
    other = (E <= 20.0) && (E >= -20.0)
    E = greater * (E - 1440.0) + less * (E + 1440.0) + other * E
    return E
end

# Function to calculate the solar position
function _get_solar_position(
    unixtime::Real,
    lat::Real,
    lon::Real,
    elev::Real,
    pressure::Real,
    temp::Real,
    delta_t::Real,
    atmos_refract::Real,
)
    jd = _julian_day(unixtime)
    jde = _julian_ephemeris_day(jd, delta_t)
    jc = _julian_century(jd)
    jce = _julian_ephemeris_century(jde)
    jme = _julian_ephemeris_millennium(jce)
    R = _heliocentric_radius_vector(jme)

    L = _heliocentric_longitude(jme)
    B = _heliocentric_latitude(jme)
    Theta = _geocentric_longitude(L)
    beta = _geocentric_latitude(B)
    x0 = _mean_elongation(jce)
    x1 = _mean_anomaly_sun(jce)
    x2 = _mean_anomaly_moon(jce)
    x3 = _moon_argument_latitude(jce)
    x4 = _moon_ascending_longitude(jce)

    delta_psi, delta_epsilon = _longitude_obliquity_nutation(jce, x0, x1, x2, x3, x4)

    epsilon0 = _mean_ecliptic_obliquity(jme)
    epsilon = _true_ecliptic_obliquity(epsilon0, delta_epsilon)
    delta_tau = _aberration_correction(R)
    lamd = _apparent_sun_longitude(Theta, delta_psi, delta_tau)
    v0 = _mean_sidereal_time(jd, jc)
    v = _apparent_sidereal_time(v0, delta_psi, epsilon)
    alpha = _geocentric_sun_right_ascension(lamd, epsilon, beta)
    delta = _geocentric_sun_declination(lamd, epsilon, beta)

    m = _sun_mean_longitude(jme)
    eot = _equation_of_time(m, alpha, delta_psi, epsilon)
    H = _local_hour_angle(v, lon, alpha)
    xi = _equatorial_horizontal_parallax(R)
    u = _uterm(lat)
    x = _xterm(u, lat, elev)
    y = _yterm(u, lat, elev)
    delta_alpha = _parallax_sun_right_ascension(x, xi, H, delta)
    delta_prime = _topocentric_sun_declination(delta, x, y, xi, delta_alpha, H)
    H_prime = _topocentric_local_hour_angle(H, delta_alpha)
    e0 = _topocentric_elevation_angle_without_atmosphere(lat, delta_prime, H_prime)
    delta_e = _atmospheric_refraction_correction(pressure, temp, e0, atmos_refract)
    e = _topocentric_elevation_angle(e0, delta_e)
    theta = _topocentric_zenith_angle(e)
    theta0 = _topocentric_zenith_angle(e0)
    gamma = _topocentric_astronomers_azimuth(H_prime, delta_prime, lat)
    phi = _topocentric_azimuth_angle(gamma)

    return theta, theta0, e, e0, phi, eot
end

"""
    get_solar_position(
        lat::Real,
        lon::Real,
        elev::Real,
        weather_data::WeatherSample,
        delta_t::Real=67.0,
        atmos_refract::Real=0.5667,
        epoch::DateTime=DateTime(1970, 1, 1, 0, 0)
    ) -> SolarPosition

    get_solar_position(
        lat::Real,
        lon::Real,
        elev::Real,
        weather_data::AbstractVector{<:WeatherSample},
        delta_t::Real=67.0,
        atmos_refract::Real=0.5667,
        epoch::DateTime=DateTime(1970, 1, 1, 0, 0)
    ) -> Vector{<:SolarPosition}

Compute solar position quantities using the National Renewable Energy Laboratory
Solar Position Algorithm (SPA) and meteorological inputs.

The input `weather_data` must provide `time` (as a `ZonedDateTime`) and the
meteorological quantities needed for refraction (`pressure`, `temp_air`).
The time is converted to UTC and then to Unix seconds relative to `epoch`.

# Arguments:
* `lat`: Observer latitude (degrees)
* `lon`: Observer longitude (degrees)
* `elev`: Observer elevation (meters)
* `weather_data`: Either a `WeatherSample` (single timestamp) or vector of `WeatherSample` (time series)
* `delta_t`: (ΔT) (seconds), difference between Terrestrial Time and UT1
* `atmos_refract`: Atmospheric refraction parameter (degrees); used in the refraction switch/threshold
* `epoch`: Reference epoch for converting `DateTime` to Unix time (default Unix epoch)

# Returns:
* `SolarPosition` for a single `WeatherSample`, or `Vector{<:SolarPosition}` for a time series.
  Angles are in degrees and equation of time is in minutes.

```jldoctest
julia> using TimeZones, Dates

julia> w = WeatherSample{Float64}(
           ZonedDateTime(2020, 6, 1, 12, 0, 0, tz"UTC"),
           1.0, 1.0, 1.0,               # ghi, dni, dhi (unused here)
           20.0,                        # temp_air (C)
           1.0,                         # relative_humidity (unused here)
           101325.0,                    # pressure (Pa)
           1.0, 1.0                     # wind_speed, wind_direction (unused here)
       );

julia> sp = get_solar_position(35.1, -106.6, 1600.0, w);

julia> (0.0 <= sp.azimuth < 360.0) && (0.0 <= sp.apparent_zenith <= 180.0)
true

julia> isapprox(sp.apparent_elevation, 0.8095246595531291; atol=1e-6)
true
```
"""
# Function to calculate the solar position for float inputs
function get_solar_position(
    lat::Real,
    lon::Real,
    elev::Real,
    weather_data::WeatherSample,
    delta_t::Real = 67.0,
    atmos_refract::Real = 0.5667,
    epoch::DateTime = DateTime(1970, 1, 1, 0, 0),
)

    utc_time = DateTime(weather_data.time, UTC)
    unixtime = (utc_time - epoch) / Second(1)

    theta, theta0, e, e0, phi, eot = _get_solar_position(
        unixtime,
        lat,
        lon,
        elev,
        weather_data.pressure / 100,
        weather_data.temp_air,
        delta_t,
        atmos_refract,
    )

    return SolarPosition(utc_time, theta, theta0, e, e0, phi, eot)
end

# Function to calculate the solar position for vector inputs
function get_solar_position(
    lat::Real,
    lon::Real,
    elev::Real,
    weather_data::AbstractVector{<:WeatherSample},
    delta_t::Real = 67.0,
    atmos_refract::Real = 0.5667,
    epoch::DateTime = DateTime(1970, 1, 1, 0, 0),
)

    utc_time = [DateTime(w.time, UTC) for w in weather_data]
    unixtime = [(t - epoch) / Second(1) for t in utc_time]

    # Infer a concrete element type for the vector (helps performance)
    T = promote_type(
        typeof(lat),
        typeof(lon),
        typeof(elev),
        typeof(delta_t),
        typeof(atmos_refract),
        Float64,
    )  # unixtime is typically Float64

    solar_position = Vector{SolarPosition{T}}(undef, length(unixtime))

    for ii in eachindex(unixtime)
        θ, θ0, e, e0, ϕ, eot = _get_solar_position(
            unixtime[ii],
            lat,
            lon,
            elev,
            weather_data[ii].pressure / 100,
            weather_data[ii].temp_air,
            delta_t,
            atmos_refract,
        )
        solar_position[ii] = SolarPosition(utc_time[ii], θ, θ0, e, e0, ϕ, eot)
    end

    return solar_position
end

# SolarPosition pretty printer (shared implementation)
function _show_solpos_pretty(io::IO, sp::SolarPosition)
    println(io, "SolarPosition at ", sp.time)
    println(io, "────────────────────────────")
    println(io, "apparent_zenith    │ ", sp.apparent_zenith)
    println(io, "zenith             │ ", sp.zenith)
    println(io, "apparent_elevation │ ", sp.apparent_elevation)
    println(io, "elevation          │ ", sp.elevation)
    println(io, "azimuth            │ ", sp.azimuth)
    println(io, "equation_of_time   │ ", sp.equation_of_time)
end

function Base.show(io::IO, sp::SolarPosition)
    _show_solpos_pretty(io, sp)
end

function Base.show(io::IO, ::MIME"text/plain", sp::SolarPosition)
    _show_solpos_pretty(io, sp)
end

# SolarPosition vector pretty printer (shared implementation)
function _show_solpos_vec_pretty(io::IO, v::Vector{<:SolarPosition})
    header = [
        "time",
        "apparent_zenith",
        "zenith",
        "apparent_elevation",
        "elevation",
        "azimuth",
        "equation_of_time",
    ]
    _show_table(
        io,
        "SolarPosition",
        header,
        i -> begin
            x = v[i]
            [
                string(x.time),
                string(x.apparent_zenith),
                string(x.zenith),
                string(x.apparent_elevation),
                string(x.elevation),
                string(x.azimuth),
                string(x.equation_of_time),
            ]
        end,
        length(v),
    )
end

function Base.show(io::IO, v::Vector{<:SolarPosition})
    _show_solpos_vec_pretty(io, v)
end

function Base.show(io::IO, ::MIME"text/plain", v::Vector{<:SolarPosition})
    _show_solpos_vec_pretty(io, v)
end
