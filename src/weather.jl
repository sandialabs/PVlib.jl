# Set url for data access
const EUROPA_URL = "https://re.jrc.ec.europa.eu/api/"

# Set of inputs for NSRDB
const NSRDB_API_BASE = "https://developer.nlr.gov/api/nsrdb/v2/solar/"
const PSM4_TMY_ENDPOINT = "nsrdb-GOES-tmy-v4-0-0-download.csv"
const PSM4_TMY_URL = NSRDB_API_BASE * PSM4_TMY_ENDPOINT
const NSRDB_DEFAULT_PARAMETERS = (
    "air_temperature",
    "dew_point",
    "dhi",
    "dni",
    "ghi",
    "surface_albedo",
    "surface_pressure",
    "wind_direction",
    "wind_speed",
)
const NSRDB_REQUEST_VARIABLE_MAP = Dict(
    "ghi" => "ghi",
    "dhi" => "dhi",
    "dni" => "dni",
    "temp_air" => "air_temperature",
    "temp_dew" => "dew_point",
    "pressure" => "surface_pressure",
    "wind_speed" => "wind_speed",
    "wind_direction" => "wind_direction",
    "albedo" => "surface_albedo",
)

"""
    WeatherSample{T}

Meteorological inputs for a single timestamp.

Typically retrieved from the PVGIS TMY endpoint using `get_meteorological_data`.

# Fields:
* `time`: Timestamp as `ZonedDateTime`
* `ghi`: Global horizontal irradiance (W/m^2)
* `dni`: Direct normal irradiance (W/m^2)
* `dhi`: Diffuse horizontal irradiance (W/m^2)
* `temp_air`: Ambient air temperature (°C)
* `relative_humidity`: Relative humidity (%)
* `temp_dew`: Dew point temperature (°C)
* `pressure`: Surface pressure (Pa)
* `wind_speed`: Wind speed (m/s)
* `wind_direction`: Wind direction (degrees)
* `albedo`: Surface albedo (-)
"""
struct WeatherSample{T}
    time::ZonedDateTime
    ghi::Union{Missing,T}
    dni::Union{Missing,T}
    dhi::Union{Missing,T}
    temp_air::Union{Missing,T}
    temp_dew::Union{Missing,T}
    relative_humidity::Union{Missing,T}
    pressure::Union{Missing,T}
    wind_speed::Union{Missing,T}
    wind_direction::Union{Missing,T}
    albedo::Union{Missing,T}
end

# Constructor for WeatherSample that allows keyword arguments and converts to type T, with missing defaults.
function WeatherSample{T}(;
    time::ZonedDateTime,
    ghi = missing,
    dni = missing,
    dhi = missing,
    temp_air = missing,
    temp_dew = missing,
    relative_humidity = missing,
    pressure = missing,
    wind_speed = missing,
    wind_direction = missing,
    albedo = missing,
) where {T}

    return WeatherSample{T}(
        time,
        _toT(T, ghi),
        _toT(T, dni),
        _toT(T, dhi),
        _toT(T, temp_air),
        _toT(T, temp_dew),
        _toT(T, relative_humidity),
        _toT(T, pressure),
        _toT(T, wind_speed),
        _toT(T, wind_direction),
        _toT(T, albedo),
    )
end

# PVGIS sometimes returns numbers as Int/Float; treat anything Real as numeric.
_toT(::Type{T}, x) where {T} = x isa Real ? T(x) : missing

# inclusive date-range check
_in_range(t::ZonedDateTime, start_date::Date, end_date::Date) =
    (Date(t) >= start_date) && (Date(t) <= end_date)

# inclusive date-range check based on month and day
_in_monthday_range(t::ZonedDateTime, start_md::Tuple{Int,Int}, end_md::Tuple{Int,Int}) =
    begin
        md = (month(Date(t)), day(Date(t)))
        start_md <= md <= end_md
    end

# PVGIS "time(UTC)" looks like "yyyymmddd:HHMM" in your current code
function _parse_pvgis_time_utc(s::AbstractString)
    DateTime(s, dateformat"yyyymmddd:HHMM")
end

# Internal helper: call the PVGIS TMY endpoint and return parsed JSON.
function _fetch_tmy_data(
    latitude::Real,
    longitude::Real,
    require_ssl_verification::Bool = true,
)

    params = Dict("lat" => latitude, "lon" => longitude, "outputformat" => "json")
    response = get(
        EUROPA_URL * "tmy",
        query = params,
        require_ssl_verification = require_ssl_verification,
    )
    return parse(String(response.body))
end

# Internal helper: convert PVGIS JSON payload to a cleaned DataFrame with local_time and renamed columns.
function _process_tmy_data(
    src,
    start_date::Date,
    end_date::Date,
    timezone::TimeZone;
    T::Type = Float64,
)

    hourly = src["outputs"]["tmy_hourly"]  # array of dict-like objects

    out = WeatherSample{T}[]
    sizehint!(out, length(hourly))

    for rec in hourly
        # --- time ---
        t_utc = _parse_pvgis_time_utc(rec["time(UTC)"])
        t_loc = ZonedDateTime(t_utc, timezone; from_utc = true)

        _in_range(t_loc, start_date, end_date) || continue

        # --- read needed variables directly from PVGIS keys ---
        # We DO NOT need to "rename"; just pull from the PVGIS key names.
        ghi = _toT(T, Base.get(rec, "G(h)", missing))
        dni = _toT(T, Base.get(rec, "Gb(n)", missing))
        dhi = _toT(T, Base.get(rec, "Gd(h)", missing))
        Ta = _toT(T, Base.get(rec, "T2m", missing))
        rh = _toT(T, Base.get(rec, "RH", missing))
        sp = _toT(T, Base.get(rec, "SP", missing))
        ws = _toT(T, Base.get(rec, "WS10m", missing))
        wd = _toT(T, Base.get(rec, "WD10m", missing))

        push!(
            out,
            WeatherSample{T}(
                time = t_loc,
                ghi = ghi,
                dni = dni,
                dhi = dhi,
                temp_air = Ta,
                relative_humidity = rh,
                pressure = sp,
                wind_speed = ws,
                wind_direction = wd,
            ),
        )
    end

    return out
end

function _format_nsrdb_wkt(latitude::Real, longitude::Real)
    lon = round(longitude, digits = 4)
    lat = round(latitude, digits = 4)
    return "POINT($lon $lat)"
end

function _get_nsrdb_request_attributes(parameters)
    mapped =
        [Base.get(NSRDB_REQUEST_VARIABLE_MAP, String(p), String(p)) for p in parameters]
    return join(mapped, ",")
end

function _fetch_nsrdb_psm4_tmy(
    latitude::Real,
    longitude::Real,
    api_key::AbstractString,
    email::AbstractString;
    name::AbstractString = "tmy",
    time_step::Int = 60,
    parameters = NSRDB_DEFAULT_PARAMETERS,
    leap_day::Bool = false,
    full_name::AbstractString = "PVlib.jl",
    affiliation::AbstractString = "PVlib.jl",
    utc::Bool = false,
    require_ssl_verification::Bool = true,
    timeout::Real = 30,
    url::AbstractString = PSM4_TMY_URL,
)

    params = Dict(
        "api_key" => api_key,
        "full_name" => full_name,
        "email" => email,
        "affiliation" => affiliation,
        "reason" => "PVlib.jl",
        "mailing_list" => "false",
        "wkt" => _format_nsrdb_wkt(latitude, longitude),
        "names" => name,
        "attributes" => _get_nsrdb_request_attributes(parameters),
        "leap_day" => lowercase(string(leap_day)),
        "utc" => lowercase(string(utc)),
        "interval" => string(time_step),
    )

    response = get(
        url,
        query = params,
        require_ssl_verification = require_ssl_verification,
        timeout = timeout,
    )

    return String(response.body)
end

function _parse_nsrdb_value(
    row::AbstractDict{<:AbstractString,<:AbstractString},
    key::AbstractString,
)
    if !haskey(row, key)
        return missing
    end

    value = strip(String(row[key]))
    isempty(value) && return missing

    parsed = tryparse(Float64, value)
    return isnothing(parsed) ? missing : parsed
end

function _nsrdb_row_to_weather_sample(
    row::AbstractDict{<:AbstractString,<:AbstractString},
    timezone::TimeZone;
    T::Type = Float64,
)
    year = Base.parse(Int, strip(row["Year"]))
    month = Base.parse(Int, strip(row["Month"]))
    day = Base.parse(Int, strip(row["Day"]))
    hour = Base.parse(Int, strip(row["Hour"]))
    minute = Base.parse(Int, strip(row["Minute"]))

    t = ZonedDateTime(DateTime(year, month, day, hour, minute), timezone) - Minute(30)

    return WeatherSample{T}(
        time = t,
        ghi = _parse_nsrdb_value(row, "GHI"),
        dni = _parse_nsrdb_value(row, "DNI"),
        dhi = _parse_nsrdb_value(row, "DHI"),
        temp_air = _parse_nsrdb_value(row, "Temperature"),
        temp_dew = _parse_nsrdb_value(row, "Dew Point"),
        pressure = _parse_nsrdb_value(row, "Pressure"),
        wind_speed = _parse_nsrdb_value(row, "Wind Speed"),
        wind_direction = _parse_nsrdb_value(row, "Wind Direction"),
        albedo = _parse_nsrdb_value(row, "Surface Albedo"),
    )
end

function _parse_nsrdb_psm4_csv(
    csv_text::AbstractString;
    start_monthday::Union{Nothing,Tuple{Int,Int}} = nothing,
    end_monthday::Union{Nothing,Tuple{Int,Int}} = nothing,
    T::Type = Float64,
)

    lines = split(replace(csv_text, "\r\n" => "\n", '\r' => '\n'), '\n')
    lines = filter(line -> !isempty(strip(line)), lines)

    length(lines) < 4 && error("NSRDB CSV did not contain enough lines to parse.")

    columns = String.(filter(col -> !isempty(strip(col)), split(lines[3], ",")))

    fields = String.(split(lines[1], ","))
    values = String.(split(lines[2], ","))
    raw = Dict(strip(fields[i]) => strip(values[i]) for i in eachindex(fields))
    offset = Base.parse(Int, raw["Time Zone"])
    tz = TimeZone(offset == 0 ? "UTC" : "UTC$(offset)")

    out = WeatherSample{T}[]
    sizehint!(out, max(length(lines) - 3, 0))

    for line in lines[4:end]
        values = String.(split(line, ","))
        length(values) < length(columns) && continue

        row = Dict(columns[i] => values[i] for i in eachindex(columns))
        sample = _nsrdb_row_to_weather_sample(row, tz; T = T)

        if !isnothing(start_monthday) && !isnothing(end_monthday)
            _in_monthday_range(sample.time, start_monthday, end_monthday) || continue
        end

        push!(out, sample)
    end

    return out
end

"""
    get_meteorological_data_pvgis(
        latitude::Real,
        longitude::Real,
        start_date::Date,
        end_date::Date,
        timezone::TimeZone,
        require_ssl_verification::Bool=true
    ) -> Vector{<:WeatherSample}

Retrieve Typical Meteorological Year (TMY) data from the PVGIS API.

This function queries the PVGIS endpoint `https://re.jrc.ec.europa.eu/api/tmy`
for the specified location.

# Arguments:
* `latitude`: Site latitude (degrees)
* `longitude`: Site longitude (degrees)
* `start_date`: Start date (inclusive)
* `end_date`: End date (inclusive)
* `timezone`: Local timezone used to construct `time`
* `require_ssl_verification`: If `false`, disables TLS/SSL certificate verification for the HTTP request

# Returns:
* `weather_data`: Vector of `WeatherSample` containing the requested meteorological time series

# Notes:
* PVGIS timestamps are provided in UTC; this function converts them to `timezone`.

```markdown
```julia
julia> using Dates, TimeZones

julia> weather = get_meteorological_data_pvgis(35.1, -106.6, Date(2023, 1, 1), Date(2023, 1, 1), TimeZone("America/Denver"), false);

julia> fieldnames(eltype(weather)) == (:time,:ghi,:dni,:dhi,:temp_air,:temp_dew,:relative_humidity,:pressure,:wind_speed,:wind_direction,:albedo)
true

julia> isapprox(weather[11].ghi, 193.05; atol=1e-6)
true
```
"""
function get_meteorological_data_pvgis(
    latitude::Real,
    longitude::Real,
    start_date::Date,
    end_date::Date,
    timezone,
    require_ssl_verification::Bool = true,
)

    src = _fetch_tmy_data(latitude, longitude, require_ssl_verification)
    return _process_tmy_data(src, start_date, end_date, timezone; T = Float64)
end

"""
    get_meteorological_data_nsrdb(
        latitude::Real,
        longitude::Real,
        api_key::AbstractString,
        email::AbstractString,
        start_monthday::Union{Nothing,Tuple{Int,Int}}=(1, 1),
        end_monthday::Union{Nothing,Tuple{Int,Int}}=(12, 31),
        require_ssl_verification::Bool=true,
    ) -> Vector{<:WeatherSample}

Retrieve Typical Meteorological Year (TMY) data from the NSRDB PSM4 API.

This function queries the NSRDB GOES TMY endpoint for the specified location
and returns a vector of `WeatherSample`.

# Arguments:
* `latitude`: Site latitude (degrees)
* `longitude`: Site longitude (degrees)
* `api_key`: NSRDB API key from NLR
* `email`: Email address associated with the NSRDB request
* `start_monthday`: Start month/day as `(month, day)`; inclusive
* `end_monthday`: End month/day as `(month, day)`; inclusive
* `require_ssl_verification`: If `false`, disables TLS/SSL certificate verification for the HTTP request

# Returns:
* `weather_data`: Vector of `WeatherSample` containing the requested meteorological time series

# Notes:
* This function uses the NSRDB PSM4 TMY dataset, not a historical year-specific time series.
* Returned timestamps use the fixed local timezone reported by the NSRDB metadata.
* Filtering is performed by month/day, not by full calendar date.

```markdown
```julia
julia> weather_data = get_meteorological_data_nsrdb(
            35.1,
            -106.6,
            "E52b7mqeTWLigj2xF5Bn4n6Mm87ecm5LFFeYh4US",
            "jtgrasb@sandia.gov",
            (1, 1),
            (1, 1),
            false,
        );
        
julia> isapprox(weather_data[11].ghi, 490.0; atol=1e-6)
true
```
"""
function get_meteorological_data_nsrdb(
    latitude::Real,
    longitude::Real,
    api_key::AbstractString,
    email::AbstractString,
    start_monthday::Union{Nothing,Tuple{Int,Int}} = (1, 1),
    end_monthday::Union{Nothing,Tuple{Int,Int}} = (12, 31),
    require_ssl_verification::Bool = true,
)

    csv_text = _fetch_nsrdb_psm4_tmy(
        latitude,
        longitude,
        api_key,
        email;
        require_ssl_verification = require_ssl_verification,
    )

    return _parse_nsrdb_psm4_csv(
        csv_text;
        start_monthday = start_monthday,
        end_monthday = end_monthday,
        T = Float64,
    )
end

# Internal helper: compute relative airmass from zenith angle (scalar).
function _get_relative_airmass(zenith::Real)
    z = isnan(zenith) ? NaN : zenith
    z = z > 90 ? NaN : z  # Replace values greater than 90 with NaN
    relative_airmass = 1.0 / (cosd(z) + 0.50572 * ((6.07995 + (90 - z))^-1.6364))
    return relative_airmass
end

# Internal helper: compute relative airmass from zenith angle(s) (vector).
function _get_relative_airmass(zenith::AbstractVector{<:Real})
    relative_airmass = [_get_relative_airmass(z) for z in zenith]
    return relative_airmass
end

# Internal helper: estimate pressure (Pa) from altitude (m).
function _get_pressure(altitude::Real)
    pressure = 100 * ((44331.514 - altitude) / 11880.516)^(1 / 0.1902632)
    return pressure
end

# Internal helper: compute absolute airmass from zenith (deg) and altitude (m) (scalar).
function _get_absolute_airmass(zenith::Real, altitude::Real)
    relative_airmass = _get_relative_airmass(zenith)
    pressure = _get_pressure(altitude)
    absolute_airmass = relative_airmass * (pressure / 101325.0)
    return absolute_airmass
end

# Internal helper: compute absolute airmass from zenith vector (deg) and altitude (m).
function _get_absolute_airmass(zenith::AbstractVector{<:Real}, altitude::Real)
    absolute_airmass = [_get_absolute_airmass(z, altitude) for z in zenith]
    return absolute_airmass
end

# WeatherSample pretty printer (shared implementation)
function _show_weather_pretty(io::IO, w::WeatherSample)
    println(io, "WeatherSample")
    println(io, "────────────")
    println(io, "time               │ ", w.time)
    println(io, "ghi                │ ", w.ghi)
    println(io, "dni                │ ", w.dni)
    println(io, "dhi                │ ", w.dhi)
    println(io, "temp_air           │ ", w.temp_air)
    println(io, "temp_dew           │ ", w.temp_dew)
    println(io, "relative_humidity  │ ", w.relative_humidity)
    println(io, "pressure           │ ", w.pressure)
    println(io, "wind_speed         │ ", w.wind_speed)
    println(io, "wind_direction     │ ", w.wind_direction)
    println(io, "albedo             │ ", w.albedo)
end

function Base.show(io::IO, w::WeatherSample)
    _show_weather_pretty(io, w)
end

function Base.show(io::IO, ::MIME"text/plain", w::WeatherSample)
    _show_weather_pretty(io, w)
end

# WeatherSample vector pretty printer (shared)
function _show_weather_vec_pretty(io::IO, v::AbstractVector{<:WeatherSample})
    header = [
        "time",
        "ghi",
        "dni",
        "dhi",
        "temp_air",
        "temp_dew",
        "relative_humidity",
        "pressure",
        "wind_speed",
        "wind_direction",
        "albedo",
    ]

    _show_table(
        io,
        "WeatherSample",
        header,
        i -> begin
            x = v[i]
            [
                string(x.time),
                string(x.ghi),
                string(x.dni),
                string(x.dhi),
                string(x.temp_air),
                string(x.temp_dew),
                string(x.relative_humidity),
                string(x.pressure),
                string(x.wind_speed),
                string(x.wind_direction),
                string(x.albedo),
            ]
        end,
        length(v),
    )
end

function Base.show(io::IO, v::AbstractVector{<:WeatherSample})
    _show_weather_vec_pretty(io, v)
end

function Base.show(io::IO, ::MIME"text/plain", v::AbstractVector{<:WeatherSample})
    _show_weather_vec_pretty(io, v)
end
