"""
    SolarModule{T}

Container for SAM/Sandia PV module parameters.

These parameters are typically loaded from the SAM/Sandia module library CSV via
`read_solar_module`.

# Fields:
* `name`: Module name
* `vintage`: Module vintage/year string
* `area`: Module area (m^2)
* `material`: Cell material description
* `cells_in_series`: Number of cells in series
* `parallel_strings`: Number of parallel strings
* `isco`, `voco`, `impo`, `vmpo`: Reference current/voltage points (A, V)
* `aisc`, `aimp`: Temperature coefficients for current
* `c0`, `c1`, `c2`, `c3`: SAPM empirical coefficients
* `bvoco`, `mbvoc`, `bvmpo`, `mbvmp`: Temperature/irradiance voltage coefficients
* `n`: Diode factor used in SAPM voltage terms
* `a0`–`a4`: Spectral-loss polynomial coefficients
* `b0`–`b5`: Incidence-angle modifier (IAM) polynomial coefficients
* `dtc`: Cell temperature difference parameter (SAPM)
* `fd`: Fraction of diffuse irradiance used by SAPM (unitless)
* `a`, `b`: SAPM temperature model coefficients (when applicable)
* `c4`, `c5`, `c6`, `c7`: Coefficients for optional `i_x` and `i_xx` terms
* `ixo`, `ixxo`: Reference values for optional `i_x` and `i_xx` terms
* `notes`: Freeform notes from the library
"""
struct SolarModule{T}
    name::String
    vintage::String
    area::T
    material::String
    cells_in_series::Int
    parallel_strings::Int
    isco::T
    voco::T
    impo::T
    vmpo::T
    aisc::T
    aimp::T
    c0::T
    c1::T
    bvoco::T
    mbvoc::T
    bvmpo::T
    mbvmp::T
    n::T
    c2::T
    c3::T
    a0::T
    a1::T
    a2::T
    a3::T
    a4::T
    b0::T
    b1::T
    b2::T
    b3::T
    b4::T
    b5::T
    dtc::T
    fd::T
    a::T
    b::T
    c4::T
    c5::T
    ixo::T
    ixxo::T
    c6::T
    c7::T
    notes::String
end

"""
    EffectiveIrradiance{T}

SAPM effective irradiance for a single timestamp.

Typically computed from POA irradiance components, module parameters, and solar
position using `sapm_effective_irradiance`.

# Fields:
* `time`: Timestamp
* `effective_irradiance`: Effective irradiance (SAPM scaling), typically in W/m^2
"""
struct EffectiveIrradiance{T}
    time::Any
    effective_irradiance::T
end

function rrule(::Type{EffectiveIrradiance}, time, effective_irradiance)
    y = EffectiveIrradiance{typeof(effective_irradiance)}(time, effective_irradiance)

    function pullback(ȳ)
        return NoTangent(), NoTangent(), getproperty(ȳ, :effective_irradiance)
    end

    return y, pullback
end

"""
    DCComponents{T}

SAPM DC operating point summary for a single timestamp.

Typically computed using `sapm_dc_components`.

# Fields:
* `time`: Timestamp
* `i_sc`: Short-circuit current (A)
* `i_mp`: Current at maximum power (A)
* `v_oc`: Open-circuit voltage (V)
* `v_mp`: Voltage at maximum power (V)
* `p_mp`: Power at maximum power (W)
* `i_x`: Optional current at (V = 0.5 V_{oc}) (A) (may be `nothing`)
* `i_xx`: Optional current at (V = 0.5(V_{oc} + V_{mp})) (A) (may be `nothing`)
"""
struct DCComponents{T}
    time::Any
    i_sc::T
    i_mp::T
    v_oc::T
    v_mp::T
    p_mp::T
    i_x::Union{Nothing,T}
    i_xx::Union{Nothing,T}
end

function rrule(::Type{DCComponents}, time, i_sc, i_mp, v_oc, v_mp, p_mp, i_x, i_xx)
    y = DCComponents(time, i_sc, i_mp, v_oc, v_mp, p_mp, i_x, i_xx)

    function pullback(ȳ)
        di_x = isnothing(i_x) ? NoTangent() : getproperty(ȳ, :i_x)
        di_xx = isnothing(i_xx) ? NoTangent() : getproperty(ȳ, :i_xx)

        return (
            NoTangent(),
            NoTangent(),
            getproperty(ȳ, :i_sc),
            getproperty(ȳ, :i_mp),
            getproperty(ȳ, :v_oc),
            getproperty(ȳ, :v_mp),
            getproperty(ȳ, :p_mp),
            di_x,
            di_xx,
        )
    end

    return y, pullback
end

"""
    read_solar_module(
        module_name="Canadian Solar CS5P-220M [ 2009]",
        module_filename="sam-library-sandia-modules-2015-6-30.csv",
        directory=joinpath(@__DIR__, "..", "data"),
        header_rows=1,
        skip_rows=4,
        T=Float64
    ) -> SolarModule{T}

Load module parameters from a SAM/Sandia module library CSV and return the entry
matching `module_name`.

# Arguments:
* `module_name`: Name of the module in the local CSV file (matched against the `Name` column)
* `module_filename`: CSV filename
* `directory`: Directory containing the CSV `module_filename`
* `header_rows`: Number of header rows in the CSV file
* `skip_rows`: Number of rows to skip at the beginning of the CSV file
* `T::Type`: Numeric type for stored module parameters

# Returns:
* `SolarModule{T}`: Solar module parameters

```jldoctest
julia> pv_module = read_solar_module("Canadian Solar CS5P-220M [ 2009]");

julia> cols = [:name,:vintage,:area,:material,:cells_in_series,:parallel_strings,:isco,:voco,:impo,:vmpo,:aisc,:aimp,:c0,:c1,:bvoco,:mbvoc,:bvmpo,:mbvmp,:n,:c2,:c3,:a0,:a1,:a2,:a3,:a4,:b0,:b1,:b2,:b3,:b4,:b5,:dtc,:fd,:a,:b,:c4,:c5,:ixo,:ixxo,:c6,:c7,:notes];

julia> fieldnames(typeof(pv_module)) == Tuple(cols)
true

julia> isapprox(pv_module.isco, 5.09115; atol=1e-6)
true
```
"""
function read_solar_module(
    module_name::AbstractString = "Canadian Solar CS5P-220M [ 2009]",
    module_filename::AbstractString = "sam-library-sandia-modules-2015-6-30.csv",
    directory::AbstractString = joinpath(@__DIR__, "..", "data"),
    header_rows::Integer = 1,
    skip_rows::Integer = 4,
    T::Type = Float64,
)

    path = joinpath(directory, module_filename)
    file = File(path; header = header_rows, skipto = skip_rows)

    for r in file
        String(r[Symbol("Name")]) == module_name || continue

        return SolarModule{T}(
            String(r[Symbol("Name")]),
            String(r[Symbol("Vintage")]),
            T(r[Symbol("Area")]),
            String(r[Symbol("Material")]),
            Int(r[Symbol("Cells in Series")]),
            Int(r[Symbol("Parallel Strings")]),
            T(r[Symbol("Isco")]),
            T(r[Symbol("Voco")]),
            T(r[Symbol("Impo")]),
            T(r[Symbol("Vmpo")]),
            T(r[Symbol("Aisc")]),
            T(r[Symbol("Aimp")]),
            T(r[Symbol("C0")]),
            T(r[Symbol("C1")]),
            T(r[Symbol("Bvoco")]),
            T(r[Symbol("Mbvoc")]),
            T(r[Symbol("Bvmpo")]),
            T(r[Symbol("Mbvmp")]),
            T(r[Symbol("N")]),
            T(r[Symbol("C2")]),
            T(r[Symbol("C3")]),
            T(r[Symbol("A0")]),
            T(r[Symbol("A1")]),
            T(r[Symbol("A2")]),
            T(r[Symbol("A3")]),
            T(r[Symbol("A4")]),
            T(r[Symbol("B0")]),
            T(r[Symbol("B1")]),
            T(r[Symbol("B2")]),
            T(r[Symbol("B3")]),
            T(r[Symbol("B4")]),
            T(r[Symbol("B5")]),
            T(r[Symbol("DTC")]),
            T(r[Symbol("FD")]),
            T(r[Symbol("A")]),
            T(r[Symbol("B")]),
            T(r[Symbol("C4")]),
            T(r[Symbol("C5")]),
            T(r[Symbol("IXO")]),
            T(r[Symbol("IXXO")]),
            T(r[Symbol("C6")]),
            T(r[Symbol("C7")]),
            String(r[Symbol("Notes")]),
        )
    end

    throw(KeyError("Module name not found: $module_name"))
end

# Internal helper: evaluate a polynomial at a scalar x.
function _polyval(coeffs::AbstractVector{<:Real}, x::Real)

    if isnan(x)
        poly_result = NaN  # Keep NaN if the input is NaN
    else
        poly_result = sum(c * x^(length(coeffs) - j) for (j, c) in enumerate(coeffs))  # Evaluate polynomial
    end

    return poly_result
end

# Internal helper: evaluate a polynomial at vector x (broadcasted via comprehension).
function _polyval(coeffs::AbstractVector{<:Real}, x::AbstractVector{<:Real})

    poly_result = [_polyval(coeffs, x[i]) for i in eachindex(x)]

    return poly_result
end

# Internal helper: compute SAPM spectral loss modifier using absolute airmass and module coefficients.
function _get_spectral_loss(pv_module::SolarModule, zenith::Real, altitude::Real)

    airmass_absolute = _get_absolute_airmass(zenith, altitude)
    airmass_coeff = [pv_module.a4, pv_module.a3, pv_module.a2, pv_module.a1, pv_module.a0]
    spectral_loss = _polyval(airmass_coeff, airmass_absolute)
    spectral_loss =
        (isnan(spectral_loss) || !isfinite(spectral_loss)) ? zero(spectral_loss) :
        pv_smooth_max(spectral_loss, zero(spectral_loss))

    return spectral_loss
end

# Internal helper: vectorized spectral loss modifier over zenith angles.
function _get_spectral_loss(
    pv_module::SolarModule,
    zenith::AbstractVector{<:Real},
    altitude::Real,
)

    spectral_loss = [_get_spectral_loss(pv_module, z, altitude) for z in zenith]

    return spectral_loss
end

# Internal helper: compute incidence angle modifier (IAM) using module coefficients.
function _incidence_angle_modifier(pv_module::SolarModule, aoi::Real)

    aoi_coeff =
        [pv_module.b5, pv_module.b4, pv_module.b3, pv_module.b2, pv_module.b1, pv_module.b0]
    iam = _polyval(aoi_coeff, aoi)
    iam = pv_smooth_max(iam, zero(iam))

    if isnan(iam) || !isfinite(iam)
        return 0.0
    end
    if aoi < 0
        return 0.0
    end

    return iam
end

# Internal helper: vectorized IAM over angles of incidence.
function _incidence_angle_modifier(pv_module::SolarModule, aoi::AbstractVector{<:Real})

    iam = [_incidence_angle_modifier(pv_module, a) for a in aoi]

    return iam
end

"""
    sapm_effective_irradiance(
        total_irradiance::TotalIrradiance,
        pv_module::SolarModule,
        solar_position::SolarPosition,
        surface_tilt::Real,
        surface_azimuth::Real,
        altitude::Real
    ) -> EffectiveIrradiance

    sapm_effective_irradiance(
        total_irradiance::Vector{<:TotalIrradiance},
        pv_module::SolarModule,
        solar_position::Vector{<:SolarPosition},
        surface_tilt::Real,
        surface_azimuth::Real,
        altitude::Real
    ) -> Vector{<:EffectiveIrradiance}

    sapm_effective_irradiance(
        total_irradiance::Vector{<:TotalIrradiance},
        pv_module::SolarModule,
        solar_position::Vector{<:SolarPosition},
        surface_tilt::Vector{<:Real},
        surface_azimuth::Real,
        altitude::Real
    ) -> Vector{<:EffectiveIrradiance}

    sapm_effective_irradiance(
        total_irradiance::Vector{<:TotalIrradiance},
        pv_module::SolarModule,
        solar_position::Vector{<:SolarPosition},
        surface_tilt::Real,
        surface_azimuth::Vector{<:Real},
        altitude::Real
    ) -> Vector{<:EffectiveIrradiance}

    sapm_effective_irradiance(
        total_irradiance::Vector{<:TotalIrradiance},
        pv_module::SolarModule,
        solar_position::Vector{<:SolarPosition},
        surface_tilt::Vector{<:Real},
        surface_azimuth::Vector{<:Real},
        altitude::Real
    ) -> Vector{<:EffectiveIrradiance}

    sapm_effective_irradiance(
        total_irradiance::Vector{<:TotalIrradiance},
        pv_module::SolarModule,
        solar_position::SolarPosition,
        surface_tilt::Real,
        surface_azimuth::Real,
        altitude::Real
    ) -> Vector{<:EffectiveIrradiance}

    sapm_effective_irradiance(
        total_irradiance::Vector{<:TotalIrradiance},
        pv_module::SolarModule,
        solar_position::SolarPosition,
        surface_tilt::Vector{<:Real},
        surface_azimuth::Real,
        altitude::Real
    ) -> Vector{<:EffectiveIrradiance}

    sapm_effective_irradiance(
        total_irradiance::Vector{<:TotalIrradiance},
        pv_module::SolarModule,
        solar_position::SolarPosition,
        surface_tilt::Real,
        surface_azimuth::Vector{<:Real},
        altitude::Real
    ) -> Vector{<:EffectiveIrradiance}

    sapm_effective_irradiance(
        total_irradiance::Vector{<:TotalIrradiance},
        pv_module::SolarModule,
        solar_position::SolarPosition,
        surface_tilt::Vector{<:Real},
        surface_azimuth::Vector{<:Real},
        altitude::Real
    ) -> Vector{<:EffectiveIrradiance}

Compute Sandia Array Performance Model (SAPM) effective irradiance from
plane-of-array (POA) irradiance components, module parameters, and solar position.

# Arguments:
* `total_irradiance`: POA irradiance components, including `poa_direct` and `poa_diffuse`
* `pv_module`: Sandia module parameters (including `a0:a4`, `b0:b5`, and `fd`)
* `solar_position`: Solar position containing at least `apparent_zenith` and `azimuth`
* `surface_tilt`: Surface tilt angle (degrees)
* `surface_azimuth`: Surface azimuth angle (degrees)
* `altitude`: Site altitude (meters)

# Returns:
* Effective irradiance value (W/m^2-equivalent scaling consistent with the SAPM formulation)

```jldoctest
julia> using TimeZones

julia> total_irradiance = TotalIrradiance(ZonedDateTime(2023, 1, 1, tz"America/Denver"), 0.0, 500.0, 100.0, 0.0, 0.0);

julia> solar_pos = SolarPosition(
           ZonedDateTime(2020, 6, 1, 12, 0, 0, tz"UTC"),
           30.0,   # apparent_zenith
           30.0,   # zenith
           60.0,   # apparent_elevation
           60.0,   # elevation
           180.0,  # azimuth
           1.0,    # equation_of_time
       );

julia> pv_module = read_solar_module("Canadian Solar CS5P-220M [ 2009]");

julia> effective_irradiance = sapm_effective_irradiance(total_irradiance, pv_module, solar_pos, 30.0, 180.0, 1600.0);

julia> isapprox(effective_irradiance.effective_irradiance, 586.1025479220853; atol=1e-6)
true
```
"""
function sapm_effective_irradiance(
    total_irradiance::TotalIrradiance,
    pv_module::SolarModule,
    solar_position::SolarPosition,
    surface_tilt::Real,
    surface_azimuth::Real,
    altitude::Real,
)

    aoi = get_angle_of_incidence(
        surface_tilt,
        surface_azimuth,
        solar_position.apparent_zenith,
        solar_position.azimuth,
    )

    spectral_loss = _get_spectral_loss(pv_module, solar_position.apparent_zenith, altitude)
    iam = _incidence_angle_modifier(pv_module, aoi)

    eff =
        spectral_loss *
        (total_irradiance.poa_direct * iam + pv_module.fd * total_irradiance.poa_diffuse)

    return EffectiveIrradiance(total_irradiance.time, eff)
end

function sapm_effective_irradiance(
    total_irradiance::Vector{<:TotalIrradiance},
    pv_module::SolarModule,
    solar_position::Vector{<:SolarPosition},
    surface_tilt::Real,
    surface_azimuth::Real,
    altitude::Real,
)

    n = length(total_irradiance)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but total_irradiance has $n rows",
        ),
    )

    return map(1:n) do ii
        sapm_effective_irradiance(
            total_irradiance[ii],
            pv_module,
            solar_position[ii],
            surface_tilt,
            surface_azimuth,
            altitude,
        )
    end
end

function sapm_effective_irradiance(
    total_irradiance::Vector{<:TotalIrradiance},
    pv_module::SolarModule,
    solar_position::Vector{<:SolarPosition},
    surface_tilt::Vector{<:Real},
    surface_azimuth::Real,
    altitude::Real,
)

    n = length(total_irradiance)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but total_irradiance has $n rows",
        ),
    )

    length(surface_tilt) == n || throw(
        DimensionMismatch(
            "surface_tilt has $(length(surface_tilt)) rows but total_irradiance has $n rows",
        ),
    )

    return map(1:n) do ii
        sapm_effective_irradiance(
            total_irradiance[ii],
            pv_module,
            solar_position[ii],
            surface_tilt[ii],
            surface_azimuth,
            altitude,
        )
    end
end

function sapm_effective_irradiance(
    total_irradiance::Vector{<:TotalIrradiance},
    pv_module::SolarModule,
    solar_position::Vector{<:SolarPosition},
    surface_tilt::Real,
    surface_azimuth::Vector{<:Real},
    altitude::Real,
)

    n = length(total_irradiance)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but total_irradiance has $n rows",
        ),
    )

    length(surface_azimuth) == n || throw(
        DimensionMismatch(
            "surface_azimuth has $(length(surface_azimuth)) rows but total_irradiance has $n rows",
        ),
    )

    return map(1:n) do ii
        sapm_effective_irradiance(
            total_irradiance[ii],
            pv_module,
            solar_position[ii],
            surface_tilt,
            surface_azimuth[ii],
            altitude,
        )
    end
end

function sapm_effective_irradiance(
    total_irradiance::Vector{<:TotalIrradiance},
    pv_module::SolarModule,
    solar_position::Vector{<:SolarPosition},
    surface_tilt::Vector{<:Real},
    surface_azimuth::Vector{<:Real},
    altitude::Real,
)

    n = length(total_irradiance)

    length(solar_position) == n || throw(
        DimensionMismatch(
            "solar_position has $(length(solar_position)) rows but total_irradiance has $n rows",
        ),
    )

    length(surface_tilt) == n || throw(
        DimensionMismatch(
            "surface_tilt has $(length(surface_tilt)) rows but total_irradiance has $n rows",
        ),
    )

    length(surface_azimuth) == n || throw(
        DimensionMismatch(
            "surface_azimuth has $(length(surface_azimuth)) rows but total_irradiance has $n rows",
        ),
    )

    return map(1:n) do ii
        sapm_effective_irradiance(
            total_irradiance[ii],
            pv_module,
            solar_position[ii],
            surface_tilt[ii],
            surface_azimuth[ii],
            altitude,
        )
    end
end

function sapm_effective_irradiance(
    total_irradiance::Vector{<:TotalIrradiance},
    pv_module::SolarModule,
    solar_position::SolarPosition,
    surface_tilt::Real,
    surface_azimuth::Real,
    altitude::Real,
)

    n = length(total_irradiance)

    return map(1:n) do ii
        sapm_effective_irradiance(
            total_irradiance[ii],
            pv_module,
            solar_position,
            surface_tilt,
            surface_azimuth,
            altitude,
        )
    end
end

function sapm_effective_irradiance(
    total_irradiance::Vector{<:TotalIrradiance},
    pv_module::SolarModule,
    solar_position::SolarPosition,
    surface_tilt::Vector{<:Real},
    surface_azimuth::Real,
    altitude::Real,
)

    n = length(total_irradiance)

    return map(1:n) do ii
        sapm_effective_irradiance(
            total_irradiance[ii],
            pv_module,
            solar_position,
            surface_tilt[ii],
            surface_azimuth,
            altitude,
        )
    end
end

function sapm_effective_irradiance(
    total_irradiance::Vector{<:TotalIrradiance},
    pv_module::SolarModule,
    solar_position::SolarPosition,
    surface_tilt::Real,
    surface_azimuth::Vector{<:Real},
    altitude::Real,
)

    n = length(total_irradiance)

    return map(1:n) do ii
        sapm_effective_irradiance(
            total_irradiance[ii],
            pv_module,
            solar_position,
            surface_tilt,
            surface_azimuth[ii],
            altitude,
        )
    end
end

function sapm_effective_irradiance(
    total_irradiance::Vector{<:TotalIrradiance},
    pv_module::SolarModule,
    solar_position::SolarPosition,
    surface_tilt::Vector{<:Real},
    surface_azimuth::Vector{<:Real},
    altitude::Real,
)

    n = length(total_irradiance)

    return map(1:n) do ii
        sapm_effective_irradiance(
            total_irradiance[ii],
            pv_module,
            solar_position,
            surface_tilt[ii],
            surface_azimuth[ii],
            altitude,
        )
    end
end

"""
    sapm_dc_components(
        pv_module::SolarModule,
        effective_irradiance::EffectiveIrradiance,
        cell_temperature::CellTemperature,
        temperature_ref::Real=25.0,
        irradiance_ref::Real=1000.0,
        q::Real=1.602176634e-19,
        kb::Real=1.380649e-23
    ) -> DCComponents

    sapm_dc_components(
        pv_module::SolarModule,
        effective_irradiance::Vector{<:EffectiveIrradiance},
        cell_temperature::Vector{<:CellTemperature},
        temperature_ref::Real=25.0,
        irradiance_ref::Real=1000.0,
        q::Real=1.602176634e-19,
        kb::Real=1.380649e-23
    ) -> DCComponents

Compute IV curve summary points using the Sandia Array Performance Model (SAPM).

# Arguments:
* `pv_module`: Sandia module parameters
* `effective_irradiance`: Effective irradiance (typically from `sapm_effective_irradiance`)
* `cell_temperature`: Cell temperature (°C)
* `temperature_ref`: Reference cell temperature (°C)
* `irradiance_ref`: Reference irradiance (W/m^2)
* `q`: Elementary charge (C)
* `kb`: Boltzmann constant (J/K)

# Returns:
DC power components of the PV module in the given irradiance conditions

```jldoctest
julia> using TimeZones

julia> pv_module = read_solar_module("Canadian Solar CS5P-220M [ 2009]");

julia> effective_irradiance = EffectiveIrradiance(ZonedDateTime(2023, 1, 1, tz"America/Denver"), 800);  # effective irradiance (W/m^2)

julia> cell_temp = CellTemperature(ZonedDateTime(2023, 1, 1, tz"America/Denver"), 25);   # cell temperature (C)

julia> dc_components = sapm_dc_components(pv_module, effective_irradiance, cell_temp);

julia> fieldnames(typeof(dc_components)) == (:time,:i_sc,:i_mp,:v_oc,:v_mp,:p_mp,:i_x,:i_xx)
true

julia> dc_components.i_sc
4.07292
```
"""
function sapm_dc_components(
    pv_module::SolarModule,
    effective_irradiance::EffectiveIrradiance,
    cell_temperature::CellTemperature,
    temperature_ref::Real = 25.0,
    irradiance_ref::Real = 1000.0,
    q::Real = 1.602176634e-19,
    kb::Real = 1.380649e-23,
)

    effective_irradiance.time == cell_temperature.time || throw(
        ArgumentError(
            "Time mismatch: effective_irradiance=$(effective_irradiance.time) cell_temperature=$(cell_temperature.time)",
        ),
    )

    Ee = effective_irradiance.effective_irradiance / irradiance_ref
    Ee_pos = pv_smooth_max(Ee, oftype(Ee, 1.0e-9))

    # Bvmpo/Bvoco terms
    Bvmpo = pv_module.bvmpo + pv_module.mbvmp * (1 - Ee)
    Bvoco = pv_module.bvoco + pv_module.mbvoc * (1 - Ee)

    delta = pv_module.n * kb * (cell_temperature.cell_temperature + 273.15) / q

    logEe = log(Ee_pos)

    Ns = pv_module.cells_in_series

    i_sc =
        pv_module.isco *
        Ee *
        (1 + pv_module.aisc * (cell_temperature.cell_temperature - temperature_ref))
    i_mp =
        pv_module.impo *
        (pv_module.c0 * Ee + pv_module.c1 * Ee^2) *
        (1 + pv_module.aimp * (cell_temperature.cell_temperature - temperature_ref))

    v_oc = pv_smooth_max(
        pv_module.voco +
        Ns * delta * logEe +
        Bvoco * (cell_temperature.cell_temperature - temperature_ref),
        zero(Ee),
    )
    v_mp = pv_smooth_max(
        pv_module.vmpo +
        pv_module.c2 * Ns * delta * logEe +
        pv_module.c3 * Ns * (delta * logEe)^2 +
        Bvmpo * (cell_temperature.cell_temperature - temperature_ref),
        zero(Ee),
    )

    p_mp = i_mp * v_mp

    # Optional IXO/IXXO parts
    i_x =
        (
            pv_module.ixo === missing ||
            pv_module.c4 === missing ||
            pv_module.c5 === missing
        ) ? nothing :
        pv_module.ixo *
        (pv_module.c4 * Ee + pv_module.c5 * Ee^2) *
        (1 + pv_module.aisc * (cell_temperature.cell_temperature - temperature_ref))

    i_xx =
        (
            pv_module.ixxo === missing ||
            pv_module.c6 === missing ||
            pv_module.c7 === missing
        ) ? nothing :
        pv_module.ixxo *
        (pv_module.c6 * Ee + pv_module.c7 * Ee^2) *
        (1 + pv_module.aimp * (cell_temperature.cell_temperature - temperature_ref))

    T = promote_type(typeof(i_sc), typeof(i_mp), typeof(v_oc), typeof(v_mp), typeof(p_mp))

    return DCComponents(effective_irradiance.time, i_sc, i_mp, v_oc, v_mp, p_mp, i_x, i_xx)
end

function sapm_dc_components(
    pv_module::SolarModule,
    effective_irradiance::Vector{<:EffectiveIrradiance},
    cell_temperature::Vector{<:CellTemperature},
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

    return map(1:n) do ii
        sapm_dc_components(
            pv_module,
            effective_irradiance[ii],
            cell_temperature[ii],
            temperature_ref,
            irradiance_ref,
            q,
            kb,
        )
    end
end

# SolarModule pretty printer (shared implementation)
function _show_module_pretty(io::IO, sm::SolarModule)
    println(io, "SolarModule")
    println(io, "────────────")
    println(io, "name               │ ", sm.name)
    println(io, "vintage            │ ", sm.vintage)
    println(io, "area               │ ", sm.area)
    println(io, "material           │ ", sm.material)
    println(io, "cells_in_series    │ ", sm.cells_in_series)
    println(io, "parallel_strings   │ ", sm.parallel_strings)
    println(io, "isco               │ ", sm.isco)
    println(io, "voco               │ ", sm.voco)
    println(io, "impo               │ ", sm.impo)
    println(io, "vmpo               │ ", sm.vmpo)
    println(io, "aisc               │ ", sm.aisc)
    println(io, "aimp               │ ", sm.aimp)
    println(io, "c0                 │ ", sm.c0)
    println(io, "c1                 │ ", sm.c1)
    println(io, "bvoco              │ ", sm.bvoco)
    println(io, "mbvoc              │ ", sm.mbvoc)
    println(io, "bvmpo              │ ", sm.bvmpo)
    println(io, "mbvmp              │ ", sm.mbvmp)
    println(io, "n                  │ ", sm.n)
    println(io, "c2                 │ ", sm.c2)
    println(io, "c3                 │ ", sm.c3)
    println(io, "a0                 │ ", sm.a0)
    println(io, "a1                 │ ", sm.a1)
    println(io, "a2                 │ ", sm.a2)
    println(io, "a3                 │ ", sm.a3)
    println(io, "a4                 │ ", sm.a4)
    println(io, "b0                 │ ", sm.b0)
    println(io, "b1                 │ ", sm.b1)
    println(io, "b2                 │ ", sm.b2)
    println(io, "b3                 │ ", sm.b3)
    println(io, "b4                 │ ", sm.b4)
    println(io, "b5                 │ ", sm.b5)
    println(io, "dtc                │ ", sm.dtc)
    println(io, "fd                 │ ", sm.fd)
    println(io, "a                  │ ", sm.a)
    println(io, "b                  │ ", sm.b)
    println(io, "c4                 │ ", sm.c4)
    println(io, "c5                 │ ", sm.c5)
    println(io, "ixo                │ ", sm.ixo)
    println(io, "ixxo               │ ", sm.ixxo)
    println(io, "c6                 │ ", sm.c6)
    println(io, "c7                 │ ", sm.c7)
    println(io, "notes              │ ", sm.notes)
end

function Base.show(io::IO, sm::SolarModule)
    _show_module_pretty(io, sm)
end

function Base.show(io::IO, ::MIME"text/plain", sm::SolarModule)
    _show_module_pretty(io, sm)
end

# EffectiveIrradiance pretty printer (shared implementation)
function _show_effective_irradiance_pretty(io::IO, ei::EffectiveIrradiance)
    println(io, "EffectiveIrradiance")
    println(io, "──────────────────")
    println(io, "time                │ ", ei.time)
    println(io, "effective_irradiance │ ", ei.effective_irradiance)
end

function Base.show(io::IO, ei::EffectiveIrradiance)
    _show_effective_irradiance_pretty(io, ei)
end

function Base.show(io::IO, ::MIME"text/plain", ei::EffectiveIrradiance)
    _show_effective_irradiance_pretty(io, ei)
end

# EffectiveIrradiance vector pretty printer (shared implementation)
function _show_effective_irradiance_vec_pretty(io::IO, v::Vector{<:EffectiveIrradiance})
    header = ["time", "effective_irradiance"]
    _show_table(
        io,
        "EffectiveIrradiance",
        header,
        i -> begin
            x = v[i]
            [string(x.time), string(x.effective_irradiance)]
        end,
        length(v),
    )
end

function Base.show(io::IO, v::Vector{<:EffectiveIrradiance})
    _show_effective_irradiance_vec_pretty(io, v)
end

function Base.show(io::IO, ::MIME"text/plain", v::Vector{<:EffectiveIrradiance})
    _show_effective_irradiance_vec_pretty(io, v)
end

# DCComponents pretty printer (shared implementation)
function _show_dc_pretty(io::IO, x::DCComponents)
    println(io, "DC Components (SAPM)")
    println(io, "──────")
    println(io, "time │ ", x.time)
    println(io, "i_sc │ ", x.i_sc)
    println(io, "i_mp │ ", x.i_mp)
    println(io, "v_oc │ ", x.v_oc)
    println(io, "v_mp │ ", x.v_mp)
    println(io, "p_mp │ ", x.p_mp)
    println(io, "i_x  │ ", x.i_x)
    println(io, "i_xx │ ", x.i_xx)
end

function Base.show(io::IO, x::DCComponents)
    _show_dc_pretty(io, x)
end

function Base.show(io::IO, ::MIME"text/plain", x::DCComponents)
    _show_dc_pretty(io, x)
end

# DCComponents vector pretty printer (shared implementation)
function _show_dc_vec_pretty(io::IO, v::Vector{<:DCComponents})
    header = ["time", "i_sc", "i_mp", "v_oc", "v_mp", "p_mp", "i_x", "i_xx"]
    _show_table(
        io,
        "DCComponents",
        header,
        i -> begin
            x = v[i]
            [
                string(x.time),
                string(x.i_sc),
                string(x.i_mp),
                string(x.v_oc),
                string(x.v_mp),
                string(x.p_mp),
                string(x.i_x),
                string(x.i_xx),
            ]
        end,
        length(v),
    )
end

# Pretty table for a vector of rows
function Base.show(io::IO, v::Vector{<:DCComponents})
    _show_dc_vec_pretty(io, v)
end

function Base.show(io::IO, ::MIME"text/plain", v::Vector{<:DCComponents})
    _show_dc_vec_pretty(io, v)
end
