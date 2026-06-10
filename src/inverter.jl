"""
    SolarInverter{T}

Container for CEC/SAM inverter parameters.

These parameters are typically loaded from the CEC/SAM inverter library CSV via
`read_solar_inverter`.

# Fields:
* `name`: Inverter name
* `vac`: Nominal AC voltage (V)
* `pso`: Self-consumption / startup power (W)
* `paco`: Rated AC output power (W)
* `pdco`: Rated DC input power (W)
* `vdco`: DC voltage at rated power (V)
* `c0`, `c1`, `c2`, `c3`: Sandia inverter model coefficients
* `pnt`: Night tare / standby power (W)
* `vdcmax`: Maximum DC voltage (V)
* `idcmax`: Maximum DC current (A)
* `mppt_low`: MPPT minimum voltage (V)
* `mppt_high`: MPPT maximum voltage (V)
* `cec_date`: CEC library date string
* `cec_type`: CEC type string
"""
struct SolarInverter{T}
    name::String
    vac::T
    pso::T
    paco::T
    pdco::T
    vdco::T
    c0::T
    c1::T
    c2::T
    c3::T
    pnt::T
    vdcmax::T
    idcmax::T
    mppt_low::T
    mppt_high::T
    cec_date::String
    cec_type::String
end

"""
Sandia inverter-model AC power for a single timestamp.

Typically computed from DC operating point(s) and an inverter model using
`sandia_ac_power`.

# Fields:
* `time`: Timestamp
* `ac_power`: AC power output (W)
"""
struct ACPower{T}
    time::Any
    ac_power::T
end

function rrule(::Type{ACPower}, time, ac_power)
    y = ACPower{typeof(ac_power)}(time, ac_power)

    function pullback(ȳ)
        return NoTangent(), NoTangent(), getproperty(ȳ, :ac_power)
    end

    return y, pullback
end

"""
    read_solar_inverter(
        inverter_name="ABB: MICRO-0.25-I-OUTD-US-208 [208V]",
        inverter_filename="sam-library-cec-inverters-2019-03-05.csv",
        directory=joinpath(@__DIR__, "..", "data"),
        header_rows=1, skip_rows=4,
        T=Float64
    ) -> SolarInverter{T}

Load inverter parameters from a CEC/SAM inverter library CSV and
return the entry matching `inverter_name`.

# Arguments:
* `inverter_name`: Name of the inverter in the local CSV file
* `inverter_filename`: CSV filename
* `directory`: Directory containing the CSV `inverter_filename`
* `header_rows::Integer`: Number of header rows in the CSV file
* `skip_rows::Integer`: Number of rows to skip at the beginning of the CSV file
* `T::Type`: Numeric type for stored inverter parameters

# Returns:
* `SolarInverter{T}`: Solar inverter parameters

```jldoctest
julia> pv_inverter = read_solar_inverter("ABB: MICRO-0.25-I-OUTD-US-208 [208V]");

julia> cols = [:name,:vac,:pso,:paco,:pdco,:vdco,:c0,:c1,:c2,:c3,:pnt,:vdcmax,:idcmax,:mppt_low,:mppt_high,:cec_date,:cec_type];

julia> fieldnames(typeof(pv_inverter)) == Tuple(cols)
true

julia> pv_inverter.pdco
259.588593
```
"""
function read_solar_inverter(
    inverter_name::AbstractString = "ABB: MICRO-0.25-I-OUTD-US-208 [208V]",
    inverter_filename::AbstractString = "sam-library-cec-inverters-2019-03-05.csv",
    directory::AbstractString = joinpath(@__DIR__, "..", "data"),
    header_rows::Integer = 1,
    skip_rows::Integer = 4,
    T::Type = Float64,
)

    path = joinpath(directory, inverter_filename)
    file = File(path; header = header_rows, skipto = skip_rows)

    for r in file
        String(r[Symbol("Name")]) == inverter_name || continue

        cec_date_raw = r[Symbol("CEC_Date")]
        cec_date = begin
            s = strip(String(cec_date_raw))
            (isempty(s) || lowercase(s) in ("n/a", "na", "null")) ? missing : s
        end

        return SolarInverter{T}(
            String(r[Symbol("Name")]),
            Base.parse(T, String(r[Symbol("Vac")])),
            T(r[Symbol("Pso")]),
            T(r[Symbol("Paco")]),
            T(r[Symbol("Pdco")]),
            T(r[Symbol("Vdco")]),
            T(r[Symbol("C0")]),
            T(r[Symbol("C1")]),
            T(r[Symbol("C2")]),
            T(r[Symbol("C3")]),
            T(r[Symbol("Pnt")]),
            T(r[Symbol("Vdcmax")]),
            T(r[Symbol("Idcmax")]),
            T(r[Symbol("Mppt_low")]),
            T(r[Symbol("Mppt_high")]),
            String(r[Symbol("CEC_Date")]),
            String(r[Symbol("CEC_Type")]),
        )
    end

    throw(KeyError("Inverter name not found: $inverter_name"))
end

"""
    sandia_ac_power(
        inverter, dc_power
    ) -> ac_power

Compute the AC power from DC operating point(s) using the Sandia
inverter model.

# Arguments:
* `inverter::SolarInverter`: Inverter properties
* `dc_power::DCComponents` or `AbstractVector{<:DCComponents}`: DC operating point(s)
  for a photovoltaic system over time

# Returns:
* `ac_power::ACPower` or `Vector{<:ACPower}`: AC power over time (W)

```jldoctest
julia> using TimeZones

julia> dc_power = DCComponents(
           ZonedDateTime(2023, 1, 1, tz"America/Denver"),
           5.0, # i_sc (A)
           2.0, # i_mp (A)
           12.0, # v_oc (V)
           10.0, # v_mp (V)
           20.0, # p_mp (W)
           nothing, # i_x (A)
           nothing, # v_x (V)
       );

julia> pv_inverter = read_solar_inverter("ABB: MICRO-0.25-I-OUTD-US-208 [208V]");

julia> ac_power = sandia_ac_power(pv_inverter, dc_power);

julia> isapprox(ac_power.ac_power, 17.615580191975287; atol=1e-6)
true
```
"""
function sandia_ac_power(pv_inverter::SolarInverter, dc_components::DCComponents)

    A = pv_inverter.pdco * (1 + pv_inverter.c1 * (dc_components.v_mp - pv_inverter.vdco))
    B = pv_inverter.pso * (1 + pv_inverter.c2 * (dc_components.v_mp - pv_inverter.vdco))
    C = pv_inverter.c0 * (1 + pv_inverter.c3 * (dc_components.v_mp - pv_inverter.vdco))

    ac_model =
        (pv_inverter.paco / (A - B) - C * (A - B)) * (dc_components.p_mp - B) +
        C * (dc_components.p_mp - B)^2
    ac_limited = pv_smooth_min(ac_model, pv_inverter.paco)
    standby = -pv_smooth_abs(pv_inverter.pnt)
    startup_gate = pv_smooth_step(dc_components.p_mp - pv_inverter.pso)
    ac = startup_gate * ac_limited + (one(startup_gate) - startup_gate) * standby

    T = promote_type(typeof(ac))
    return ACPower{T}(dc_components.time, T(ac))
end

function sandia_ac_power(
    pv_inverter::SolarInverter,
    dc_components::AbstractVector{<:DCComponents},
)

    n = length(dc_components)

    return map(1:n) do ii
        sandia_ac_power(pv_inverter, dc_components[ii])
    end
end

# SolarInverter pretty printer (shared implementation)
function _show_inverter_pretty(io::IO, inv::SolarInverter)
    println(io, "SolarInverter")
    println(io, "────────────")
    println(io, "name       │ ", inv.name)
    println(io, "vac        │ ", inv.vac)
    println(io, "pso        │ ", inv.pso)
    println(io, "paco       │ ", inv.paco)
    println(io, "pdco       │ ", inv.pdco)
    println(io, "vdco       │ ", inv.vdco)
    println(io, "c0         │ ", inv.c0)
    println(io, "c1         │ ", inv.c1)
    println(io, "c2         │ ", inv.c2)
    println(io, "c3         │ ", inv.c3)
    println(io, "pnt        │ ", inv.pnt)
    println(io, "vdcmax     │ ", inv.vdcmax)
    println(io, "idcmax     │ ", inv.idcmax)
    println(io, "mppt_low   │ ", inv.mppt_low)
    println(io, "mppt_high  │ ", inv.mppt_high)
    println(io, "cec_date   │ ", inv.cec_date)
    println(io, "cec_type   │ ", inv.cec_type)
end

function Base.show(io::IO, inv::SolarInverter)
    _show_inverter_pretty(io, inv)
end

function Base.show(io::IO, ::MIME"text/plain", inv::SolarInverter)
    _show_inverter_pretty(io, inv)
end

# ACPower pretty printer (shared implementation)
function _show_ac_pretty(io::IO, ac::ACPower)
    println(io, "ACPower")
    println(io, "───────")
    println(io, "time     │ ", ac.time)
    println(io, "ac_power │ ", ac.ac_power)
end

function Base.show(io::IO, x::ACPower)
    _show_ac_pretty(io, x)
end

function Base.show(io::IO, ::MIME"text/plain", x::ACPower)
    _show_ac_pretty(io, x)
end

# ACPower vector pretty printer (shared implementation)
function _show_ac_vec_pretty(io::IO, ac::Vector{<:ACPower})
    header = ["time", "ac_power"]
    _show_table(
        io,
        "ACPower",
        header,
        i -> begin
            x = ac[i]
            [string(x.time), string(x.ac_power)]
        end,
        length(ac),
    )
end

function Base.show(io::IO, v::Vector{<:ACPower})
    _show_ac_vec_pretty(io, v)
end

function Base.show(io::IO, ::MIME"text/plain", v::Vector{<:ACPower})
    _show_ac_vec_pretty(io, v)
end
