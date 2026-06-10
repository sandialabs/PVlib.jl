# Function to display a table for timeseries vectors.
function _show_table(
    io::IO,
    title::AbstractString,
    header::Vector{String},
    rowfun::Function,
    n::Int;
    maxrows::Int = 10,
)

    println(io, title, " × ", n)
    n == 0 && return

    k = min(n, maxrows)

    # compute widths using header + only printed rows
    w = [length(h) for h in header]
    cached = Vector{Vector{String}}(undef, k)
    for i = 1:k
        r = rowfun(i)
        cached[i] = r
        for j in eachindex(w)
            w[j] = max(w[j], length(r[j]))
        end
    end

    printrow(cells) =
        (println(io, join((rpad(cells[j], w[j]) for j in eachindex(cells)), "  ")))

    printrow(header)
    printrow([repeat("─", wj) for wj in w])
    for i = 1:k
        printrow(cached[i])
    end
    if n > maxrows
        println(io, "… (", n - maxrows, " more rows)")
    end
end

# Recipe for plotting timeseries vectors.
@recipe function timeseries_plot(v::AbstractVector, field::Symbol, xmode::Symbol = :time)
    isempty(v) && throw(ArgumentError("Cannot plot an empty vector."))

    # Require a `time` field
    hasfield(typeof(first(v)), :time) ||
        throw(ArgumentError("Elements must have a `time` field."))

    # Require requested field
    hasfield(typeof(first(v)), field) ||
        throw(ArgumentError("Unknown field $(field) for element type $(typeof(first(v)))."))

    field == :time && throw(ArgumentError("Choose a field other than :time"))

    x = getfield.(v, :time)
    y = getfield.(v, field)

    if xmode == :time
        x = x
        xguide --> "time"
        xrotation --> -30 # rotate x tick labels
    elseif xmode == :elapsed
        t0 = first(x)
        x = (x .- t0) ./ Millisecond(1) ./ 1_000  # seconds since start (Float64)
        xguide --> "time (s)"
    end

    yguide --> String(field)
    label --> String(field)

    x, y
end
