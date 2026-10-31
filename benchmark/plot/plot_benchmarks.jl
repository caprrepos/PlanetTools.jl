## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-27 to 2026-10-31

# Plot AirspeedVelocity.jl benchmark results of PlanetTools: run time, run time relative to
# a baseline revision, and allocated memory, one row per benchmark. Also writes the
# numbers as a Markdown table next to the figure.
#
# Usage (from the repository root):
#     julia --project=benchmark/plot -e 'using Pkg; Pkg.instantiate()'
#     julia --project=benchmark/plot benchmark/plot/plot_benchmarks.jl [options] inputs...
#
# Inputs are result files `results_PlanetTools@<rev>.json` written by `benchpkg`, or
# directories containing them. A revision is labeled by the `<rev>` part of its file name;
# write `label=path` to choose another label, e.g. to compare thread counts:
#     ... plot_benchmarks.jl "1 thread=t1/results_PlanetTools@dirty.json" \
#         "4 threads=t4/results_PlanetTools@dirty.json"
#
# Options:
#     --output=<file>   figure file, `.png`, `.svg`, or `.pdf` (default
#                       `benchmark_results.png`); the table is written to `<file>.md`.
#     --rev=<a,b,...>   revisions to show, in this order (default: all, in input order).
#     --baseline=<rev>  revision the others are compared with (default: the first).
#     --filter=<text>   show only benchmarks whose name contains `text`.
using CairoMakie
using JSON: JSON
using Printf: @sprintf
using Statistics: median, quantile

const PACKAGE = "PlanetTools"

# Colors of the revisions in fixed order (validated for color-vision deficiencies); each
# revision also gets its own marker shape, so that color is never the only distinction.
const COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
const MARKERS = [:circle, :rect, :diamond, :utriangle, :dtriangle, :xcross]
const INK, MUTED, GRID = "#222222", "#6b6b6b", "#e6e6e6"

"""
    Result

Benchmark result of one benchmark and revision: sample times in ns, allocated memory in
bytes, and number of allocations.
"""
struct Result
    times::Vector{Float64}
    memory::Float64
    allocs::Int
end

# Leaves of the (possibly nested) result dictionary, keyed by "group / … / name".
function flatten_results(d, prefix="")
    out = Dict{String,Result}()
    for (k, v) in d
        name = isempty(prefix) ? string(k) : string(prefix, " / ", k)
        if v isa AbstractDict && haskey(v, "times")
            out[name] = Result(Float64.(v["times"]), v["memory"], v["allocs"])
        elseif v isa AbstractDict && haskey(v, "data") && haskey(v, "tags")
            merge!(out, flatten_results(v["data"], name))  # A nested BenchmarkGroup.
        elseif v isa AbstractDict
            merge!(out, flatten_results(v, name))
        end
    end
    return out
end

function read_results(path)
    d = JSON.parsefile(path)
    return flatten_results(haskey(d, "data") ? d["data"] : d)
end

# Label from `results_PlanetTools@<rev>.json`; full commit hashes are shortened to 7 digits.
function rev_label(path)
    m = match(Regex("^results_$(PACKAGE)@(.+)\\.json\$"), basename(path))
    label = m === nothing ? splitext(basename(path))[1] : m.captures[1]
    return occursin(r"^[0-9a-f]{40}$", label) ? label[1:7] : label
end

function collect_inputs(args)
    inputs = Pair{String,String}[]
    for a in args
        label, path = occursin('=', a) ? split(a, '='; limit=2) : ("", a)
        files = if isdir(path)
            sort(
                filter(
                    f -> occursin(Regex("^results_$(PACKAGE)@.+\\.json\$"), f),
                    readdir(path),
                ),
            )
        else
            [basename(path)]
        end
        dir = isdir(path) ? path : dirname(path)
        isempty(files) && error("no result files in $path")
        for f in files
            push!(inputs, (isempty(label) ? rev_label(f) : label) => joinpath(dir, f))
        end
    end
    return inputs
end

function parse_args(args)
    opts = Dict(
        "output" => "benchmark_results.png", "rev" => "", "baseline" => "", "filter" => ""
    )
    positional = String[]
    for a in args
        m = match(r"^--(\w+)=(.*)$", a)
        if m === nothing
            push!(positional, a)
        elseif haskey(opts, m.captures[1])
            opts[m.captures[1]] = m.captures[2]
        else
            error("unknown option --$(m.captures[1])")
        end
    end
    isempty(positional) && error("no result files given; see the header of this script")
    return opts, positional
end

fmt_time(ns) =
    if ns < 1e3
        @sprintf("%.3g ns", ns)
    elseif ns < 1e6
        @sprintf("%.3g μs", ns / 1e3)
    elseif ns < 1e9
        @sprintf("%.3g ms", ns / 1e6)
    else
        @sprintf("%.3g s", ns / 1e9)
    end

fmt_bytes(b) =
    if b < 1024
        @sprintf("%.3g B", b)
    elseif b < 1024^2
        @sprintf("%.3g KiB", b / 1024)
    elseif b < 1024^3
        @sprintf("%.3g MiB", b / 1024^2)
    else
        @sprintf("%.3g GiB", b / 1024^3)
    end

function main(args)
    opts, positional = parse_args(args)
    inputs = collect_inputs(positional)
    if !isempty(opts["rev"])
        wanted = split(opts["rev"], ',')
        inputs = [first(filter(p -> p.first == r, inputs)) for r in wanted]
    end
    labels = first.(inputs)
    results = Dict(l => read_results(p) for (l, p) in inputs)
    baseline = isempty(opts["baseline"]) ? labels[1] : opts["baseline"]
    baseline in labels || error("baseline $baseline is not among the revisions $labels")
    # Without --rev, the baseline comes first.
    isempty(opts["rev"]) && (labels = [baseline; filter(!=(baseline), labels)])

    # Benchmarks present in all revisions, sorted by name; package load time last.
    names = sort(collect(intersect((keys(results[l]) for l in labels)...)))
    isempty(opts["filter"]) || filter!(n -> occursin(opts["filter"], n), names)
    i = findfirst(==("time_to_load"), names)
    i === nothing || push!(deleteat!(names, i), "time_to_load")
    isempty(names) && error("no benchmark is common to all revisions")
    nb, nr = length(names), length(labels)

    # One row per benchmark (top to bottom); revisions are offset within a row.
    # The first revision is drawn at the top of a row, as in the legend.
    ypos(ib, ir) = (nb - ib + 1) - 0.5 * ((ir - 1) / max(nr - 1, 1) - 0.5) * (nr > 1)
    ylabels = [n == "time_to_load" ? "package load time" : n for n in names]

    fig = Figure(;
        size=(1500, 120 + 34 * nb * max(1, nr / 2)), fontsize=14, backgroundcolor=:white
    )
    kw = (;
        xgridcolor=GRID,
        ygridvisible=false,
        xminorgridvisible=false,
        leftspinevisible=false,
        rightspinevisible=false,
        topspinevisible=false,
        bottomspinecolor=MUTED,
        xtickcolor=MUTED,
        ytickcolor=MUTED,
        xticklabelcolor=INK,
        yticklabelcolor=INK,
        titlealign=:left,
        titlecolor=INK,
        titlesize=16,
        yticks=(nb:-1:1, ylabels),
        xlabelcolor=MUTED,
    )
    ax_t = Axis(
        fig[1, 1];
        title="Run time (median; bar: 25–75 % of samples)",
        xscale=log10,
        xlabel="time per call",
        xtickformat=v -> fmt_time.(v .* 1e9),
        kw...,
    )
    ax_r = Axis(
        fig[1, 2];
        title="Minimum run time relative to $baseline",
        xscale=log2,
        xlabel="ratio (< 1 faster, > 1 slower)",
        xtickformat=v -> [@sprintf("%.3g×", x) for x in v],
        kw...,
    )
    ax_m = Axis(
        fig[1, 3];
        title="Allocated memory per call",
        xscale=log10,
        xlabel="memory per call",
        xtickformat=v -> fmt_bytes.(v),
        kw...,
    )
    hideydecorations!(ax_r; grid=false)
    hideydecorations!(ax_m; grid=false)
    linkyaxes!(ax_t, ax_r, ax_m)
    colsize!(fig.layout, 1, Relative(0.45))

    # Light bands behind every other row, and the band of ±10 % around a ratio of 1,
    # within which differences are usually measurement noise.
    for ax in (ax_t, ax_r, ax_m), ib in 2:2:nb
        hspan!(ax, nb - ib + 0.5, nb - ib + 1.5; color=(:black, 0.035))
    end
    vspan!(ax_r, 1 / 1.1, 1.1; color=(:black, 0.06))
    vlines!(ax_r, 1; color=MUTED, linewidth=1)

    for (ir, l) in enumerate(labels)
        c, mk = COLORS[mod1(ir, length(COLORS))], MARKERS[mod1(ir, length(MARKERS))]
        for (ib, n) in enumerate(names)
            r, y = results[l][n], ypos(ib, ir)
            t = r.times ./ 1e9
            lines!(
                ax_t, [quantile(t, 0.25), quantile(t, 0.75)], [y, y]; color=c, linewidth=2
            )
            scatter!(
                ax_t,
                [median(t)],
                [y];
                color=c,
                marker=mk,
                markersize=10,
                strokecolor=:white,
                strokewidth=1,
            )
            ratio = minimum(r.times) / minimum(results[baseline][n].times)
            l == baseline || scatter!(
                ax_r,
                [ratio],
                [y];
                color=c,
                marker=mk,
                markersize=10,
                strokecolor=:white,
                strokewidth=1,
            )
            scatter!(
                ax_m,
                [max(r.memory, 1.0)],
                [y];
                color=c,
                marker=mk,
                markersize=10,
                strokecolor=:white,
                strokewidth=1,
            )
        end
    end
    # Symmetric ratio axis that always shows the ±10 % band.
    ratios = [
        minimum(results[l][n].times) / minimum(results[baseline][n].times) for
        l in labels, n in names
    ]
    rmax = max(1.25, maximum(max.(ratios, 1 ./ ratios)) * 1.15)
    xlims!(ax_r, 1 / rmax, rmax)
    ax_r.xticks = filter(
        x -> 1 / rmax <= x <= rmax, [1 / 8, 1 / 4, 1 / 2, 0.8, 1, 1.25, 2, 4, 8]
    )
    ax_m.xticks = [1024.0^k for k in 0:4]  # 1 B, 1 KiB, 1 MiB, 1 GiB, 1 TiB.
    ylims!(ax_t, 0.4, nb + 0.6)

    elems = [
        MarkerElement(;
            color=COLORS[mod1(i, end)], marker=MARKERS[mod1(i, end)], markersize=12
        ) for i in eachindex(labels)
    ]
    Legend(
        fig[0, 1:3],
        elems,
        labels,
        "Revision";
        orientation=:horizontal,
        framevisible=false,
        titleposition=:left,
        labelcolor=INK,
        titlecolor=INK,
    )

    out = opts["output"]
    save(out, fig)
    write_table(out * ".md", names, labels, results, baseline)
    println("Wrote $out and $out.md")
    return nothing
end

function write_table(path, names, labels, results, baseline)
    open(path, "w") do io
        println(
            io,
            "| Benchmark | Revision | Minimum | Median | Ratio to $baseline ",
            "(minimum) | Memory | Allocations |",
        )
        println(io, "|:--|:--|--:|--:|--:|--:|--:|")
        for n in names, l in labels
            r = results[l][n]
            ratio = minimum(r.times) / minimum(results[baseline][n].times)
            println(
                io,
                "| ",
                n,
                " | ",
                l,
                " | ",
                fmt_time(minimum(r.times)),
                " | ",
                fmt_time(median(r.times)),
                " | ",
                @sprintf("%.2f", ratio),
                " | ",
                fmt_bytes(r.memory),
                " | ",
                r.allocs,
                " |",
            )
        end
    end
    return nothing
end

main(ARGS)
