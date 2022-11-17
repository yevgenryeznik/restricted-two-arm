# A simulation study done for the publication:
# -   title: A roadmap to using randomization in clinical trials
# - authors: Vance W. Berger, Louis Joseph Bour, Kerstine Carter, Jonathan J. Chipman, Colin C. Everett, Nicole Heussen, 
#            Catherine Hewitt, Ralf-Dieter Hilgers, Yuqun Abigail Luo, Jone Renteria, Yevgen Ryeznik, Oleksandr Sverdlov & 
#            Diane Uschner
# - journal: BMC Medical Research Methodology 
# -  volume: 21
# -  number: 168
# -    year: 2021

using ColorSchemes
using CSV
using DataFrames
using Distributions
using Gadfly
using HypothesisTests
using LaTeXStrings
using LinearAlgebra
using Pipe
using Random: seed!
using StatsBase
import Cairo, Fontconfig

include("simulated-rct.jl")
include("types.jl")
include("set-label.jl")
include("probabilities.jl")
include("simulation.jl")
include("calculate-op.jl")
include("generate-rsp.jl")

# sample size
nsbj = 50

# number of simulations
nsim = 10000

# seed
seed = 3141592

# randomization procedures
rnd = [
    RAND(nsbj), TBD(nsbj), PBD(2), PBD(4), BSD(3), BCDWIT(2//3, 3), 
    EBCD(2//3), ABCD(2), GBCD(1), GBCD(2), GBCD(5), CRD()
]

# setting labels
label = [set_label(item) for item in rnd] 

# running simulations
@time rct = [simulate(item, nsbj, nsim, seed) for item in rnd];

# calculating operational characteristics
@time op = vcat([calculate_op(item) for item in rct]...)

# saving output
CSV.write("output/case-study01-example1-op.csv", op)

shapes = [
    Shape.circle, Shape.square, Shape.diamond, Shape.cross,
    Shape.utriangle, Shape.dtriangle, Shape.ltriangle, Shape.rtriangle,
    Shape.xcross, Shape.hexagon, Shape.star1, Shape.star2
]

# expected absolute imbalance plot (fig01)
fig01 = Gadfly.plot(
    sort(op, [:design, :sbj]),
    x = :sbj, 
    y = :Imb, 
    color = :design,
    shape = :design,    
    Geom.point,
    Geom.line,
    Scale.color_discrete_manual(
        "#8c0808", "#224966", "#001b58", "#e49b14", 
        "#0d660d", "#000000", "#b1698c", "#76a5af",
        "#0b5394", "#134f5c", "#808080", "#b45f06"
    ),
    Guide.xlabel("Allocation step"), 
    Guide.ylabel("Expected absolute imbalance"),
    Guide.shapekey(title=""),
    Theme(
        point_size = 7pt,
        point_shapes = shapes, 
        minor_label_font = "Serif",
        minor_label_font_size = 12pt,
        major_label_font = "Serif",
        major_label_font_size = 14pt,  
        key_label_font = "Serif",
        key_label_font_size = 12pt,
        key_position = :bottom,
        background_color = "white"
    )
)

# saving absolute imbalance plot:
# pdf format
fig01_pdf = "output/case-staudy01-example1-fig01.pdf"
draw(PDF(fig01_pdf, 12inch, 8inch, dpi=300), fig01)

# png format
fig01_png = "output/case-staudy01-example1-fig01.png"
draw(PNG(fig01_png, 12inch, 8inch, dpi=300), fig01)


# expected proportion of correct guesses (fig02)
fig02 = Gadfly.plot(
    sort(op, [:design, :sbj]),
    x = :sbj, 
    y = :PCG, 
    color = :design,
    shape = :design,    
    Geom.point,
    Geom.line,
    Scale.color_discrete_manual(
        "#8c0808", "#224966", "#001b58", "#e49b14", 
        "#0d660d", "#000000", "#b1698c", "#76a5af",
        "#0b5394", "#134f5c", "#808080", "#b45f06"
    ),
    Guide.xlabel("Allocation step"), 
    Guide.ylabel("Expected proportion of correct guesses"),
    Guide.shapekey(title=""),
    Theme(
        point_size = 7pt,
        point_shapes = shapes, 
        minor_label_font = "Serif",
        minor_label_font_size = 12pt,
        major_label_font = "Serif",
        major_label_font_size = 14pt,  
        key_label_font = "Serif",
        key_label_font_size = 12pt,
        key_position = :bottom,
        background_color = "white"
    )
)

# saving PCG plot:
# pdf format
fig02_pdf = "output/case-staudy01-example1-fig02.pdf"
draw(PDF(fig02_pdf, 12inch, 8inch, dpi=300), fig02)

# png format
fig02_png = "output/case-staudy01-example1-fig02.png"
draw(PNG(fig02_png, 12inch, 8inch, dpi=300), fig02)


# aggregated expected loss (Imb(n)) vs forcing index (FI(n)),
# where n is a sample size (fig03)
op_ = @pipe op |>
    groupby(_, :design) |>
    combine(_,
        :sbj,
        :FI, 
        :Loss => (x -> cumsum(x) ./ (1:length(x))) => :aLoss
    )

fig03 = Gadfly.plot(
    filter(:sbj => (x -> x == nsbj), sort(op_, [:design, :sbj])),
    x = :FI, 
    y = :aLoss, 
    color = :design,
    shape = :design,    
    Geom.point,
    Geom.line,
    Scale.color_discrete_manual(
        "#8c0808", "#224966", "#001b58", "#e49b14", 
        "#0d660d", "#000000", "#b1698c", "#76a5af",
        "#0b5394", "#134f5c", "#808080", "#b45f06"
    ),
    Guide.xticks(ticks = 0:0.25:1),
    Guide.yticks(ticks = 0:0.25:1),
    Guide.xlabel("Allocation step"), 
    Guide.ylabel("Aggregated expected loss, Imb(n)"),
    Guide.shapekey(title=""),
    Theme(
        point_size = 7pt,
        point_shapes = shapes, 
        minor_label_font = "CMU Serif",
        minor_label_font_size = 12pt,
        major_label_font = "CMU Serif",
        major_label_font_size = 14pt,  
        key_label_font = "CMU Serif",
        key_label_font_size = 12pt,
        key_position = :right,
        background_color = "white"
    )
)

# saving PCG plot:
# pdf format
fig03_pdf = "output/case-staudy01-example1-fig03.pdf"
draw(PDF(fig03_pdf, 12inch, 8inch, dpi=300), fig03)

# png format
fig03_png = "output/case-staudy01-example1-fig03.png"
draw(PNG(fig03_png, 12inch, 8inch, dpi=300), fig03)


# heatmap plot of balance/randomness tradeoff (fig04)
transform!(op_, [:aLoss, :FI] => ((x, y) -> x.^2 + y.^2) => :d)
op_ = @pipe op_ |> 
    groupby(_, :design) |> 
    combine(_, :d => (x -> x[nsbj]) => :dmax, :sbj, :d)

fig04 = Gadfly.plot(
    sort(op_, [:dmax, :sbj], rev = [true, false]),
    x = :sbj, 
    y = :design, 
    color = :d,
    Geom.rectbin,
    Scale.ContinuousColorScale(
        palette -> get(ColorSchemes.inferno, palette)
    ),
    Coord.cartesian(xmin=1, xmax=nsbj),# x_continuous(maxvalue = 50),
    Guide.xlabel("Allocation step"), 
    Guide.ylabel(""),
    Guide.colorkey(title=""),
    Theme(
        point_size = 7pt,
        point_shapes = shapes, 
        minor_label_font = "CMU Serif",
        minor_label_font_size = 12pt,
        major_label_font = "CMU Serif",
        major_label_font_size = 14pt,  
        key_label_font = "CMU Serif",
        key_label_font_size = 12pt,
        key_position = :right,
        background_color = "white"
    )
)

# saving PCG plot:
# pdf format
fig04_pdf = "output/case-staudy01-example1-fig04.pdf"
draw(PDF(fig04_pdf, 12inch, 8inch, dpi=300), fig04)

# png format
fig04_png = "output/case-staudy01-example1-fig04.png"
draw(PNG(fig04_png, 12inch, 8inch, dpi=300), fig04)


# inference for 4 different response models
function inference(μ1::Float64, μ0::Float64, σ::Float64, rct::SimulatedRCT, seed::Int64)
    nsbj, _ = size(rct.trt)

    R1m1 = Normal(μ1, σ)
    R0m1 = Normal(μ0, σ)
    m1rsp = generate_rsp(R1m1, R0m1, rct.trt, seed)
    rctm1 = set_response(rct, m1rsp)

    m2rsp = add_tt(m1rsp, 5 .*collect(1:nsbj) ./ (nsbj + 1))
    rctm2 = set_response(rct, m2rsp)
    
    R1m3 = Cauchy(μ1)
    R0m3 = Cauchy(μ0)
    m3rsp = generate_rsp(R1m3, R0m1, rct.trt, seed)
    rctm3 = set_response(rct, m3rsp)

    m4rsp = add_sb(rct.trt, m1rsp, 0.5)
    rctm4 = set_response(rct, m4rsp)

    model = ["Errors N(0, 1)", "Linear trend", "Errors Caushy", "Selection bias"]
    error1 = [t_test(rctm, 0.05) for rctm in [rctm1, rctm2, rctm3, rctm4]]
    error2 = [randomization_test(rctm, 0.05, false) for rctm in [rctm1, rctm2, rctm3, rctm4]]
    error3 = [randomization_test(rctm, 0.05, true) for rctm in [rctm1, rctm2, rctm3, rctm4]]
    
    return DataFrame(
        design = rct.label, 
        model = model, 
        t_test = error1, 
        rando_test = error2,
        rando_test_ranks = error3
    )
end

μ1_null = 0.00 # Null            (Type I error ~ 5%)  
μ1_alt1 = 0.73 # Alternative I   (power ~ 70%)   
μ1_alt2 = 0.82 # Alternative II  (power ~ 80%)
μ1_alt3 = 0.95 # Alternative III (power ~ 90%)

@time null = vcat([inference(μ1_null, μ0, 1.0, rct[s], 314159) for s in eachindex(rct)]...)
insertcols!(null, :scenario => "Null", :xintercept => 0.05)

@time alt1 = vcat([inference(μ1_alt1, μ0, 1.0, rct[s], 314159) for s in eachindex(rct)]...)
insertcols!(alt1, :scenario => "Alternative 1", :xintercept => 0.70)

@time alt2 = vcat([inference(μ1_alt2, μ0, 1.0, rct[s], 314159) for s in eachindex(rct)]...)
insertcols!(alt2, :scenario => "Alternative 2", :xintercept => 0.80)

@time alt3 = vcat([inference(μ1_alt3, μ0, 1.0, rct[s], 314159) for s in eachindex(rct)]...)
insertcols!(alt3, :scenario => "Alternative 3", :xintercept => 0.90)

inference_df = @pipe vcat([null, alt1, alt2, alt3]...) |>
    stack(_, [:t_test, :rando_test, :rando_test_ranks]) |>
    transform(_, 
        :value => (x -> x*100) => :value,
        :xintercept => (x -> x*100) => :xintercept,
    )

CSV.write("output/case-study01-example1-inference.csv", inference_df)

fig05 = Gadfly.plot(
    inference_df, 
    x = :value, 
    y = :design, 
    xgroup = :scenario,
    ygroup = :model,
    shape = :variable,
    xintercept = :xintercept,
    Geom.subplot_grid(
        Geom.point, 
        Geom.vline(color="black", size = 2pt, style = :dash),
        Guide.xticks(ticks = [[1, 5]; 20:10:90])
    ),
    Guide.xlabel("Type I error/power (%)"), 
    Guide.ylabel(""),
    Guide.shapekey(title = ""),
    Theme(
        default_color = "#004d7f",
        point_size = 4pt,
        point_shapes = shapes,
        minor_label_font = "CMU Serif",
        minor_label_font_size = 10pt,
        major_label_font = "CMU Serif",
        major_label_font_size = 12pt,  
        key_label_font = "CMU Serif",
        key_label_font_size = 12pt,
        key_position = :bottom,
        background_color = "white"
    )
)

# saving PCG plot:
# pdf format
fig05_pdf = "output/case-staudy01-example1-fig05.pdf"
draw(PDF(fig05_pdf, 12inch, 8inch, dpi=300), fig05)

# png format
fig05_png = "output/case-staudy01-example1-fig05.png"
draw(PNG(fig05_png, 12inch, 8inch, dpi=300), fig05)

