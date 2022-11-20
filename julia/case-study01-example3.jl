# Example 3 of the simulation study done for the publication:
# -   title: A roadmap to using randomization in clinical trials
# - authors: Vance W. Berger, Louis Joseph Bour, Kerstine Carter, 
#            Jonathan J. Chipman, Colin C. Everett, Nicole Heussen, 
#            Catherine Hewitt, Ralf-Dieter Hilgers, Yuqun Abigail Luo, 
#            Jone Renteria, Yevgen Ryeznik, Oleksandr Sverdlov & Diane Uschner
# - journal: BMC Medical Research Methodology 
# -  volume: 21
# -  number: 168
# -    year: 2021

using CSV
using DataFrames 
using Distributions  
using FLoops  
using Gadfly
using LinearAlgebra
using Random: seed! 
using SharedArrays
using Statistics 
using StatsModels
using Survival
import Cairo, Fontconfig

include("types.jl")
include("set-label.jl")
include("probabilities.jl")
include("simulated-rct.jl")
include("simulation.jl")

# baseline log serum bilirubin level data from the azathioprine trial, 
# reproduced from Fig. 1 of Altman DG, Royston JP (1988) 
# "The hidden effect of time". SiM;7:629–37.
tt_file = CSV.File("./data/case-study01-example3.csv"; header = 1, delim = ",");
tt_data = DataFrame(tt_file)
u = tt_data.Si

fig06 = Gadfly.plot(
    tt_data,
    x = :sbj, 
    y = :Si, 
    Geom.point,
    Guide.xlabel("Subject"), 
    Guide.ylabel("Cusum"),
    Theme(
        default_color = "#004d7f",
        point_size = 4pt,
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

# saving plot:
# pdf format
fig06_pdf = "output/case-study01-example3-fig06.pdf"
draw(PDF(fig06_pdf, 12inch, 8inch, dpi=300), fig06)

# png format
fig06_png = "output/case-study01-example3-fig06.png"
draw(PNG(fig06_png, 12inch, 8inch, dpi=300), fig06)

# function to generate rasponse from the Cox PH model
function response(trt::Array{Int64}, u::Vector{Float64}, logHR::Float64)
    nsbj, nsim = size(trt)
    rsp = zeros(Float64, nsbj, nsim)

    for s in 1:nsim
        δ = trt[:, s]
        λ = exp.(δ.*logHR .+ 0.1*u)
        rsp[:, s] = rand.(Exponential.(1 ./ λ))
    end

    return rsp
end

# function to test HO: HR = 1 vs. H1: HR ≠ 1 using population test
function coxph_pop_test(R::Vector{Float64}, δ::Vector{Int64}, u::Vector{Float64}, adjust::Bool)
    obs = DataFrame(T = R, event = 1, delta = δ, u = u)
    obs.event = EventTime.(obs.T, obs.event .== 1);
            
    if adjust
        model_fit = coxph(@formula(event ~ delta + u), obs)
    else
        model_fit = coxph(@formula(event ~ delta), obs)
    end

    est = coef(model_fit)[1]
    se = stderror(model_fit)[1]
    zstat = est/se
    
    return 2*ccdf(Normal(), abs(zstat))
end


# simulation of inference using two approaches
nsbj = length(u)
nsim = 10000
seed = 314159

rnd = [CRD(), RAND(nsbj), TBD(nsbj), PBD(2), PBD(4), BSD(3), GBCD(2)]
rct = [simulate(rnd[i], nsbj, nsim, seed) for i in eachindex(rnd)];
trt = [rct[i].trt for i in eachindex(rnd)];

seed!(seed);
# generating responses under H0: HR = 1
rsp_h0 = [response(trt[i], u, log(1)) for i in eachindex(trt)];

# generating responses under H1: HR = 0.6
rsp_h1 = [response(trt[i], u, log(0.6)) for i in eachindex(trt)];

# testing H0 when H0 is true
# adjusted for u
type_I_error_adj = zeros(Float64, length(rnd));

# non-adjusted for u
type_I_error_nadj = zeros(Float64, length(rnd));

# significance level α
α = 0.05
@time for i in eachindex(rnd)
    # adjusted for u
    pvalue = vcat([coxph_pop_test(rsp_h0[i][:, s], trt[i][:, s], u, true) for s in 1:nsim]...);
    type_I_error_adj[i] = mean(pvalue .< α)
    
    # non-adjusted for u
    pvalue = vcat([coxph_pop_test(rsp_h0[i][:, s], trt[i][:, s], u, false) for s in 1:nsim]...);
    type_I_error_nadj[i] = mean(pvalue .< α)
end

# testing H0 when H1 is true
# adjusted for u
power_adj = zeros(Float64, length(rnd));

# non-adjusted for u
power_nadj = zeros(Float64, length(rnd));

@time for i in eachindex(rnd)
    # adjusted for u
    pvalue = vcat([coxph_pop_test(rsp_h1[i][:, s], trt[i][:, s], u, true) for s in 1:nsim]...);
    power_adj[i] = mean(pvalue .< α)
    
    # non-adjusted for u
    pvalue = vcat([coxph_pop_test(rsp_h1[i][:, s], trt[i][:, s], u, false) for s in 1:nsim]...);
    power_nadj[i] = mean(pvalue .< α)
end

population_test = DataFrame(
    design = [rct[i].label for i in eachindex(rct)],
    type_I_error_nadj = type_I_error_nadj,
    type_I_error_adj = type_I_error_adj,
    power_nadj = power_nadj,
    power_adj = power_adj
)

# Below, we are applying randomization based test.
# we want to parallelize the process 
using Distributed
addprocs(4)
nprocs()

@everywhere begin
    using SharedArrays
    using Survival
    using StatsModels
    using DataFrames
end

@eval @everywhere rsp_h0 = $(rsp_h0)
@eval @everywhere rsp_h1 = $(rsp_h1)
@eval @everywhere trt = $(trt)
@eval @everywhere u = $(u)

@everywhere function fit_coxph_adj(obs::DataFrame)
    model_fit = coxph(@formula(event ~ delta + u), obs)

    est = coef(model_fit)[1]
    se = stderror(model_fit)[1]
    
    return est/se
end

@everywhere function fit_coxph_nadj(obs::DataFrame)
    model_fit = coxph(@formula(event ~ delta), obs)

    est = coef(model_fit)[1]
    se = stderror(model_fit)[1]
    
    return est/se
end

function coxph_rando_test(nseq::Int64, rsp::Array{Float64}, trt::Array{Int64}, u::Vector{Float64}, fit_fcn::Function)
    @eval @everywhere drsp = $(rsp[:, 1:nseq])
    @eval @everywhere dtrt = $(trt[:, 1:nseq])
    @eval @everywhere du = $(u)

    zstat = SharedArray{Float64}(nseq, nseq)
    @sync @distributed for s in axes(z, 2)
        R = drsp[:, s]
        @sync @distributed for j in 1:nseq
            δ = dtrt[:, j]
            obs = DataFrame(T = R, event = 1, delta = δ, u = u)
            obs.event = EventTime.(obs.T, obs.event .== 1);

            zstat[j, s] = fit_fcn(obs) 
        end
    end

    # here, we shift zstat colums such that the first row of the zstat
    # matrix will contain zstat[s, s] (s=1:nseq) after shifting
    zstat_ = hcat([circshift(zstat[:, s], 1-s) for s in axes(zstat, 2)]...);
    
    # calculating rando-based p-values
    reject = hcat([abs.(zstat_[2:end, s]) .>= abs.(zstat_[1, s]) for s in 1:nseq]...);
    pvalue = mean(reject, dims = 1);
    
    return pvalue
end

n = 1000 # number of randomization sequences used in rando-based test

# testing H0 when H0 is true
# adjusted for u
type_I_error_adj = zeros(Float64, length(rnd));

# non-adjusted for u
type_I_error_nadj = zeros(Float64, length(rnd));

# significance level α
α = 0.05
@time for i in eachindex(rnd)
    # adjusted for u
    pvalue = coxph_rando_test(n, rsp_h0[i], trt[i], u, fit_coxph_adj);
    type_I_error_adj[i] = mean(pvalue .< α)
    
    # non-adjusted for u
    pvalue = coxph_rando_test(n, rsp_h0[i], trt[i], u, fit_coxph_nadj);
    type_I_error_nadj[i] = mean(pvalue .< α)
end

# testing H0 when H1 is true
# adjusted for u
power_adj = zeros(Float64, length(rnd));

# non-adjusted for u
power_nadj = zeros(Float64, length(rnd));

@time for i in eachindex(rnd)
    # adjusted for u
    pvalue = coxph_rando_test(n, rsp_h1[i], trt[i], u, fit_coxph_adj);
    power_adj[i] = mean(pvalue .< α)
    
    # non-adjusted for u
    pvalue = coxph_rando_test(n, rsp_h1[i], trt[i], u, fit_coxph_nadj);
    power_nadj[i] = mean(pvalue .< α)
end

randomization_test = DataFrame(
    design = [rct[i].label for i in eachindex(rct)],
    type_I_error_nadj = type_I_error_nadj,
    type_I_error_adj = type_I_error_adj,
    power_nadj = power_nadj,
    power_adj = power_adj
)

# saving inference output
CSV.write("output/case-study01-example3-pop-test.csv", population_test)
CSV.write("output/case-study01-example3-rando-test.csv", randomization_test)