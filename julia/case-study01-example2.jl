# Example 2 of the simulation study done for the publication:
# -   title: A roadmap to using randomization in clinical trials
# - authors: Vance W. Berger, Louis Joseph Bour, Kerstine Carter, Jonathan J. Chipman, Colin C. Everett, Nicole Heussen, 
#            Catherine Hewitt, Ralf-Dieter Hilgers, Yuqun Abigail Luo, Jone Renteria, Yevgen Ryeznik, Oleksandr Sverdlov & 
#            Diane Uschner
# - journal: BMC Medical Research Methodology 
# -  volume: 21
# -  number: 168
# -    year: 2021

using CSV
using DataFrames
using Distributions
using Pipe
using Random: seed!

include("simulated-rct.jl")
include("types.jl")
include("set-label.jl")
include("probabilities.jl")
include("simulation.jl")
include("calculate-op.jl")

# sample size
nsbj = 96

# number of simulations
nsim = 10000

# seed
seed = 3141592

# randomization procedures
rnd = [PBD(2), PBD(4), PBD(6), BSD(1), BSD(2), BSD(3)]

# setting labels
label = [set_label(item) for item in rnd] 

# running simulations
@time rct = [simulate(item, nsbj, nsim, seed) for item in rnd];

# calculating operational characteristics
op = vcat([calculate_op(rct[i]) for i in eachindex(rct)]...)

# utility functions to create summary table of predictability
minus50(x) = [round(100*x[i] - 50, digits = 1) for i in eachindex(x)]
extract_design_name(x) = [replace(x[i], match(r"\W\d\W", x[i]).match => "") for i in eachindex(x)]
function extract_mti(x)
    out = zeros(Int64, length(x))
    design = extract_design_name(x)
    for i in eachindex(x)
        out[i] = parse(Int64, match(r"\d", x[i]).match)
        if design[i] == "PBD"
            out[i] = out[i]/2
        end
    end
    return out
end

# function to calculate proportion of deterministic assignments
function calculate_pda(rct::SimulatedRCT)
    prb = rct.prb
    nsbj = rct.nsbj

    # expected number of deterministic assignments
    nda = mean(sum(prb .== 0, dims = 1) + sum(prb .== 1, dims = 1))

    # proportion of deterministic assignments
    return round(nda/nsbj*100, digits = 1)
end

# summary: Excess Correct Guess Probability
ecgp = @pipe op |> 
    filter(:sbj => (x -> x == nsbj), _) |> 
    select(_, :design, :PCG) |> 
    transform(_, :PCG => minus50  => :ECGP) |> # Excess Correct Guess Probability
    transform(_, 
        :design => extract_mti => :MTI,    
        :design => extract_design_name  => :design
    ) |>
    select(_, Not([:PCG])) |>
    unstack(_, :MTI, :design, :ECGP) |>
    rename(_, :PBD => :ECGP_PBD, :BSD => :ECGP_BSD)

    # summary: Proportion of Deterministic Assignments
pda = @pipe DataFrame(
        design = label,
        PDA = [calculate_pda(rct[i]) for i in eachindex(rct)]
        ) |>
        transform(_, 
            :design => extract_mti => :MTI,    
            :design => extract_design_name  => :design
        ) |>
        unstack(_, :MTI, :design, :PDA) |>
        rename(_, :PBD => :PDA_PBD, :BSD => :PDA_BSD) 

# output 
predictability = innerjoin(pda, ecgp, on = :MTI)        

# saving output
CSV.write("output/case-study01-example2-predictability.csv", predictability)