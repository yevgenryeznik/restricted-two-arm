# A user-defined type to represent a simulated Randomized Clinical Trial (RCT)
struct SimulatedRCT{T}
    label::String       # name of RR procedure
    nsbj::Int64         # sample size  
    nsim::Int64         # number of trial simulations 
    trt::Array{Int64} # treatment sequences
    prb::Array{Float64} # probabilities of treatment assignments
    rsp::Array{T}       # subjects' responses 
end


# Ovveriding a `Base` function `show` to print out an instance of `SimulatedRCT``
function Base.show(io::IO, rct::SimulatedRCT)
    nsbj, nsim, label = rct.nsbj, rct.nsim, rct.label

    println("Simulated RCT, targeting 1:1 allocation:")
    println("   - $(label) randomization;")
    println("   - $(nsbj) subjects;")
    println("   - $(nsim) simulations.")
    println("\nTreatment assignments:")
    print("trt = ")
    display(rct.trt)
    println("\nProbabilities of treatment assignments:")
    print("prb = ")
    display(rct.prb)
    println("\nSubjects' responses:")
    print("rsp = ")
    if rct.rsp == []
        println("[ ] (not set)")
    else
        display(rct.rsp)
    end
end