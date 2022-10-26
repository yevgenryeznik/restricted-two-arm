function assign_trt(prb::Float64)
    coin = Binomial(1, prb)
    assignment = rand(coin)

    return assignment
end 

function randomize(nsbj::Int64, rnd::Function; kwargs...)
    trt = zeros(Int64, nsbj)
    prb = zeros(Float64, nsbj)
    
    # number of subjects allocated to treatment 1
    N1 = 0

    # number of subjects allocated to treatment 2
    N2 = 0

    for j ∈ 1:nsbj
        # calculating probability of treatment assignment, 
        # given the current value of imbalance
        prb[j] = rnd(N1, N2; kwargs...)

        # here, the treatment assignment is made
        trt[j] = assign_trt(prb[j])

        # updating treatment numbers and the value of imbalance
        N1 += trt[j]
        N2 = j - N1
    end

    return trt, prb
end
