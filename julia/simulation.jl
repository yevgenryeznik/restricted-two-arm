# function to make a treatment assignment
function assign_trt(prb::Float64)
    coin = Binomial(1, prb)
    assignment = rand(coin)

    return assignment
end 


# function to simulate restricted randomization
function simulate(rnd::RestrictedRandomization, nsbj::Int64, nsim::Int64, seed::Int64)    
    seed!(seed)
    if typeof(rnd) == CRD
        trt = rand(Binomial(), nsbj, nsim)
        prb = 0.5 .* ones(Float64, nsbj, nsim)
    else
        trt = zeros(Int64, nsbj, nsim)
        prb = zeros(Float64, nsbj, nsim)
    
        for s ∈ 1:nsim
            # numbers of subjects allocated to treatments
            N1, N2 = 0, 0

            for j ∈ 1:nsbj
                # calculating probability of treatment assignment, 
                # given the current value of imbalance
                prb[j, s] = calculate_prb(rnd, N1, N2)
        
                # here, the treatment assignment is made
                trt[j, s] = assign_trt(prb[j, s])
        
                # updating treatment numbers and the value of imbalance
                N1 += trt[j, s]
                N2 = j - N1
            end
        end
    end

    return SimulatedRCT(set_label(rnd), nsbj, nsim, trt, prb, [])
end
