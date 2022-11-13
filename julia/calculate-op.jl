# function to calculate expected imbalance at each allocation step
function calc_expected_imb(rct::SimulatedRCT)
    # getting simulated treatment sequences
    trt = rct.trt

    # sample size & number of simulations used
    nsbj, nsim = size(trt)

    # imbalance at each allocation step per simulation
    imb = zeros(Int64, nsbj, nsim)

    for s ∈ 1:nsim
        N1 = 0
        for j ∈ 1:nsbj
            N1 += trt[j, s]
            imb[j, s] = 2*N1 - j
        end
    end

    # returning vector of expected imbalance values 
    # at each allocation step
    return vec(mean(imb, dims = 2))
end


# function to calculate expected absolute imbalance at each allocation step
function calc_expected_abs_imb(rct::SimulatedRCT)
    # getting simulated treatment sequences
    trt = rct.trt

    # sample size & number of simulations used
    nsbj, nsim = size(trt)

    # imbalance at each allocation step per simulation
    imb = zeros(Int64, nsbj, nsim)

    for s ∈ 1:nsim
        N1 = 0
        for j ∈ 1:nsbj
            N1 += trt[j, s]
            imb[j, s] = abs(2*N1 - j)
        end
    end

    # returning vector of expected absolute imbalance values 
    # at each allocation step
    return vec(mean(imb, dims = 2))
end


# function to calculate expected loss at each allocation step
function calc_expected_loss(rct::SimulatedRCT)
    # getting simulated treatment sequences
    trt = rct.trt

    # sample size & number of simulations used
    nsbj, nsim = size(trt)

    # loss at each allocation step per simulation
    loss = zeros(Float64, nsbj, nsim)

    for s ∈ 1:nsim
        N1 = 0
        for j ∈ 1:nsbj
            N1 += trt[j, s]
            loss[j, s] = abs(2*N1 - j)^2/j
        end
    end

    # returning vector of expected loss values 
    # at each allocation step
    return vec(mean(loss, dims = 2))
end


# function to calculate forcing index at each allocation step
function calc_fi(rct::SimulatedRCT)
    # getting probabilities of treatment assignments
    prb = rct.prb

    # sample size & number of simulations used
    nsbj, nsim = size(prb)

    # absolute difference between probability of treatment assinments
    # and CRD probability (0.5) at each allocation step per simulation
    prb_diff = zeros(Float64, nsbj, nsim)

    for s ∈ 1:nsim
        for j ∈ 1:nsbj
            prb_diff[j, s] = abs(prb[j, s] - 0.5)
        end
    end

    # expected absolute difference between probability of treatment 
    # assinments and CRD probability (0.5) at each allocation step
    expected_prb_diff = vec(mean(prb_diff, dims = 2))

    # returning vector of FI values at each allocation step
    return 4 .* cumsum(expected_prb_diff) ./(1:nsbj)
end


# function to calculate proportions of correct guesses up to 
# given allocation step
function calc_pcg(rct::SimulatedRCT)
    # getting simulated treatment sequences
    trt = rct.trt

    # sample size & number of simulations used
    nsbj, nsim = size(trt)

    # treatment guess, given treatment imbalance,
    # at each allocation step per simulation
    guess = zeros(Int64, nsbj, nsim)

    binomial_rv = rand(Binomial(1, 0.5), nsbj, nsim)
    for s ∈ 1:nsim
        N1 = 0 # no treatment assignment yet
        guess[1, s] = binomial_rv[1, s]
        for j ∈ 1:(nsbj-1)
            N1 += trt[j, s] 
            imb = 2*N1 - j # imbalance afer j-th allocation step 
            guess[j+1, s] = imb == 0 ?  binomial_rv[j+1, s] : (imb < 0 ? 1 : 0)
        end
    end

    # probabilities of correct guesses at each allocation step
    prb_guess == vec(mean(guess .== trt, dims = 2))

    # returning proportion of correct guesses up to given allocation step
    return cumsum(prb_guess) ./ (1:nsbj)
end


# function to calculate operational characteristics, given simulation output
function calculate_op(rct::SimulatedRCT)
    design = rct.label                # randomization procedure used
    sbj = collect(1:rct.nsbj)         # subjects' IDs (allocation step)
    Imb = calc_expected_abs_imb(rct)  # expected (average) absolute imbalance
    Loss = calc_expected_loss(rct)    # expected (average) loss
    FI = calc_fi(rct)                 # expected (average) forcing index
    PCG = calc_pcg(rct)               # expected (average) proportions of correct guesses

    # an output as a DataFrame
    rct_summary = DataFrame(design = design, sbj = sbj, Imb = Imb, Loss = Loss, FI = FI, PCG = PCG)

    return rct_summary
end