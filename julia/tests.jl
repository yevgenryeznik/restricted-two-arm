# function to perform two-sample t-test, given simulated RCT and significance α
function t_test(rct::SimulatedRCT, α::Float64)
    trt = rct.trt
    rsp = rct.rsp

    _, nsim = size(trt)

    reject = map(1:nsim) do s
        δ, y = trt[:, s], rsp[:, s]
        Y1 = y[δ .== 1]
        Y0 = y[δ .== 0]
        test = EqualVarianceTTest(Y1, Y0)

        pvalue(test) < α
    end

    return mean(reject)
end


function randomization_test(rct::SimulatedRCT, α::Float64, use_ranks::Bool)
    trt = rct.trt
    rsp = rct.rsp

    # getting the number of treatment sequences
    nseq = size(rsp, 2)
    
    # initializing a vector of randomization-based p-values
    pvalue = zeros(Float64, nseq)

    # if the linear rank test-based statistics is used
    if use_ranks
        # here, we calculate the linear rank test statistics in the form:
        # 1. rsp -> R (each responses is associated with a correpsonding rank)
        # 2. a = R - mean(R) 
        # 3. S = a' * trt
        a = hcat(map(1:nseq) do col
            Y = rsp[:, col]
            R = ordinalrank(Y)
            return R .- mean(R)
        end...)

        # in the following two lines of code, the multiplication of
        # two matrices is performed: a' and trt.
        # for some reason (!), multiplication by the unit matrix D improves
        # the performance drammatically
        D = Diagonal(ones(nseq))
        S = a' * (trt * D)

    # if the test statistics is not based on ranks
    else
        # here, we calculate the mean response in each treatment group 
        # (YE and YC), and the test statistics is just a difference between 
        # tretament effects: S = YE - YC
        DE = Diagonal(1 ./ sum(trt, dims = 1)[1,:])
        YE = rsp' * (trt * DE)
        
        DC = Diagonal(1 ./ sum(1 .- trt, dims = 1)[1, :])
        YC = rsp' * ((1 .- trt) * DC)
        
        S = YE - YC
    end

    # the randomization-based p-value 
    idx = 1:nseq
    @simd for i ∈ idx
        pvalue[i] = mean(abs.(S[i, i .!= idx]) .>= abs.(S[i, i]))
    end

    return mean(pvalue .< α)
end