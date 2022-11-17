# Function to simulate responses, given responses' distributions and treatment assignments
function generate_rsp(R1::Distribution, R0::Distribution, trt::Array{Int64}, seed::Int64)
    nsbj, nsim = size(trt)
    rsp = zeros(nsbj, nsim)

    seed!(seed)
    rsp1 = rand(R1, nsbj, nsim)
    rsp0 = rand(R0, nsbj, nsim)

    for s ∈ 1:nsim
        for j ∈ 1:nsbj
            rsp[j, s] = trt[j, s] == 1 ? rsp1[j, s] : rsp0[j, s]
        end
    end

    return rsp
end


# Function to add a time trend to the responses
function add_tt(rsp::Array{T}, tt::Vector{Float64}) where {T <: Number}
    nsbj, nsim = size(rsp)
    rsp_tt = zeros(nsbj, nsim)

    for s ∈ 1:nsim
        for j ∈ 1:nsbj
            rsp_tt[j, s] = rsp[j, s] + tt[j] 
        end
    end

    return rsp_tt
end


# Function to simulate selection bias
function add_sb(trt::Array{Int64}, rsp::Array{T}, sb::Number) where {T <: Number}
    nsim = size(trt, 2)
    imb = hcat(map(1:nsim) do col
        δ = trt[:, col]
        imb_ = cumsum(2 .*δ .- 1)
        
        return [0; imb_[1:end-1]]
    end...)
    
    return rsp - sign.(imb).*sb
end