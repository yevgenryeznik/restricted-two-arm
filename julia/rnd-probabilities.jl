# Here, we define a set of functions used to calculate probabilities of 
# treatment assignments, given restricted radomization procedure.
# Input parameters indicate:
#   - rnd -- a type of rundomization procedure used.
#   -  N1 -- a number of subjects allocated to treatment 1.
#   -  N2 -- a number of subjects allocated to treatment 2.


# Completely Randomized Design
function calculate_prb(rnd::CRD, N1::Int64, N2::Int64)
    # getting probability of a coin: 
    # it is a probability of treatment assignment
    p = rnd.param

    return p
end


# Truncated Binomial Design
function calculate_prb(rnd::TBD, N1::Int64, N2::Int64)
    # getting probability of a coin
    p = rnd.param1

    # getting total sample size
    nsbj = rnd.param2

    # calculating probability of treatment assignment
    prb = max(N1, N2) < nsbj/2 ? p : (N1 < N2 ? 1 : 0)

    return prb
end


# Permuted Block Design
function calculate_prb(rnd::PBD, N1::Int64, N2::Int64)
    # getting a block size
    bs = rnd.param

    # indicating a current allocation step (current subject's ID)
    j = N1 + N2 + 1

    # indicating a current subject's position in a block:
    # 1 <= k <= bs
    k = (j % bs == 0)*bs + (j % bs)

    # calculating the number of subjects in a block 
    # allocated to the 1st treatment arm (n1)
    n1 = k == 1 ? 0 : N1 - ((j - 1) ÷ bs)*(bs/2) 

    # calculating probability of treatment assignment
    prb = (0.5*bs-n1)/(bs-k+1)

    return prb
end


# Random Allocation Rule (Rand)
function calculate_prb(rnd::RAND, N1::Int64, N2::Int64)
    # getting a block size (for Rand, it is equal to the total sample size)
    nsbj = rnd.param

    # indicating a current allocation step (current subject's ID) --
    # it is equal to a current subject's position in a block.
    # 1 <= j <= nsbj (nsbj == block size)
    j = N1 + N2 + 1

    # calculating probability of treatment assignment
    prb = (0.5*nsbj-N1)/(nsbj-j+1)

    return prb
end


# Efron's Biased Coin Design
function calculate_prb(rnd::EBCD, N1::Int64, N2::Int64)
    # getting a probability of a biased coin
    p = rnd.param

    # calculating current treatment imbalance
    d = N1 - N2
 
    # calculating probability of tretament asssgnment, given imbalance
    prb = abs(d) == 0 ? 0.5 : (d < 0 ? p : 1-p)    

    return prb
end


# Adjustable Biased Coin Design
function calculate_prb(rnd::ABCD, N1::Int64, N2::Int64)
    # probability function proposed by A. Baldi Antognini and A. Giovagnoli
    F(x, a) = abs(x)^a/(abs(x)^a + 1)

    # getting parameter of the randomization procedure
    a = rnd.param

    # calculating current treatment imbalance
    d = N1 - N2

    # calculating probability of tretament asssgnment, given imbalance
    prb = abs(d) <= 1 ? 0.5 : (d < -1 ? F(d, a) : 1-F(d, a))    

    return prb
end


# Generalized Biased Coin Design
function calculate_prb(rnd::GBCD, N1::Int64, N2::Int64)
    # getting a parameter value of the randomization procedure
    γ = rnd.param

    # indicating a current allocation step (current subject's ID)
    j = N1 + N2 + 1

    # calculating probability of tretament asssgnment
    prb = j == 1 ? 0.5 : N2^γ/(N1^γ + N2^γ) 

    return prb
end


# Big Stick Design
function calculate_prb(rnd::BSD, N1::Int64, N2::Int64)
    # getting MTI parameter
    mti = rnd.param

    # calculating current treatment imbalance
    d = N1 - N2
 
    # calculating probability of tretament asssgnment, given imbalance
    prb = abs(d) < mti ? 0.5 : (d == mti ? 0 : 1)    

    return prb
end


# Biased Coin Design With Imbalance Tolerance
function calculate_prb(rnd::BCDWIT, N1::Int64, N2::Int64)
    # getting a probability of a biased coin
    p = rnd.param1

    # getting MTI parameter
    mti = rnd.param2

    # calculating current treatment imbalance
    d = N1 - N2
 
    # calculating probability of tretament asssgnment, given imbalance
    prb = abs(d) < mti ? (d == 0 ? 0.5 : (d < 0 ? p : 1-p)) : (d == mti ? 0 : 1)  

    return prb
end






