# Here, we define different types of restricted
# randomization procedures


abstract type RestrictedRandomization end


# Completely Randomized Design (like tossing a coin)
struct CRD <: RestrictedRandomization 
end


# Truncated Binomial Design
# subjects are allocated acoring to CRD (tossing a coin ) until one of 
# treatments gets enough number of subjects (half of the total sample size). 
# Then, the rest of the subjects is deterministically allocated to another 
# tretament. 
#   - `nsbj` represents a total sample size in a study. 
struct TBD <: RestrictedRandomization 
    nsbj::Int64
end


# Permuted Block Design
# randomization is conducted in blocks such that at the end of each block, 
# there is an exact balance in treatment assignments .
# param reperesents a block size
struct PBD <: RestrictedRandomization
    param::Int64
end


# Random Allocation Rule (Rand)
# A type of PBD, where block size equals to the total sample size.
#   - `nsbj` represents a total sample size in a study. 
struct RAND <: RestrictedRandomization
    nsbj::Int64
end


# Efron's Biased Coin Design
#   - `param` represents a probability of a biased coin (p).
#     p = 2/3 is recommended.
struct EBCD <: RestrictedRandomization
    param::Number
end


# Adjustable Biased Coin Design
#   - `param`represents a parameter (a) of a probability function used 
#     in treatment assignments. 
struct ABCD <: RestrictedRandomization
    param::Number
end


# Generalized Biased Coin Design
#   - `param` represents a parameter (γ) of a probability function used 
#     in treatment assignments. 
struct GBCD <: RestrictedRandomization
    param::Number
end


# Big Stick Design
#   - `param`` represents maximum tolerated imbalance (MTI), a parameter that 
#     controls treatment imbalance. Thus, if at some allocation step, an absolute 
#     treatment imbalance equals to MTI, then next allocation is deterministically 
#     forced toward balance.
struct BSD <: RestrictedRandomization
    param::Int64
end


# Biased Coin Design With Imbalance Tolerance
# A combination of Efron's Biased Coin Design and MTI procedure.
#   - `param1` represents a probability of BCD.
#   - `param2` represents a MTI parameter.
struct BCDWIT <: RestrictedRandomization
    param1::Number
    param2::Int64
end