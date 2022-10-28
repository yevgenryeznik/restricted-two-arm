# Here, we define a set of utility functions to create descriptive labels
# for randomization procedures used.


# Completely Randomized Design
function set_label(rnd::CRD)
    label = "CRD"

    return label
end


# Truncated Binomial Design
function set_label(rnd::TBD)
    label = "TBD"

    return label
end


# Permuted Block Design
function set_label(rnd::PBD)
    param = rnd.param
    label = "PBD($param)"

    return label
end


# Random Allocation Rule (Rand)
function set_label(rnd::RAND)
    label = "Rand"
    
    return label
end


# Efron's Biased Coin Design
function set_label(rnd::EBCD)
    param = rnd.param
    if typeof(param) == Rational{Int64} # param is a Rational number
        n = numerator(param)            # get numerator
        d = denominator(param)          # get denominator
        label = "BCD($n/$d)"
    else
        p = round(param, digits = 2)    # otherwise, round it to 2 disgits
        label = "BCD($p)"
    end

    return label
end 


# Adjustable Biased Coin Design
function set_label(rnd::ABCD)
    param = rnd.param
    label = "ABCD($param)"

    return label
end


# Generalized Biased Coin Design
function set_label(rnd::GBCD)
    param = rnd.param
    label = "GBCD($param)"

    return label
end


# Big Stick Design
function set_label(rnd::BSD)
    param = rnd.param
    label = "BSD($param)"

    return label
end


# Biased Coin Design With Imbalance Tolerance
function set_label(rnd::BCDWIT)
    param1 = rnd.param1
    param2 = rnd.param2

    if typeof(param1) == Rational{Int64} # param1 is a Rational number
        n = numerator(param1)            # get numerator
        d = denominator(param1)          # get denominator
        label = "BCDWIT($n/$d, $param2)"
    else
        p = round(param1, digits = 2)    # otherwise, round it to 2 disgits
        label = "BCDWIT($p, $param2)"
    end

    return label
end