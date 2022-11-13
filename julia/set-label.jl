# A utility function to create a descriptive label for a given randomization procedures.
function set_label(rnd::RestrictedRandomization)
    label = typeof(rnd).name.wrapper

    if hasfield(typeof(rnd), :param)
        param = rnd.param
        if typeof(param) == Rational{Int64}  # param is a Rational number
            n = numerator(param)             # get numerator
            d = denominator(param)           # get denominator
            label = "$label($n/$d)"
        elseif typeof(param) == Int64        # param is Integer number
            label = "$label($param)"
        else
            param = round(param, digits = 2) # otherwise, round it to 2 disgits
            label = "$label($param)"
        end
    else
        if hasfield(typeof(rnd), :param1) & hasfield(typeof(rnd), :param2) 
            param = (rnd.param1, rnd.param2)
            param_ = ["", ""]
            for i in eachindex(param)
                if typeof(param[i]) == Rational{Int64}  # param is a Rational number
                    n = numerator(param[i])             # get numerator
                    d = denominator(param[i])           # get denominator
                    param_[i] = "$n/$d"
                elseif typeof(param[i]) == Int64        # param is Integer number
                    param_[i] = "$(param[i])"
                else
                    param[i] = round(param[i], digits = 2) # otherwise, round it to 2 disgits
                    param_[i] = "$(param[i])"
                end
            end
            label = "$label($(param_[1]), $(param_[2]))"
        else
            label = "$label"
        end
    end
    
    return label
end