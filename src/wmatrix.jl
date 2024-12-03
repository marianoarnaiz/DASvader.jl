"This is the Whole Matrix operations part of DAS VADER V1.0"


export nan2zero!, inf2zero!, There, squaredas!, squaredas, absdas!, absdas, chop!, chop

"""
nan2zero!: Remove NaN dDAS IN PLACE

Input:
 - dDAS: DAS data structure

Outputs:
- dDAS: All NaN are trasnformed to zero (0.0)

Notes:
 - No figure is created.
 - There is not-in-place version of this function
# Example: Save one channel and its time vector to a file.
```
julia> nan2zero!(dDAS)
```
"""
function nan2zero!(dDAS)
    for i in eachindex(dDAS.data)
        @inbounds dDAS.data[i] = ifelse(isnan(dDAS.data[i]), 0.0, dDAS.data[i])
    end
end
#############################################################################################

"""
inf2zero!: Remove Inf dDAS IN PLACE

Input:
 - dDAS: DAS data structure

Outputs:
- dDAS: All Inf are trasnformed to zero (0.0)

Notes:
 - No figure is created.
 - There is not-in-place version of this function
# Example: Save one channel and its time vector to a file.
```
julia> inf2zero!(dDAS)
```
"""
function inf2zero!(dDAS)
    for i in eachindex(dDAS.data)
        @inbounds dDAS.data[i] = ifelse(isinf(dDAS.data[i]), 0.0, dDAS.data[i])
    end
end


#############################################################################################

"""
inf2zero!: Remove Inf dDAS IN PLACE

Input:
 - dDAS: DAS data structure

Outputs:
- dDAS: All Inf are trasnformed to zero (0.0)

Notes:
 - No figure is created.
 - There is not-in-place version of this function
# Example: Save one channel and its time vector to a file.
```
julia> inf2zero!(dDAS)
```
"""
function squaredas!(dDAS)
    dDAS.data = dDAS.data .^ 2
end
function squaredas(dDAS)
    unDAS = deepcopy(dDAS)
    unDAS.data = dDAS.data .^ 2
    return unDAS
end


function absdas!(dDAS)
    dDAS.data = abs.(dDAS.data)
end


function absdas(dDAS)
    unDAS = deepcopy(dDAS)
    unDAS.data = abs.(dDAS.data)
    return unDAS
end




function chop!(dDAS; type=:negative)
    if type == :negative
        idx = findall(dDAS.data .< 0.0)
        dDAS.data[idx] .= 0.0
    elseif type == :positive
        idx = findall(dDAS.data .> 0.0)
        dDAS.data[idx] .= 0.0
    end
end


function chop(dDAS; type=:negative)
    unDAS = deepcopy(dDAS)

    if type == :negative
        idx = findall(dDAS.data .< 0.0)
        unDAS.data[idx] .= 0.0
    elseif type == :positive
        idx = findall(dDAS.data .> 0.0)
        unDAS.data[idx] .= 0.0
    end
    return unDAS
end
