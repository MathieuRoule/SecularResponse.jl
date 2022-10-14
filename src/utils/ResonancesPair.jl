"""
helper functions for making resonance pair containers

"""


function GetNbResPair(lmax::Int64,n1max::Int64,ndim::Int64=3)
    if ndim == 2
        return GetNbResPair2d(lmax,n1max)
    elseif ndim == 3
        return GetNbResPair3d(lmax,n1max)
    else
        error("Unknow dimension in GetNbResPair")
    end
end

function MakeTabResPair(lmax::Int64,n1max::Int64,ndim::Int64=3)
    if ndim == 2
        return MakeTabResPair2d(lmax,n1max)
    elseif ndim == 3
        return MakeTabResPair3d(lmax,n1max)
    else
        error("Unknow dimension in MakeTabResPair")
    end
end


"""
# Function that returns the total number
# of resonances pair ((n1,n2),(n1p,n2p)) to consider
# for the harmonics lharmonic
# There a few constraints to satisfy:
# + |n2| <= lharmonic
# + (lharmonic-n2) even
# + |n1| <= n1max
# + (n1,n2) = (0,0) does not contribute
"""
function GetNbResPair3d(lmax::Int64,n1max::Int64)
    count = 0 # Initialisation of the counter
    #####
    # TO DO !!
    #####
    return count # Returning the number of resonance pairs
end


"""
Container of the ((n1,n2),(n1p,n2p)) resonance pairs to consider
"""
function MakeTabResPair3d(lmax::Int64,n1max::Int64)

    # calculate the number
    nbResPair = GetNbResPair(lmax,n1max,3)
    # Create container array
    tabResPair = zeros(Int64,4,nbResPair)
    count = 1 # Initialisation of the counter
    #####
    # TO DO !!
    #####
    return nbResPair, tabResPair
end


"""
# Function that returns the total number
# of resonances pair ((n1,n2),(n1p,n2p)) to consider
# for lharmonic for discs
# There a few constraints to satisfy:
# + n2 = lharmonic
# + |n1| <= n1max
# + (n1,n2) = (0,0) does not contribute
"""
function GetNbResPair2d(lharmonic::Int64,n1max::Int64)
    count = 0 # Initialisation of the counter
    #####
    if lharmonic == 0 # n1 == 0 not contributing
        for n1 = 1:n1max
            for n1p = 1:n1max
                count += 4
            end
        end
    else
        for n1 = -n1max:n1max
            for n1p = -n1max:n1max
                count += 1
            end
        end
    end
    return count # Returning the number of resonance pairs
end


"""
Container of the ((k1,k2),(k1p,k2p)) resonance pair to consider

@IMPROVE it would be best to use the same code as in get_nbResVec()
"""
function MakeTabResPair2d(lharmonic::Int64,k1max::Int64)
    # calculate the number
    nbResPair = GetNbResPair(lharmonic,k1max,2)
    # Create container array
    tabResPair = zeros(Int64,4,nbResPair)
    count = 1 # Initialisation of the counter
    #####
    if lharmonic == 0 # n1 == 0 not contributing
        for k1 = 1:k1max
            for k1p = 1:k1max
                tabResPair[1,count], tabResPair[2,count], tabResPair[3,count], tabResPair[4,count] = k1, lharmonic, k1p, lharmonic
                count += 1
                tabResPair[1,count], tabResPair[2,count], tabResPair[3,count], tabResPair[4,count] = -k1, lharmonic, k1p, lharmonic
                count += 1
                tabResPair[1,count], tabResPair[2,count], tabResPair[3,count], tabResPair[4,count] = k1, lharmonic, -k1p, lharmonic
                count += 1
                tabResPair[1,count], tabResPair[2,count], tabResPair[3,count], tabResPair[4,count] = -k1, lharmonic, -k1p, lharmonic
                count += 1
            end
        end
    else
        for k1 = -k1max:k1max
            for k1p = -k1max:k1max
                # lharmonic
                tabResPair[1,count], tabResPair[2,count], tabResPair[3,count], tabResPair[4,count] = k1, lharmonic, k1p, lharmonic
                count += 1
            end
        end
    end
    return nbResPair, tabResPair
end
