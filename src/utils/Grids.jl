

"""
function to make a table of (a, e) values
"""
function gridomega(Amin::Float64,Amax::Float64,nA::Int64,
                   Emin::Float64,Emax::Float64,nE::Int64)

    # compute the total number of (complex) frequencies to probe
    nAE = nA*nE

    # initialise the blank table
    tabAE = zeros(Float64,2,nAE)

    # Range in semi-major axis
    Arange = range(Amin,Amax,length=nA)
    # Range in eccentricity
    Erange = range(Emin,Emax,length=nE)

    # initialise a counter
    icount = 1
    for e in Erange
        for a in Arange
            tabAE[1,i] = a
            tabAE[2,i] = e
            # update the counter
            icount += 1
        end
    end

    return tabAE

end