


function CheckValidDirectory(dir::String)

    # check that this is in fact a directory path ("/" at the end, "" return true)
    if !isdirpath(dir)
        println("SecularResponse.IO.CheckValidDirectory: nonexistent directory $dir.")
        return false
    end

    # finally, try opening a file to check write permissions
    try
        tst = open(dir*"tst.dat", "w")
        close(tst)
        rm(dir*"tst.dat")
    catch e
        println("SecularResponse.IO.CheckValidDirectory: cannot write test file to $dir")
        return false
    end

    return true
end


"""
    SecularFilename()

"""
function SecularFilename(params::SecularParameters,couplingname::String)

    Linearparams = params.Linearparams
    return params.secdir*couplingname*"_"*Linearparams.modelname*"_df_"*Linearparams.dfname*"_l_"*string(params.lharmonic)*"_n1_"*string(params.n1max)*"_rb_"*string(Linearparams.rbasis)*"_Kv_"*string(params.Kv)*".h5"
end

"""
write all the parameters to a file
"""
function WriteParameters(filename::String,
                         params::SecularParameters,
                         mode::String="r+")

    h5open(filename, mode) do file
        WriteParameters(file,params)
    end
end

function WriteParameters(file::HDF5.File,
                         params::SecularParameters)

    group = create_group(file,"SecularParameters")
    for i = 1:fieldcount(SecularParameters)
        varname = string(fieldname(SecularParameters,i))
        if (varname == "tabResPair")
            continue
        elseif (varname == "Linearparams")
            LR.WriteParameters(file,params.Linearparams)
        else
            try write(group,varname,getfield(params,i)) catch; println("Unable to write parameter: "*varname) end
        end
    end
end