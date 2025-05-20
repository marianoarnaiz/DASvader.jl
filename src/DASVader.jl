"DAS VADER MODULE"

module DASVader


# choose the plot backend
#const VERBOSE = true  # or false, depending on your use case
#if get(ENV, "DISPLAY", "") != ""
#    include("vizGL.jl")
#    VERBOSE && println("Using GLMakie (Display detected)")
#   else
#    include("vizCairo.jl")
#    VERBOSE && println("Using CairoMakie (Headless)")
#end


using Requires

function __init__()
    if get(ENV, "DISPLAY", "") != ""
        @require GLMakie = "eaaa29a5-0c1c-538d-82f0-3b444fb5c4ed" begin
            include("vizGL.jl")
            @info "Using GLMakie (Display detected)"
        end
    else
        @require CairoMakie = "13f3f980-e62e-5c1d-bb38-64c898e36e9b" begin
            include("vizCairo.jl")
            @info "Using CairoMakie (Headless)"
        end
    end
end


include("tandf.jl")
include("wmatrix.jl")
include("filemanagement.jl")
#include("vizGL.jl")
include("xandw.jl")
include("detectiontools.jl")
include("Findpeaks.jl")
include("myspectogram.jl")
include("myspectra.jl")
include("datasrc.jl")
include("interface.jl")

end

#############################################
# Life line call!

#include("headless2.jl")
#using .HeadlessCheck
#
#if get(ENV, "DISPLAY", "") != ""
#    include("vizGL.jl")
#    println("Using GLMakie (Display detected)")
#else
#    include("vizCairo.jl")
#    println("Using CairoMakie (Headless)")
#end

#include("tandf.jl")
#include("wmatrix.jl")
#include("filemanagement.jl")
#include("myspectra.jl")
#include("xandw.jl")
#include("detectiontools.jl")
#include("Findpeaks.jl")
#include("myspectogram.jl")
#include("datasrc.jl")
#include("interface.jl")
