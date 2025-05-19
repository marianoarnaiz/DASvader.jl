"DAS VADER MODULE"

module DASVader

const VERBOSE = true  # or false, depending on your use case

# in your main module or script
if get(ENV, "DISPLAY", "") != ""
    include("vizGL.jl")
    VERBOSE && println("Using GLMakie (Display detected)")
else
    include("vizCairo.jl")
    VERBOSE && println("Using CairoMakie (Headless)")
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
