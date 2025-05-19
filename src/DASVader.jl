"DAS VADER MODULE"

module DASVader

# in your main module or script
if get(ENV, "DISPLAY", "") != ""
    include("vizGL.jl")
    println("Using GLMakie (Display detected)")
else
    include("vizCairo.jl")
    println("Using CairoMakie (Headless)")
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

include("headless2.jl")
#using .HeadlessCheck
#
if has_graphical_display()
    include("vizGL.jl")
    #using .vizGL
    @info "Working on graphical mode"
else
    include("vizCairo.jl")
    #using .vizCairo
    @info "Working on headless mode"
end

include("tandf.jl")
include("wmatrix.jl")
include("filemanagement.jl")
include("myspectra.jl")
include("xandw.jl")
include("detectiontools.jl")
include("Findpeaks.jl")
include("myspectogram.jl")
include("datasrc.jl")
include("interface.jl")
