"This is the Visualization  part of DAS VADER V1.0"

using MakieThemes, GLMakie, CairoMakie
Makie.set_theme!(ggthemr(:fresh)) #Activate the cool theme :)
GLMakie.activate!() #activate GL plottinge


## Quicker DAS viewer
function viewdas(dDAS; cm=:RdBu_11, climit=10000)
    GLMakie.activate!() #activate GL plottinge
    fig = Figure(; size=(1000, 600))
    ihm = iheatmap(fig[1, 1], dDAS.time, dDAS.offset, dDAS.data, colormap=cm, colorrange=(-climit, climit))
    Colorbar(fig[1, 2], ihm.plot, label="StrainRate [nm/m/s]")
    ihm.axis.xlabel = "Relative Time [s]"
    ihm.axis.ylabel = "Offset [m]"
    ihm.axis.title = dDAS.name
    DataInspector(fig)
    display(fig)
    return fig
end


## Quicker DAS viewer
function fig4(dDAS, specM; x=1, climit=10000)
    GLMakie.activate!() #activate GL plottinge
    chind = argmin(abs.(dDAS.offset .- x))
    cm = :RdBu_11
    fig = Figure(; size=(1200, 1000))
    ihm1 = iheatmap(fig[1:2, 1], dDAS.time, dDAS.offset, dDAS.data, colormap=cm, colorrange=(-climit, climit))
    Colorbar(fig[3, 1], ihm1.plot, label="StrainRate [nm/m/s]", vertical=false)
    ihm2 = iheatmap(fig[1:2, 2], 1:size(specM, 1), dDAS.offset, specM, colormap=:turbo)
    Colorbar(fig[3, 2], ihm2.plot, label="Spectral Amp", vertical=false)
    linkyaxes!(ihm1.axis, ihm2.axis)
    ihm1.axis.xlabel = "Relative Time [s]"
    ihm1.axis.ylabel = "Offset [m]"
    ihm2.axis.xlabel = "Frequency [Hz]"
    ihm2.axis.ylabel = "Offset [m]"
    p3 = ilines(fig[4, 1], dDAS.time, dDAS.data[:, chind])
    p4 = ilines(fig[4, 2], 1:size(specM, 1), specM[:, chind])
    p3.axis.xlabel = "Relative Time [s]"
    p3.axis.ylabel = "ηε"
    p4.axis.xlabel = "Frequency [Hz]"
    p4.axis.ylabel = "Amp"
    linkxaxes!(ihm1.axis, p3.axis)
    linkxaxes!(ihm2.axis, p4.axis)
    display(fig)

    return fig
end



#a function to safe the current figure
function savefig(fig, name)
    #using CairoMakie
    save("$name", fig, backend=CairoMakie)

    return
end


## Show a DAS Channel in time (this is constant X all times)
function viewchannel(dDAS; x, type=identity)

    GLMakie.activate!() #activate GL plottinge
    chind = argmin(abs.(dDAS.offset .- x)) # get the channel index
    printstyled("\n Plotting channel $x at index $chind \n", color=:green)
    ich = dDAS.offset[chind]
    #get the spect of the channel
    delta = dDAS.time[2] - dDAS.time[1]
    sr = Int(round(1 / delta))
    wl = Int(sr * 2)
    Σ = spectra(dDAS.data[:, chind], sr, sr; func=type)
    spectrum = Σ.y
    replace!(spectrum, NaN => 0.0)
    replace!(spectrum, Inf => 0.0)
    replace!(spectrum, -Inf => 0.0)
    #get the spectogram of the channel

    n = length(dDAS.data[:, chind])
    #nw = n ÷ 400
    #spec = spectrogram(dDAS.data[:, chind], nw, nw ÷ 2; fs=sr)

    spec = gspectrogram(dDAS.data[:, chind], sr, window_size=max(nextpow(2, n/300),512), overlap=0.5, alpha=2.5)

    spectogram = transpose(pow2db.(spec.amp))
    fig = Figure(; size=(600, 1000))

    p1 = ilines(fig[1, 1], dDAS.time, dDAS.data[:, chind])
    p2 = ilines(fig[3, 1], 1:size(spectrum, 1), spectrum)
    p3 = heatmap(fig[2, 1], spec.time, spec.freq, spectogram, colormap=:inferno, interpolate=true)
    Colorbar(fig[2, 2], p3.plot, label="dB")

    #p3 = contourf(fig[2, 1], spec.time, spec.freq[1]:(spec.freq[2]-spec.freq[1]):spec.freq[end], transpose(pow2db.(spec.power)), colormap=:inferno)

    p1.axis.xlabel = "Relative Time [s]"
    p1.axis.ylabel = "ηε"
    p1.axis.title = "Channel $chind at offset $ich"
    p2.axis.xlabel = "Frequency [Hz]"
    p2.axis.ylabel = "Amp"
    p3.axis.title = "Frequency Spectrum"
    p3.axis.ylabel = "Frequency [Hz]"
    p3.axis.xlabel = "Relative Time [s]"
    p2.axis.title = "Spectogram"
    linkxaxes!(p1.axis, p3.axis)
    #linkxaxes!(ihm2.axis, p4.axis)
    display(fig)
    return channelx, spectrum, spectogram, fig
end



function recordsection(dDAS; scale=1.0)

GLMakie.activate!() #activate GL plotting
dx = dDAS.offset[2]-dDAS.offset[1]
n = length(dDAS.time)
scaleddata = scale*dx*dDAS.data 
scaleddata = scaleddata./maximum(abs.(scaleddata))
nc = size(dDAS.data,2)
# Plot signals with positive/negative areas filled:

fig = Figure(; size=(1000, 1000))
ax = Axis(fig[1,1], title = dDAS.name, xlabel = "Time [s]", ylabel = "Offset [m]", xminorgridvisible = true, yminorgridvisible = true, limits = (nothing, nothing))
    for i in 1:nc
        #band!(ax, time, map(x-> x < 0 ? x + i*dx : i*dx, strainrate[i,:]), i*dx, color=("red", 0.3))
        #band!(time, map(x-> x >= 0 ? x + i*dx : i*dx, strainrate[i,:]), i*dx, color=("blue", 0.3))
        ilines!(ax, dDAS.time, scaleddata[:,i] .+ i*dx .+dDAS.offset[1], linewidth=0.5, color=:black)
    end
    tightlimits!(ax)
    fig
    return fig
end






export viewdas, fig4, function, savefig, viewchannel, recordsection
