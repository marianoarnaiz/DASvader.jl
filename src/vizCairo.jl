"This is the Visualization  part of DAS VADER V1.0"

using MakieThemes, CairoMakie
Makie.set_theme!(ggthemr(:fresh)) #Activate the cool theme :)

##################################################
println("Using CairoMakie (Headless)")
##################################################

## Quicker DAS viewer
function viewdas(dDAS; cm=:RdBu_11, climit=10000)
    fig = Figure(; size=(1000, 600))
    ihm = iheatmap(fig[1, 1], dDAS.time, dDAS.offset, dDAS.data, colormap=cm, colorrange=(-climit, climit))
    Colorbar(fig[1, 2], ihm.plot, label="StrainRate [nm/m/s]")
    ihm.axis.xlabel = "Relative Time [s]"
    ihm.axis.ylabel = "Offset [m]"
    DataInspector(fig)
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

    chind = argmin(abs.(dDAS.offset .- x)) # get the channel index
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
    channelx = dDAS.data[:, chind]
    #get the spectogram of the channel
    n = length(dDAS.data[:, chind])
    nw = n ÷ 400
    spec = spectrogram(dDAS.data[:, chind], nw, nw ÷ 2; fs=sr)
    spectogram = transpose(pow2db.(spec.power))
    fig = Figure(; size=(600, 1000))

    p1 = ilines(fig[1, 1], dDAS.time, dDAS.data[:, chind])
    p2 = ilines(fig[3, 1], 1:size(spectrum, 1), spectrum)
    p3 = iheatmap(fig[2, 1], spec.time, spec.freq[1]:(spec.freq[2]-spec.freq[1]):spec.freq[end], spectogram, colormap=:inferno)
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


export viewdas, viewchannel, savefig
