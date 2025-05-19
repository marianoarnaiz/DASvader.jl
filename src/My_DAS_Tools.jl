"My_DAS_Tools.jl"

## Show a DAS Channel in time (this is constant X all times)
function viewchannel(strainrate, time, offset; x)
    channel = argmin(abs.(offset .- x))
    schannel = offset[channel]
    printstyled(" Plotting Channel $channel at distance $schannel m.\n", color=:yellow)
    #fig=Figure()
    c, viz = InteractiveViz.iplot(time, strainrate[channel, :]; grid=true, cursor=true, ylabel="Strain Rate [nm/m/s]", xlabel="Time [s]")
    return c, viz
end

## DAS Channel 2 SAC format with basic information

## Show a DAS Channel in time (this is constant X all times)
function channel2sac(strainrate, time, htime, offset; x)

    if typeof(x) == Int64 || typeof(x) == Float64
        channel = argmin(abs.(offset .- x))
        schannel = offset[channel]
        printstyled(" Saving Channel $channel at distance $schannel m to SAC file\n", color=:yellow)
        signal = strainrate[channel, :]
        t = Trace(0.0, time[2] - time[1], signal)
        t.evt.time = DateTime(year(htime[1]), month(htime[1]), day(htime[1]), hour(htime[1]), minute(htime[1]), second(htime[1]), millisecond(htime[1]))
        write_sac(t, "x=$schannel.sac")
    end

    if typeof(x) == UnitRange{Int64} || typeof(x) == StepRange{Int64,Int64} || typeof(x) == StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64} || typeof(x) == StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64}

        for i = 1:size(x, 1)
            channel = argmin(abs.(offset .- x[i]))
            schannel = offset[channel]
            printstyled(" Saving Channel $channel at distance $schannel m to SAC file\n", color=:yellow)
            signal = strainrate[channel, :]
            t = Trace(0.0, time[2] - time[1], signal)
            t.evt.time = DateTime(year(htime[1]), month(htime[1]), day(htime[1]), hour(htime[1]), minute(htime[1]), second(htime[1]), millisecond(htime[1]))
            write_sac(t, "x=$schannel.sac")
        end
    end


end




## Spectogram View and Write!
function channelspectogram(strainrate, time, offset; x, cm=:inferno)
    channel = argmin(abs.(offset .- x))
    schannel = offset[channel]
    printstyled(" Spectogram of Channel $channel at distance $schannel m.\n", color=:yellow)
    channeldata = slice_strainrate[channel, :]
    n = length(channeldata)
    nw = n ÷ 25
    fs1 = 1 / (slice_time[2] - slice_time[1])
    spec = spectrogram(channeldata, nw, nw ÷ 2; fs=fs1)
    c, viz = InteractiveViz.iheatmap(transpose(pow2db.(spec.power)), minimum(spec.time), maximum(spec.time), minimum(spec.freq), maximum(spec.freq); pooling=mean, cursor=true, colormap=cm, ylabel="Frequency [Hz]", xlabel="Time [s]")
    Colorbar(viz.scene; bbox=BBox(100, 0, 10, 0), spinewidth=1, height=20, width=500, vertical=false, label="Power [dB]", colormap=cm, halign=:left, valign=:bottom)
    return c, viz
end

## Spectogram View and Write!
function dasspectogram(strainrate, time, offset; x, cm=:default)
    channel = argmin(abs.(offset .- x))
    schannel = offset[channel]
    printstyled(" Spectogram of Channel $channel at distance $schannel m.\n", color=:yellow)
    channeldata = slice_strainrate[channel, :]
    n = length(channeldata)
    nw = n ÷ 50
    fs = 1 / (slice_time[2] - slice_time[1])
    spec = spectrogram(channeldata, nw, nw ÷ 10; fs=fs)
    c, viz = InteractiveViz.iheatmap(transpose(pow2db.(spec.power)), minimum(spec.time), maximum(spec.time), minimum(spec.freq), maximum(spec.freq); pooling=maximum, cursor=true, colormap=cm, ylabel="Frequency [Hz]", xlabel="Time [s]")
    Colorbar(viz.scene; bbox=BBox(100, 0, 10, 0), spinewidth=1, height=20, width=500, vertical=false, label="Power [dB]", colormap=cm, halign=:left, valign=:bottom)
    return c, viz
end

# Plot Spectrums. TO CHECK!!!!!!!!
function spectrum(strainrate, time; fformat="png")
    tsr = reshape(strainrate, size(strainrate, 2), size(strainrate, 1))
    sampling_rate = Int64(1 / (time[2] - time[1]))
    window_lenght = sampling_rate
    srS = spectra(tsr, sampling_rate, window_lenght, tapering=riesz, smoothing=blackmanSmoother, func=decibel)
    p1 = Plots.plot(srS.flabels, srS.y, lw=0.5, xaxis=:log10, xlabel="Frequency [Hz]", ylabel="dB (All Channels)", minorgrid=true, legend=false)
    p2 = Plots.plot(srS.flabels, mean(srS.y, dims=2), lw=0.5, lc=:black, xaxis=:log10, xlabel="Frequency [Hz]", ylabel="dB (Mean Spec.)", minorgrid=true, legend=false)
    Plots.plot(p1, p2, layout=(2, 1))
    savefig("Spectrums.$fformat")
end
