"This is the Space & wavenumber (read, write, etc) part of DAS VADER V1.0"


## 1. Using: Loading the package to be used
using FourierAnalysis, Statistics, DSP, Seis

"""
xppdas: Preprocess dDAS data for each channel(remove mean, remove trend, taper)

Input:
- dDAS: DAS data structure
- rmean: remove mean flag (true/false). Default is true.
- rtrend: remove trend flag (true/false). Default is true.
- taper: taper the ends of the signal flag (true/false). Default is true.
- f (form) may be :hanning, :hamming or :cosine (include the :).
- w (width) represents the fraction (at both ends) of the trace tapered, up to 0.5.

Outputs:
- ppDAS: DAS data preprocessed

Notes:
 - You can turn off any of the processing steps.
 - No figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xppdas(dDAS; rmean=true, rtrend=true, taper=false, w=0.05, f=:hanning)
```
"""
function xppdas(dDAS; rmean=true, rtrend=true, taper=true, w=0.05, f=:hanning)
    # Begin
    ppDAS = deepcopy(dDAS)
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    #Preprocessing Loop
    @inbounds @simd for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        if rmean == true
            remove_mean!(dastrace) #remove the mean value
        end
        if rtrend == true
            remove_trend!(dastrace) #remove a linear trend from the data
        end
        if taper == true
            taper!(dastrace, w; form=f) # taper the ends with the trace
        end
        ppDAS.data[i, :] = trace(dastrace)
    end
    return ppDAS
end
#############################################################################################


"""
xppdas!: Preprocess IN PLACE dDAS data for each channel(remove mean, remove trend, taper)

Input:
- dDAS: DAS data structure
- rmean: remove mean flag (true/false). Default is true.
- rtrend: remove trend flag (true/false). Default is true.
- taper: taper the ends of the signal flag (true/false). Default is true.
?

Outputs:
- dDAS: dDAS data preprocessed

Notes:
 - You can turn off any of the processing steps.
 - No figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xppdas!(dDAS; rmean=true, rtrend=true, taper=false, width=0.05, form=:hanning)
```
"""
function xppdas!(dDAS; rmean=true, rtrend=true, taper=true, w=0.05, f=:hanning)
    # Begin
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    #Preprocessing Loop
    @inbounds @simd for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        if rmean == true
            remove_mean!(dastrace) #remove the mean value
        end
        if rtrend == true
            remove_trend!(dastrace) #remove a linear trend from the data
        end
        if taper == true
            taper!(dastrace, w; form=f) # taper the ends with the trace
        end
        dDAS.data[i, :] = trace(dastrace)
    end
    return
end
#############################################################################################

"""
xbpdas: Band Pass dDAS data for each channel.

Input:
 - dDAS: DAS data structure
 - w1: lower wavelenght of the filter. In Hz.
 - w2: upper wavelenght of the filter. In Hz.
 - poles: poles of the filter.

Outputs:
- bpDAS: band passed DAS data

Notes:
 - No figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xbpDAS = bpdas(dDAS; f1=50, f2=500, poles=4)
```
"""
function xbpdas(dDAS; w1=1, w2=500, poles=4, twopass=false)
    # Begin
    bpDAS = deepcopy(dDAS)
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    # Bandpass Filtering Loop
    for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        bandpass!(dastrace, w1, w2; poles, twopass)
        bpDAS.data[i, :] = trace(dastrace)
    end
    return bpDAS
end
#############################################################################################

"""
xbpdas!: Band Pass IN PLACE dDAS data for each channel.

Input:
 - dDAS: DAS data structure
 - w1: lower wavelenght of the filter. In Hz.
 - w2: upper wavelenght of the filter. In Hz.
 - poles: poles of the filter.

Outputs:
- dDAS filtered

Notes:
 - No figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xbpdas!(dDAS, f1=50, f2=500, poles=4)
```
"""
function xbpdas!(dDAS; w1=1, w2=500, poles=4, twopass=false)
    # Begin
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    # Bandpass Filtering Loop
    for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        bandpass!(dastrace, w1, w2; poles, twopass)
        dDAS.data[i, :] = trace(dastrace)
    end
    return
end
#############################################################################################

"""
xlpdas: low Pass dDAS data for each channel.

Input:
 - dDAS: DAS data structure
 - w1: wavelenght of the filter. In Hz.
 - poles: poles of the filter.

Outputs:
- lpDAS: low passed DAS data

Notes:
 - No figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xlpDAS = lpdas(dDAS; f1=50, poles=4)
```
"""
function xlpdas(dDAS; w1=1, poles=4, twopass=false)
    # Begin
    lpDAS = deepcopy(dDAS)
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    # Bandpass Filtering Loop
    for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        lowpass!(dastrace, w1; poles, twopass)
        lpDAS.data[i, :] = trace(dastrace)
    end
    return lpDAS
end
#############################################################################################

"""
xlpdas!: Low Pass IN PLACE dDAS data for each channel.

Input:
 - dDAS: DAS data structure
 - w1: wavelenght of the filter. In Hz.
 - poles: poles of the filter.

Outputs:
- lpDAS: low passed DAS data

Notes:
 - No figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xlpdas!(dDAS; f1=50, poles=4)
```
"""
function xlpdas!(dDAS; w1=1, poles=4, twopass=false)
    # Begin
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    # Bandpass Filtering Loop
    for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        lowpass!(dastrace, w1; poles, twopass)
        dDAS.data[i, :] = trace(dastrace)
    end
    return
end
#############################################################################################

"""
xhpdas: High Pass dDAS data for each channel.

Input:
 - dDAS: DAS data structure
 - w1: wavelenght of the filter. In Hz.
 - poles: poles of the filter.

Outputs:
- hpDAS: high passed DAS data

Notes:
 - No figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xhpDAS = hpdas(dDAS; w1=1, poles=4, twopass=false)
```
"""
function xhpdas(dDAS; w1=100, poles=4, twopass=false)
    # Begin
    hpDAS = deepcopy(dDAS)
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    # Bandpass Filtering Loop
    for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        highpass!(dastrace, w1; poles, twopass)
        hpDAS.data[i, :] = trace(dastrace)
    end
    return hpDAS
end
#############################################################################################

"""
xhpdas!: High Pass IN PLACE dDAS data for each channel.

Input:
 - dDAS: DAS data structure
 - w1: frequency of the filter. In Hz.
 - poles: poles of the filter.

Outputs:
- hpDAS: high passed DAS data

Notes:
 - No figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xhpdas!(dDAS; w1=1, poles=4, twopass=false)
```
"""
function xhpdas!(dDAS; w1=100, poles=4, twopass=false)
    # Begin
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    # Bandpass Filtering Loop
    for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        highpass!(dastrace, w1; poles, twopass)
        dDAS.data[i, :] = trace(dastrace)
    end
    return
end
#############################################################################################

"""
xnormdas: Normalize dDAS

Input:
 - dDAS: DAS data structure
 - style: normalize all the array ("all") or by each channel ("channel")

Outputs:
- nDAS: normalized DAS data

Notes:
 - No figure is created.
 - We choose to find the max(abs(data)) meaning that no value greater than 1 will be produced
 - "all" option normalizes the entire matrix to 1.
 - "channel" option normalizes each channel independently.
# Example: Save one channel and its time vector to a file.
```
julia> nDAS = xnormdas(dDAS; style="all")
```
"""
function xnormdas(dDAS; style="all")

    nDAS = deepcopy(dDAS)

    if style == "all"
        m = maximum(abs.(dDAS.data)) #get largest value in array
        nDAS.data = nDAS.data ./ m #divide by m

    elseif style == "channel"
        noc = size(dDAS.time, 1)
        for i = 1:size(dDAS.time, 1)
            printstyled(" Working on time sample $i of $noc.\n", color=43)
            m = maximum(abs.(dDAS.data[i, :])) #get largest value in array
            nDAS.data[i, :] = nDAS.data[i, :] ./ m #divide by m
        end

    else
        printstyled(" Style not recognized c.\n", color=:red)
    end
    return nDAS
end
#############################################################################################


"""
xnormdas!: Normalize dDAS IN PLACE

Input:
 - dDAS: DAS data structure
 - style: normalize all the array ("all") or by each channel ("channel")

Outputs:
- nDAS: normalized DAS data

Notes:
 - No figure is created.
 - We choose to find the max(abs(data)) meaning that no value greater than 1 will be produced
 - "all" option normalizes the entire matrix to 1.
 - "channel" option normalizes each channel independently.
# Example: Save one channel and its time vector to a file.
```
julia> xnormdas!(dDAS; style="all")
```
"""
function xnormdas!(dDAS; style="all")


    if style == "all"
        m = maximum(abs.(dDAS.data)) #get largest value in array
        dDAS.data = dDAS.data ./ m #divide by m

    elseif style == "channel"
        noc = size(dDAS.time, 1)

        for i = 1:size(dDAS.time, 1)
            printstyled(" Working on time sample $i of $noc.\n", color=43)
            m = maximum(abs.(dDAS.data[i, :])) #get largest value in array
            dDAS.data[i, :] = dDAS.data[i, :] ./ m #divide by m
        end

    else
        printstyled(" Style not recognized c.\n", color=:red)
    end
    return
end

##########################################################


## Decimate all the DAS signals in time axis
function xdecimatedas(dDAS; n=2, antialias=true)
    # Make a copy
    unDAS = deepcopy(dDAS)
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    #compute the new size of strain rate
    dastrace = Trace(b, delta, dDAS.data[1, :])
    decimate!(dastrace, n; antialias=antialias)
    unDAS.data = zeros(noc, size(trace(dastrace), 1))
    # new times
    unDAS.offset = dDAS.offset[1:n:end]

    # Decimating Loop
    for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        decimate!(dastrace, n; antialias=antialias)
        unDAS.data[i, :] = trace(dastrace)
    end
    return unDAS
end


## Decimate all the DAS signals in time axis IN PLACE
function xdecimatedas!(dDAS; n=2, antialias=true)
    #Set up some trace information
    b = 0.0
    delta = dDAS.offset[2] - dDAS.offset[1]
    noc = size(dDAS.time, 1)
    #compute the new size of strain rate
    dastrace = Trace(b, delta, dDAS.data[1, :])
    decimate!(dastrace, n; antialias=antialias)
    holdit = zeros(noc, size(trace(dastrace), 1))
    # new times
    dDAS.offset = dDAS.offset[1:n:end]

    # Decimating Loop
    for i = 1:size(dDAS.time, 1)
        printstyled(" Working on time sample $i of $noc.\n", color=43)
        dastrace = Trace(b, delta, dDAS.data[i, :]) #write the data to trace format
        decimate!(dastrace, n; antialias=antialias)
        holdit[i, :] = trace(dastrace)
    end
    dDAS.data = holdit
    return
end


#############################################################################################

"""
xdasspectogram: show the average spectogram

Input:
 - dDAS: DAS data structure

Outputs:
- av_spc: average spectogram matrix of the das

Notes:
 - A figure is created.
# Example: Save one channel and its time vector to a file.
```
julia> xdasspectogram(dDAS)
```
"""
function xspectogramdas(dDAS; nsamples=50, noverlap=2)

    meandas = mean(dDAS.data, dims=1)[:]
    delta = dDAS.offset[2] - dDAS.offset[1]
    sr = Int(round(1 / delta))
    n = length(dDAS.data[1, :])
    nw = n ÷ nsamples
    spec = spectrogram(dDAS.data[:, 1], nw, nw ÷ noverlap; fs=sr)
    spectogram = transpose(pow2db.(spec.power))
    spectdas = zeros(size(spectogram, 1), size(spectogram, 2), size(dDAS.offset, 1))

    for i = 1:size(dDAS.offset, 1)
        printstyled(" Spectogram of Channel $i.\n", color=:cyan)
        specf = spectrogram(dDAS.data[:, i], nw, nw ÷ noverlap; fs=sr)
        spectdas[:, :, i] = transpose(pow2db.(specf.power))
    end
    spectdas = mean(spectdas, dims=3)[:, :, 1]

    fig = Figure(; size=(800, 800))
    p1 = iheatmap(fig[1:2, 1], spec.time, spec.freq[1]:(spec.freq[2]-spec.freq[1]):spec.freq[end], spectogram, colormap=:inferno)
    Colorbar(fig[1:2, 2], p1.plot, label="dB")
    p2 = ilines(fig[3, 1], dDAS.time, meandas)

    p1.axis.xlabel = "Relative Time [s]"
    p1.axis.ylabel = "Frequency [Hz]"
    p1.axis.title = "Average Spectrogram"
    p2.axis.xlabel = "Relative Time [s]"
    p2.axis.ylabel = "ηε"
    p2.axis.title = "Average DAS data"
    linkxaxes!(p1.axis, p2.axis)

    #linkxaxes!(ihm2.axis, p4.axis)
    display(fig)
    return meandas, spectdas, fig

end

#############################################################################################

"""
dasfspec: compute the frequency spectrum for each channel

Input:
 - dDAS: DAS data structure
 type can be:
 nothing: power spectra are extracted;
 sqrt: will extract amplitude spectra,
 log: will extract log-spectra,
 decibel: will extract spectra in deciBels (see decibel).

Outputs:
- spc: frequency spectrum matrix of the das

Notes:
# Example: compute the spectrum for each channel.
```
julia> dasfspec(dDAS)
```
"""
function dasfspec(dDAS; type=identity)

    #Set up some trace information
    delta = dDAS.time[2] - dDAS.time[1]
    sr = Int(round(1 / delta))
    wl = Int(sr * 2)
    global specM = zeros(Int(sr / 2), size(dDAS.offset, 1))
    noc = size(dDAS.offset, 1)
    for i = 1:size(dDAS.offset, 1)
        printstyled(" Working on channel $i of $noc.\n", color=:blue)
        Σ = spectra(dDAS.data[:, i], sr, sr; func=type)
        global specM[:, i] = Σ.y
    end
    replace!(specM, NaN => 0.0)
    replace!(specM, Inf => 0.0)
    replace!(specM, -Inf => 0.0)


    GLMakie.activate!() #activate GL plottinge
    fig = Figure(; size=(600, 700))
    ihm2 = iheatmap(fig[1:2, 1], specM, colormap=:turbo)
    ihm2.axis.xlabel = "Frequency [Hz]"
    ihm2.axis.ylabel = "Offset [m]"
    Colorbar(fig[1:2, 2], ihm2.plot, label="Spectral Amp")
    p1 = ilines(fig[3, 1], mean(specM, dims=2)[:])
    p2 = ilines!(fig[3, 1], minimum(specM, dims=2)[:])
    p3 = ilines!(fig[3, 1], maximum(specM, dims=2)[:])

    smean = mean(specM, dims=2)[:]

    return smean, specM, fig
end

#############################################################################################

"""
xdasmin: find the average value for each channel and find the min

Input:
 - dDAS: DAS data structure

Outputs:
- No output

Notes:
# Example: compute the spectrum for each channel.
```
julia> xdasmin(dDAS)
```
"""
function xdasmin(dDAS)


    #get them mean and median
    mea = mean(dDAS.data, dims=2)[:]
    med = median(dDAS.data, dims=2)[:]
    #get SNR
    noise = mean(abs.(dDAS.data), dims=2)[:]
    signal = maximum(abs.(dDAS.data), dims=2)[:]
    snr = signal ./ noise
    replace!(snr, NaN => 0.0)

    #print some results
    imea = argmin(mea)
    ichan = dDAS.offset[imea]
    imed = argmin(med)
    ichad = dDAS.offset[imed]
    printstyled("\n Channel with lowest mean is $imea at $ichan m from 0.0 \n", color=:white)
    printstyled("\n Channel with lowest median is $imed at $ichad m from 0.0 \n", color=:white)

    #make a figure
    GLMakie.activate!() #activate GL plottinge
    fig = Figure()
    ax1 = Axis(fig[1, 1], xlabel="Value", ylabel="Offset [m]", title="Mean and Median")
    lines!(ax1, mea, dDAS.offset, label="Mean")
    lines!(ax1, med, dDAS.offset, label="Median")
    axislegend()
    ax2 = Axis(fig[1, 2], xlabel="SNR", title="SNR")
    lines!(ax2, snr, dDAS.offset)
    display(fig)

    return fig
end


"""
xcat: Concatenate 2 strainrate matrixes in distance assuming that the second begins right after the end of the first. Notice that sizes must be compatible!

Input:

- dDAS1: first file to cat .
- dDAS2: second file to cat .

Outputs:

- dDAS3: cat-ed file.


Notes:
 - No processing is applied to the data.
 - No figure is created.
 - Catting strainrate matrixes can easily overload your RAM memory. Proceed with caution.
# Example: Concatenate strainrate1 followed by strainrate2.
```
julia> dDAS3 = xcat(dDAS1,dDAS2);
```
"""
function xcat(dDAS1, dDAS2)
    if size(dDAS1.offset) == size(dDAS2.offset) && dDAS1.time[2] - dDAS1.time[1] == dDAS2.time[2] - dDAS2.time[1]
        #Cat the Strain Rates Matrixes
        cat_data = [dDAS1.data dDAS2.data]
        #Cat Human time
        utime = range(start=datetime2unix(dDAS1.htime[1]), step=(dDAS1.time[2] - dDAS1.time[1]), length=(size(dDAS1.time)[1] + size(dDAS2.time)[1]))
        cat_htime = unix2datetime.(utime)
        #cat_htime=[htime1 htime2]
        #Cat Time
        cat_time = range(start=dDAS1.time[1], step=(dDAS1.time[2] - dDAS1.time[1]), length=(size(dDAS1.time)[1] + size(dDAS2.time)[1]))
        cat_offset = copy(dDAS1.offset)

        atrib = attb(dDAS1.atrib.AmpliPower, dDAS1.atrib.BlockRate, dDAS1.atrib.Components, dDAS1.atrib.DataDomain, dDAS1.atrib.DerivationTime, dDAS1.atrib.Extent, dDAS1.atrib.FiberLength, dDAS1.atrib.GaugeLength, dDAS1.atrib.Hostname, dDAS1.atrib.Origin, dDAS1.atrib.Oversampling, dDAS1.atrib.PipelineTracker, dDAS1.atrib.PulseRateFreq, dDAS1.atrib.PulseWidth, dDAS1.atrib.SamplingRate, dDAS1.atrib.SamplingRes, dDAS1.atrib.Spacing, dDAS1.atrib.BlockOverlap)

        dDAS3 = iDAS(cat_data, cat_time, cat_htime, cat_offset, atrib)


    else
        printstyled("Not The same Number of Channels. Can´t xcat!", color=:red, bold=true)
    end
    return dDAS3
end


## Take the integral of  all the DAS signals in time axis
function intdas(dDAS; method=:trapezium)
    # Begin
    unDAS = deepcopy(dDAS)
    #Set up some trace information
    b = 0.0
    delta = dDAS.time[2] - dDAS.time[1]
    noc = size(dDAS.offset, 1)
    # Integration loop
    for i = 1:size(dDAS.offset, 1)
        printstyled(" Integrating channel $i of $noc.\n", color=:yellow)
        dastrace = Trace(b, delta, dDAS.data[:, i]) #write the data to trace format
        integrate!(dastrace, method)
        unDAS.data[1:size(trace(dastrace), 1), i] = trace(dastrace)
    end
    return unDAS
end


## Take the integral in place of  all the DAS signals in time axis
function intdas!(dDAS; method=:trapezium)
    # Begin
    #Set up some trace information
    b = 0.0
    delta = dDAS.time[2] - dDAS.time[1]
    noc = size(dDAS.offset, 1)
    # Integration loop
    for i = 1:size(dDAS.offset, 1)
        printstyled(" Integrating channel $i of $noc.\n", color=:yellow)
        dastrace = Trace(b, delta, dDAS.data[:, i]) #write the data to trace format
        integrate!(dastrace, method)
        dDAS.data[1:size(trace(dastrace), 1), i] = trace(dastrace)
    end
    return dDAS
end


## Take the integral of  all the DAS signals in time axis
function diffdas(dDAS; points=2)
    # Begin
    unDAS = deepcopy(dDAS)
    #Set up some trace information
    b = 0.0
    delta = dDAS.time[2] - dDAS.time[1]
    noc = size(dDAS.offset, 1)
    # Integration loop
    for i = 1:size(dDAS.offset, 1)
        printstyled(" Differentiating channel $i of $noc.\n", color=:yellow)
        dastrace = Trace(b, delta, dDAS.data[:, i]) #write the data to trace format
        differentiate!(dastrace; points=points)
        unDAS.data[1:size(trace(dastrace), 1), i] = trace(dastrace)
    end
    return unDAS
end


## Take the integral in place of  all the DAS signals in time axis
function diffdas!(dDAS; points=2)
    # Begin
    #Set up some trace information
    b = 0.0
    delta = dDAS.time[2] - dDAS.time[1]
    noc = size(dDAS.offset, 1)
    # Integration loop
    for i = 1:size(dDAS.offset, 1)
        printstyled(" Differentiating channel $i of $noc.\n", color=:yellow)
        dastrace = Trace(b, delta, dDAS.data[:, i]) #write the data to trace format
        differentiate!(dastrace; points=points)
        dDAS.data[1:size(trace(dastrace), 1), i] = trace(dastrace)
    end
    return dDAS
end

export ppdas, ppdas!, bpdas, bpdas!, lpdas, lpdas!, hpdas, hpdas!, normdas, normdas!, envdas, envdas!, decimatedas, decimatedas!, dasfspec, timecat, intdas, intdas!
