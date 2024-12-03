"This is the Detection (sta/lta, MER, HOS, etc) part of DAS VADER V1.0"


######################################################################
# FUNCTIONS FOR DETECTION IN DAS!
######################################################################

function staltadas(dDAS; stw=0.002, ltw=0.05, thr=3)
#Parameters
dt = dDAS.time[2] - dDAS.time[1]  # dt
ns = Int64(floor(stw / dt))       # no. of data points in the short time window
nl = Int64(floor(ltw / dt))       # no. of data points in the long time window
nt = size(dDAS.data, 1)           # number of time samples
nc = size(dDAS.data, 2)           # number of channels

# STA/LTA ratio matrix (transposed to [nt, nc] for column-major order)
sra = zeros(nt, nc)  # STA/LTA ratio, initiate variable

# For each channel
for c = 1:nc
    printstyled(" STA/LTA on channel $c of $nc.\n", color=32)

    # Initialize the first `sta` and `lta` windows
    sta_sum = sum(abs.(dDAS.data[1:ns, c]))  # Sum for short window
    lta_sum = sum(abs.(dDAS.data[1:nl, c]))  # Sum for long window

    # Calculate initial STA/LTA ratio
    sta = sta_sum / ns
    lta = lta_sum / nl
    sra[1, c] = sta / lta

    # Compute STA/LTA for the rest of the samples
    for k = 2:nt
        # Update `sta_sum` by adding new element and removing the old one
        if k + ns - 1 <= nt
            sta_sum += abs(dDAS.data[k + ns - 1, c]) - abs(dDAS.data[k - 1, c])
            sta = sta_sum / ns
        end

        # Update `lta_sum` similarly
        if k + nl - 1 <= nt
            lta_sum += abs(dDAS.data[k + nl - 1, c]) - abs(dDAS.data[k - 1, c])
            lta = lta_sum / nl
        end

        # Calculate STA/LTA ratio
        sra[k, c] = sta / lta
    end
end

sra[isnan.(sra)] .= 0.0; #replace any NAN with zero. Delete for speed


    #sra = sra
    #stack trhe strainrate too
    Σsra = sum(sra, dims=2)[:] ./ nc

    #find all points higher than threshold and print message to user
    ipeaks = findpeaks(Σsra, dDAS.time, min_height=Float64(thr), min_dist=0.001)
    event_count = size(ipeaks, 1)
    peaks = [dDAS.time[ipeaks] dDAS.htime[ipeaks] ipeaks]

    #stack trhe strainrate too
    Σstrainrate = sum(dDAS.data, dims=2)[:] / nc

    if event_count > 0
        printstyled(" STA/LTA detected $event_count events.\n", color=:green)
        #Make a plot
        fig = Figure()
        hm1 = iheatmap(fig[1:2, 1], dDAS.time, dDAS.offset, sra, colormap=:amp, colorrange=(0, round(maximum(sra) * 0.7)))
        hm1.axis.ylabel = "Offset [m]"
        Colorbar(fig[1:2, 2], hm1.plot, label="Sta/Lta")
        hm1.axis.title = "Sta/Lta.\n"
        p2 = ilines(fig[3, 1], dDAS.time, Σstrainrate, color=:grey)
        scatter!(fig[3, 1], dDAS.time[ipeaks], Σstrainrate[ipeaks], marker=:dtriangle, markersize=10, color=:black)
        p2.axis.ylabel = "ηε"
        p3 = ilines(fig[4, 1], dDAS.time, Σsra, color=:darkred)
        scatter!(fig[4, 1], dDAS.time[ipeaks], Σsra[ipeaks], marker=:dtriangle, markersize=10, color=:black)
        p3.axis.xlabel = "Relative Time [s]"
        p3.axis.ylabel = "Mean Sta/Lta"
        linkxaxes!(hm1.axis, p2.axis, p3.axis)
        DataInspector(fig)
        display(fig)

    else
        printstyled(" STA/LTA detected no events.\n", color=:red)

        fig = Figure()
        hm1 = iheatmap(fig[1:2, 1], dDAS.time, dDAS.offset, sra, colormap=:amp, colorrange=(0, round(maximum(sra) * 0.7)))
        hm1.axis.ylabel = "Offset [m]"
        Colorbar(fig[1:2, 2], hm1.plot, label="Sta/Lta")
        hm1.axis.title = "Sta/Lta.\n"
        p2 = ilines(fig[3, 1], dDAS.time, Σstrainrate, color=:grey)
        p2.axis.ylabel = "ηε"
        p3 = ilines(fig[4, 1], dDAS.time, Σsra, color=:darkred)
        p3.axis.xlabel = "Relative Time [s]"
        p3.axis.ylabel = "Mean Sta/Lta"
        linkxaxes!(hm1.axis, p2.axis, p3.axis)
        DataInspector(fig)
        display(fig)

    end
    return sra, Σsra, peaks, fig
end



function dasmer(dDAS; ltw=0.01, thr=0.8)
# Parameters
dt = dDAS.time[2] - dDAS.time[1]     # dt
nl = Int64(floor(ltw / dt))          # number of data points in the long time window
nt = size(dDAS.data, 1)              # number of time samples
nc = size(dDAS.data, 2)              # number of channels

# Initialize energy ratio and modified energy ratio matrices (nt x nc)
er = zeros(nt, nc)
mer = zeros(nt, nc)

# For each channel
for c = 1:nc
    printstyled("Working on channel $c of $nc.\n", color=:yellow)

    # Initialize sums for the first window
    sum_future = sum(dDAS.data[1:nl, c] .^ 2)
    sum_past = sum(dDAS.data[1:nl, c] .^ 2)

    # Calculate `er` for the remaining time samples
    for k = nl+1:nt-nl
        # Update rolling sums
        sum_future += dDAS.data[k + nl - 1, c]^2 - dDAS.data[k - 1, c]^2
        sum_past += dDAS.data[k, c]^2 - dDAS.data[k - nl, c]^2

        # Calculate energy ratio
        er[k, c] = sum_future / sum_past
    end
end

# Modified Energy Ratio (MER) calculation, normalized
for c = 1:nc
    mer[:, c] .= er[:, c] .* abs.(dDAS.data[:, c]) .^ 3
    mer[:, c] ./= maximum(mer[:, c])  # Normalization
end

mer[isnan.(mer)] .= 0.0; #replace any NAN with zero. Delete for speed

# Stack the mer results
Σmer = sum(mer, dims=2)[:, 1] / nc

     #find all points higher than threshold and print message to user
    ipeaks = findpeaks(Σmer, dDAS.time, min_height=Float64(thr), min_dist=0.001)
    event_count = size(ipeaks, 1)
    peaks = [dDAS.time[ipeaks] dDAS.htime[ipeaks] ipeaks]
    #stack trhe strainrate too
    Σstrainrate = sum(abs.(dDAS.data), dims=2)
    Σstrainrate = Σstrainrate[:] #because sometimes computational languages are dum!

    if event_count > 0
        printstyled(" MER detected $event_count events.\n", color=:green)

        fig = Figure()
        hm1 = iheatmap(fig[1:2, 1], dDAS.time, dDAS.offset, mer, colormap=Reverse(:ice), colorrange=(0, round(maximum(mer) * 0.7)))
        hm1.axis.ylabel = "Offset [m]"
        Colorbar(fig[1:2, 2], hm1.plot, label="MER")
        hm1.axis.title = "MER.\n"
        p2 = ilines(fig[3, 1], dDAS.time, Σstrainrate, color=:grey)
        scatter!(fig[3, 1], dDAS.time[ipeaks], Σstrainrate[ipeaks], marker=:dtriangle, markersize=10, color=:black)
        p2.axis.ylabel = "ηε"
        p3 = ilines(fig[4, 1], dDAS.time, Σmer, color=:darkblue)
        scatter!(fig[4, 1], dDAS.time[ipeaks], Σmer[ipeaks], marker=:dtriangle, markersize=10, color=:black)
        p3.axis.xlabel = "Relative Time [s]"
        p3.axis.ylabel = "Mean MER"
        linkxaxes!(hm1.axis, p2.axis, p3.axis)
        DataInspector(fig)
        display(fig)
   
    else
        printstyled(" MER detected no events.\n", color=:red)

        fig = Figure()
        hm1 = iheatmap(fig[1:2, 1], dDAS.time, dDAS.offset, mer, colormap=Reverse(:ice), colorrange=(0, round(maximum(mer) * 0.7)))
        hm1.axis.ylabel = "Offset [m]"
        Colorbar(fig[1:2, 2], hm1.plot, label="MER")
        hm1.axis.title = "MER.\n"
        p2 = ilines(fig[3, 1], dDAS.time, Σstrainrate, color=:grey)
        p2.axis.ylabel = "ηε"
        p3 = ilines(fig[4, 1], dDAS.time, Σmer, color=:darkblue)
        p3.axis.xlabel = "Relative Time [s]"
        p3.axis.ylabel = "Mean MER"
        linkxaxes!(hm1.axis, p2.axis, p3.axis)
        DataInspector(fig)
        display(fig)
    end
    return mer, Σmer, peaks, fig

end



function dasrechos(dDAS; C::Float64, order::Int64, var_min=-1, definition=0, thr=0.8)

    nt = size(dDAS.data,1) #number of time samples
    nc = size(dDAS.data,2) #number of channels
    #HOS initial matrix
    hos = zeros(nt,nc);
    # for each channel
   for c=1:nc
        ## Recurrent HOS Algorithm for P-Wave detection
          hos[:,c]=rec_hos2(dDAS.data[:,c]; C=C, order=order, var_min=var_min, definition=definition) 
          printstyled(" Rec HOS of channel $c of $nc.\n",color=:yellow)       
    end
    
    #Clear the matrix for NaN, 
    hos[isnan.(hos)] .= 0.0; #replace any NAN with zero. Delete for speed
    #stack the rec hos results
    Σhos=sum(hos,dims=2)[:]./nc
    #Σhos=Σhos./maximum(Σhos); # Yes we will normalize the HOZ because is hard to predic the output value!
   
    #find the events 
    ipeaks = findpeaks(Σhos, dDAS.time, min_height=Float64(maximum(Σhos)*thr), min_dist=0.001)
    event_count=size(ipeaks,1);
    peaks = [dDAS.time[ipeaks] dDAS.htime[ipeaks] ipeaks]
    #stack trhe strainrate too
    Σstrainrate = sum(abs.(dDAS.data), dims=2)
    Σstrainrate = Σstrainrate[:] #because sometimes computational languages are dum!
    
     if event_count > 0
         printstyled(" HOS detected $event_count events.\n", color=:green)
         #Make a plot
    
         fig = Figure()
         hm1 = iheatmap(fig[1:2, 1], dDAS.time, dDAS.offset, hos, colormap=:PuRd, colorrange=(0, round(maximum(hos) * 0.5)))
         hm1.axis.ylabel = "Offset [m]"
         Colorbar(fig[1:2, 2], hm1.plot, label="Rec HOS")
         hm1.axis.title = "Rec HOS.\n"
         p2 = ilines(fig[3, 1], dDAS.time, Σstrainrate, color=:grey)
         scatter!(fig[3, 1], dDAS.time[ipeaks], Σstrainrate[ipeaks], marker=:dtriangle, markersize=10, color=:black)
         p2.axis.ylabel = "ηε"
         p3 = ilines(fig[4, 1], dDAS.time, Σhos, color=:magenta)
         scatter!(fig[4, 1], dDAS.time[ipeaks], Σhos[ipeaks], marker=:dtriangle, markersize=10, color=:black)
         p3.axis.xlabel = "Relative Time [s]"
         p3.axis.ylabel = "Mean Rec HOS"
         linkxaxes!(hm1.axis, p2.axis, p3.axis)
         DataInspector(fig)
         display(fig)
    
     else
         printstyled(" HOS detected no events.\n", color=:red)
    
    
         fig = Figure()
         hm1 = iheatmap(fig[1:2, 1], dDAS.time, dDAS.offset, hos, colormap=:PuRd	, colorrange=(0, round(maximum(hos) * 0.5)))
         hm1.axis.ylabel = "Offset [m]"
         Colorbar(fig[1:2, 2], hm1.plot, label="Rec HOS")
         hm1.axis.title = "Rec HOS.\n"
         p2 = ilines(fig[3, 1], dDAS.time, Σstrainrate, color=:grey)
         p2.axis.ylabel = "ηε"
         p3 = ilines(fig[4, 1], dDAS.time, Σhos, color=:magenta)
         p3.axis.xlabel = "Relative Time [s]"
         p3.axis.ylabel = "Mean Rec HOS"
         linkxaxes!(hm1.axis, p2.axis, p3.axis)
         DataInspector(fig)
         display(fig)
    
    
     end
    
    
          return hos, Σhos, peaks, fig
    end
    



##########################################################################
# apply rec_hos
##########################################################################


function channelrechos(strainrate,time,offset; x, C::Float64, order::Int64, var_min=-1, definition=0, thr=0.8) 
    #Select the channel and make sure the usern knows. 
    channel=argmin(abs.(offset.-x));
    schannel=offset[channel];
    printstyled(" Rec HOS of $channel at distance $schannel m.\n",color=:yellow)
## Recurrent HOS Algorithm for P-Wave detection
    signal = strainrate[channel,:]    
    hos=rec_hos2(signal; C=C, order=order, var_min=var_min, definition=definition)
    
    peaks = findpeaks(hos, time, min_height=Float64(thr), min_dist=0.05)
    event_count=size(peaks,1);
     
    
      if event_count > 0
          printstyled(" REC HOS detected $event_count events.\n",color=:green)
          #Make a plot
          # p1=Plots.plot(time,signal,lc=:black, lw=0.5, label="StrainRate",legend=:topright)
          # Plots.plot!(time[peaks], seriestype="vline", lw=0.5,lc=:red, label="Events?")
          # p2=Plots.plot(time,hos, lw=1, label="REC HOS")
          # Plots.scatter!(time[peaks],hos[peaks],mc=:black,markershape=:dtriangle,label="Events?" )
          # p = Plots.plot(p1,p2,layout=(2,1))
          # Plots.savefig("REC_HOS_Channel.pdf")  
      else
          printstyled(" REC HOS detected no events.\n",color=:red)
      end
      return peaks, hos
end

##




######################################################################
# FUNCTIONS FOR DETECTION IN CHANNELS
######################################################################


"""
channelstalta: Run Sta/Lta in one channel

Input:
 - dDAS: DAS data structure
 - x: distance of the channel
 - stw: short time window lenght
 - ltw: long time window lenght
 - thr: threshold
Outputs:
- sra: Sta/Lta function
- peaks: time and index of the Sta/Lta detections
- fig: A figure showing the detected events

Notes:
# Example: Run Sta/Lta in channel x=1600.
```
julia> sra, peaks, fig = channelstalta(dDAS; x=1600, stw=0.002, ltw=0.05, thr=4)

```
"""
function channelstalta(dDAS; x, stw=0.002, ltw=0.05, thr=4)
    #Select the channel and make sure the usern knows.
    channel = argmin(abs.(dDAS.offset .- x))
    schannel = dDAS.offset[channel]
    printstyled(" STA/LTA of $channel at distance $schannel m.\n", color=32)

    # strainrate to data variable
    data = dDAS.data[:, channel]
    dt = dDAS.time[2] - dDAS.time[1] #dt

    ns = Int64(floor(stw / dt)) #no. of data points in the short time window
    nl = Int64(floor(ltw / dt)) #no. of data points in the long time window
    nt = size(data, 1) #safe size for later
    sta = zeros(nt) #initiate variable
    lta = zeros(nt) #initiate variable
    sra = zeros(nt) #initiate variable



    # for every time window
    #Up wind 0 to nl
    for k = 1:nl
        sta[k] = (1 / ns) * sum(abs.(data[k:k+ns]))
        #sta[k] = (1/ns)* trapz.(abs.(data[k-ns:k]));
        lta[k] = (1 / nl) * sum(abs.(data[k:k+nl]))
        #lta[k] = (1/nl)* trapz.(abs.(data[k-ns:k]));
        sra[k] = sta[k] ./ lta[k]
    end
    #from nl onwards
    for k = nl+1:nt

        sta[k] = (1 / ns) * sum(abs.(data[k-ns:k]))
        #sta[k] = (1/ns)* trapz.(abs.(data[k-ns:k]));
        lta[k] = (1 / nl) * sum(abs.(data[k-nl:k]))
        #lta[k] = (1/nl)* trapz.(abs.(data[k-ns:k]));
        sra[k] = sta[k] ./ lta[k]
    end

    #sra[1:nl].=mean(sra[nl+1:end])
    #find all points higher than threshold and print message to user
    ipeaks = findpeaks(sra, dDAS.time, min_height=Float64(thr), min_dist=0.001)
    event_count = size(ipeaks, 1)
    peaks = [dDAS.time[ipeaks] ipeaks]

    if event_count > 0
        printstyled(" STA/LTA detected $event_count events.\n", color=:green)
        #Make a plot
        fig = Figure()
        p1 = ilines(fig[1, 1], dDAS.time, data)
        scatter!(fig[1, 1], dDAS.time[ipeaks], data[ipeaks], marker=:dtriangle, markersize=20, color=:black)
        p1.axis.ylabel = "ηε"
        p1.axis.title = "Sta/Lta of $channel at distance $schannel m.\n"

        p2 = ilines(fig[2, 1], dDAS.time, sra)
        scatter!(fig[2, 1], dDAS.time[ipeaks], sra[ipeaks], marker=:dtriangle, markersize=20, color=:black)
        p2.axis.xlabel = "Relative Time [s]"
        p2.axis.ylabel = "Sta/Lta"

        linkxaxes!(p1.axis, p2.axis)
        display(fig)

    else
        printstyled(" STA/LTA detected no events.\n", color=:red)

        #Make a plot
        fig = Figure()
        p1 = ilines(fig[1, 1], dDAS.time, data)
        p1.axis.ylabel = "ηε"
        p1.axis.title = "Sta/Lta of $channel at distance $schannel m.\n"

        p2 = ilines(fig[2, 1], dDAS.time, sra)
        p2.axis.xlabel = "Relative Time [s]"
        p2.axis.ylabel = "Sta/Lta"

        linkxaxes!(p1.axis, p2.axis)
        display(fig)

    end
    return sra, peaks, fig
end





function channelmer(strainrate, time, offset; x, L=0.01, thr=0.8)
    #Select the channel and make sure the usern knows.
    channel = argmin(abs.(offset .- x))
    schannel = offset[channel]
    printstyled(" MER of $channel at distance $schannel m.\n", color=:yellow)
    ## STA-LTA Algorithm gor P-Wave detection
    grm = strainrate[channel, :]
    dt = time[2] - time[1]
    nt = size(grm, 1) #safe size for later
    nl = Int64(floor(L / dt)) #no. of data points in the long time window

    #Energy ratio. Check formula
    er = zeros(nt)
    for k = nl+1:nt-nl
        #er[k]=sum(grm[k:k+nl] .^2)./ sum(grm[k-nl:k] )
        er[k] = sum(grm[k:k+nl] .^ 2) ./ sum(grm[k-nl:k] .^ 2)
    end

    # modified energy ratio (Crewes). We Normalize to 1.... fo
    mer = er .* abs.(grm) .^ 3
    mer = mer ./ maximum(mer)

    #sra[1:nl].=mean(sra[nl+1:end])
    #find all points higher than threshold and print message to user
    peaks = findpeaks(mer, time, min_height=Float64(thr), min_dist=0.001)
    event_count = size(peaks, 1)

    if event_count > 0
        printstyled(" MER detected $event_count events.\n", color=:green)
        #Make a plot
        p1 = Plots.plot(time, grm, lc=:black, lw=0.5, label="StrainRate")
        Plots.plot!(time[peaks], seriestype="vline", lw=0.5, lc=:red, label="Events?")
        p2 = Plots.plot(time, mer, lw=1, label="MER")
        Plots.scatter!(time[peaks], mer[peaks], mc=:black, markershape=:dtriangle, label="Events?")
        p = Plots.plot(p1, p2, layout=(2, 1))
        Plots.savefig("MER_Channel.pdf")
    else
        printstyled(" MER detected no events.\n", color=:red)
    end
    return peaks, mer
end





######################################################################
# FUNCTIONS FOR RECURRENT HIGH ORDER STATISTICS
######################################################################


"rec_hos.jl"

"""
rec_hos:
NOTE: This is a translation of Claudio Satriano's original rec_hos_py code!

Recursive high order statistics (hos) of a signal, pure Python
implementation.

Defined by Poiata2018 (definition 0) as:

    hos[i] = C·(signal[i]-μ[i-1])ⁿ / (σ²[i])ⁿᐟ² + (1-C)·hos[i-1]

with

    σ²[i] = C·(signal[i]-μ[i-1])² + (1-C)·σ²[i-1]

Or, defined as in Langet2014 (definition 1) as:

    hos[i] = C·(signal[i]-μ[i])ⁿ / (σ²[i])ⁿᐟ² + (1-C)·hos[i-1]

with

    σ²[i] = C·(signal[i]-μ[i])² + (1-C)·σ²[i-1]

For both definitions:

    μ[i] = C·signal[i] + (1-C)·μ[i-1]

Input:
signal : signal to compute recursive hos for
C : decay constant, in the [0, 1] interval
order : hos order
var_min : values of variance σ² (hos denominator) smaller than `var_min` will be replaced by `var_min`
definition : which formula to use

Output:
hos: the recursive hos, with the same length than signal
"""

using StatsBase

function rec_hos(signal; C::Float64, order::Int64, var_min= 1e-9, definition=0)

ss=size(signal,1) #signal size
mean = zeros(ss)
var = ones(ss)
hos = zeros(ss)
n_win=Int64(round(1/C))
# initialize:
@inbounds @fastmath for i=1:n_win
    mean[end] = C * signal[i] + (1 - C) * mean[end]
    var[end] = C * (signal[i] - mean[end])^2 + (1 - C) * var[end]
if definition == 0
    for j = 2:ss
        mean[j] = C * signal[j] + (1 - C) * mean[j-1]
        var[j] = C * (signal[j] - mean[j-1])^2 + (1 - C) * var[j-1]
        _var = min(var[j],var_min)       
        hos[j] = C * ((signal[j] - mean[j-1])^order / _var^(order*0.5)) + (1 - C) * hos[j-1]
       
    end

elseif definition == 1
    for j = 2:ss
        mean[j] = C * signal[j] + (1 - C) * mean[j-1]
        var[j] = C * (signal[j] - mean[j])^2 + (1 - C) * var[j-1]
        _var = min(var[j],var_min)
        hos[j] = C * ((signal[j] - mean[i])^order / _var^(order*0.5)) +(1 - C) * hos[j-1]
    end
end
end

return hos
end

######################## V2


function rec_hos2(signal; C::Float64, order::Int64, var_min= 1e-9, definition=0)

    ss=size(signal,1) #signal size
    mean = zeros(ss)
    var = ones(ss)
    hos = zeros(ss)
    n_win=Int64(round(1/C))
    
    if definition == 0
        @inbounds @fastmath for i=1:n_win
            mean[end] = C * signal[i] + (1 - C) * mean[end]
            var[end] = C * (signal[i] - mean[end])^2 + (1 - C) * var[end]
                for j = 2:ss
                    mean[j] = C * signal[j] + (1 - C) * mean[j-1]
                    var[j] = C * (signal[j] - mean[j-1])^2 + (1 - C) * var[j-1]
                    _var = min(var[j],var_min)       
                    hos[j] = C * ((signal[j] - mean[j-1])^order / _var^(order*0.5)) + (1 - C) * hos[j-1]  
                end
            end
    end

     if definition == 1
        @inbounds @fastmath for i=1:n_win
            mean[end] = C * signal[i] + (1 - C) * mean[end]
            var[end] = C * (signal[i] - mean[end])^2 + (1 - C) * var[end]
            for j = 2:ss
                mean[j] = C * signal[j] + (1 - C) * mean[j-1]
                var[j] = C * (signal[j] - mean[j])^2 + (1 - C) * var[j-1]
                _var = min(var[j],var_min)
                hos[j] = C * ((signal[j] - mean[i])^order / _var^(order*0.5)) +(1 - C) * hos[j-1]
            end
        end
    end
    
    return hos
end

export staltadas, dasmer, dasrechos, channelrechos, sra:, channelstalta, channelmer, rec_hos, rec_hos2
    
