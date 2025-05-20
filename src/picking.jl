"Pickers for dDAS and DAVader!"

# Using
using StatsBase:mean 
using DelimitedFiles
using DSP:conv
using StatsBase:mode

"""
filterpicker: This is a reinterpretation of the well known functions in Matlab and Python of the FiltPicker code. By Mariano Arnaiz (mararnai@ipgp.fr). Mars 2024. 

FilterPicker (Lomax et al., 2011) is a general purpose, broad-band, phase detector and picker which is applicable to real-time seismic monitoring and earthquake early-warning. FilterPicker uses an efficient algorithm which operates stably on continuous, real-time, broadband signals, avoids excessive picking during large events, and produces onset timing, realistic timing uncertainty, onset polarity and amplitude information. In practice, it operates on a pre-defined number of frequency bands by generating a set of band-passed time-series with different center frequencies. Characteristic functions are determined for each frequency band and a pick is declared if and when, within a window of predefined time width, the integral of the maximum of the characteristic functions exceeds a pre-defined threshold.

References:

Lomax, A., C. Satriano and M. Vassallo (2012), Automatic picker developments and optimization: FilterPicker - a robust, broadband picker for real-time seismic monitoring and earthquake early-warning, Seism. Res. Lett. , 83, 531-540, doi: 10.1785/gssrl.83.3.531.

Vassallo, M., C. Satriano and A. Lomax, (2012), Automatic picker developments and optimization: A strategy for improving the performances of automatic phase pickers, Seism. Res. Lett. , 83, 541-554, doi: 10.1785/gssrl.83.3.541.

Inputs: 
 - signal: Event data. We will try to pick the P wave and other arrival on this signal.
 - sr: Sampling rate (in Hz) of signal.
 - tflt: in samples, determines how far back in time the previous samples are examined. The filter window will be adjusted upwards to be an integer N power of 2 times the sample interval (deltaTime).  Then numRecursive = N + 1 "filter bands" are created.  For each filter band n = 0,N  the data samples are processed through a simple recursive filter backwards from the current sample, and picking statistics and characteristic function are generated.  Picks are generated based on the maximum of the characteristic funciton values over all filter bands relative to the threshold values threshold1 and threshold2.(default=300)
 - tlong: the long term window (longTermWindow) determines: a) a stabilisation delay time after the beginning of data; before this delay time picks will not be generated. b) the decay constant of a simple recursive filter to accumlate/smooth all picking statistics and characteristic functions for all filter bands.(default=500)
 - tup: determines the maximum time the integral of the (clipped) characteristic function is accumlated after threshold1 is reached (pick event triggered) to check for this integral exceeding threshold2 * tUpEvent (pick declared). (default=20)
 - thr1: sets the threshold to trigger a pick event (potential pick).  This threshold is reached when the (clipped) characteristic function for any filter band exceeds threshold1.(default=10)
 - thr2: sets the threshold to declare a pick (pick will be accepted when tUpEvent reached).  This threshold is reached when the integral of the (clipped) characteristic function for any filter band over the window tUpEvent exceeds threshold2 * tUpEvent (i.e. the average (clipped) characteristic function over tUpEvent is greater than threshold2). (default=10)

Outputs: All te outputs are written to a singles structure array (e.g. picks) that holds all the outputs. To call uppon them use:
 - picks.idx: index (sample) of the probable picks
 - picks.unc: picks uncertanty
 - picks.f: index of the frequency used for each pick
 - picks.t: time [s] after 0.0 (begining of signal) for each pick
 - picks.signal_lp: low passed filtered signal on each frequency band
 - picks.Fn: characteristic function at each frequency band
 - picks.FnS: characteristic function with the max of each frequency band
 - picks.FnMAvg: convolved FnS with the clip set by tup

Notes:
 - I have written most of the resulst to output for later use.
 - No plot is created by the main code.
 - No processing is done to the signal! If you think this is necessary do it before.
 - A fast convolution algorithm is used. See functions at the end of the file for details. 

# Example 1: Run over a signal with sr=48000 with default parameters.
```
julia> picks=filterpicker(signal,48000)
```
# Example 2: Run over a signal with sr=48000 with different parameters.
```
julia> picks=filterpicker(signal,48000; tflt=200, tlong=400, tup=10, thr1=8, thr2=8)
```
"""
function filterpicker(signal,sr; tflt=300, tlong=500, tup=20, thr1=10, thr2=10)

#tflt=300; tlong=500; tup=20; thr1=6; thr2=6; sr=48000


dt=1/sr # dt from sampling rate
numSMP=size(signal,1) #number of samples

#Initial parameters
tflt = 2.0^ceil(log(tflt)/log(2.0))   #Filter window in dt
dc_lng = 1-(1/tlong)  #Decay long constant
tup = max(1,tup)    #tUpEvent in dt
min_sig = eps(0.0) # min value of signal
Navg = min(round(tlong),numSMP) 
y0   = mean(signal[1:Navg])  # set y(0) to mean of y(i) over Tlong
numBnd  = Int64(ceil(log(tflt)/log(2.0)))
#Iniriate variables
Fn  = zeros(numBnd,numSMP) 
FnL = zeros(numBnd,numSMP) 
signal_lp = zeros(numBnd,numSMP) 
## Loop over bands
@inbounds for n=1:numBnd 
    w   = (2.0^n*dt)/(2.0*pi)   #Band omega
    cHP = w/(w+dt)          #High pass filter coefficient
    cLP = dt/(w+dt)         #Low pass filter coefficient
    yHP1p = 0.0
    yHP2p = 0.0
    yLPp  = 0.0
    avg  = 0.0
    vrn  = 0.0
    sig  = 0.0
    tmpFNL=thr1/2
    
   @inbounds for i=1:numSMP
        if i != 1  
            yHP1 = cHP*(yHP1p + signal[i]-signal[i-1])    #First high-pass
        else
            yHP1 = cHP*(yHP1p + signal[i]-y0)    #First high-pass
        end
        
        yHP2 = cHP*(yHP2p + yHP1 - yHP1p)   #Second high-pass
        yLP  = yLPp + cLP*(yHP2 - yLPp)     #Low-pass
        En   = yLP^2                        #Envelope
        yHP1p = yHP1
        yHP2p = yHP2
        yLPp  = yLP
        signal_lp[n,i] = yLP #only for plots
        if sig > min_sig
            if (En-avg)/sig > 5*thr1
                En = 5*thr1*sig + avg
                Fn[n,i] = 5*thr1
            else
                Fn[n,i] = (En-avg)/sig  #Characteristic function
                if Fn[n,i] < 1
                    Fn[n,i] = 0.0
                end
            end
        end
        
        avg = dc_lng*avg + (1-dc_lng)*En     #Mean
        vrn = dc_lng*vrn + (1-dc_lng)*(En-avg)^2  #Variance
        sig = sqrt(vrn) #Standard deviation
        
        tmpFNL = dc_lng*tmpFNL + (1-dc_lng)*Fn[n,i]
        tmpFNL = min(thr1/2,tmpFNL)
        tmpFNL = max(0.5,tmpFNL)
        FnL[n,i]=tmpFNL
    end
end

FnS, iFnS =findmax(Fn,dims=1) #max(Fn,[],1)     #Maximum of all bands
FnS[ FnS .<0.0 ] .= 0.0
FnMAvg = fastconv(FnS,[0 ones(1,tup-1) 0]/tup)
PotTrg = findall(FnS .> thr1) #Vector of potential triggers

## Make an automatic pick    
ff1=findall(PotTrg .< CartesianIndex(1,(tlong+1)) .|| PotTrg .> CartesianIndex(1,numSMP-tup))
deleteat!(PotTrg,ff1)

flagTrg= ones(Int64,size(PotTrg))
pickIDX= zeros(size(PotTrg)) #time of the pick
pickUNC= zeros(size(PotTrg)) #uncertainty ...
pickFRQ= zeros(size(PotTrg)) #frequency ...

@inbounds for i=1:length(flagTrg)
        tmpTrg  = PotTrg[i]
        if FnMAvg[tmpTrg] > thr2
            bandtrg = findfirst(Fn[:,tmpTrg[2]].>thr1)
            firstPot = findfirst(Fn[bandtrg,tmpTrg[2]:-1:1]-FnL[bandtrg,tmpTrg[2]:-1:1].<0)
            if(isempty(firstPot))
                firstPot=0.0
            end
            pickFRQ[i]  = bandtrg
            pickIDX[i]  = tmpTrg[2] - firstPot
            limT        = (2.0^(bandtrg))/20 #(2.0^(bandtrg-1))/20
            if firstPot < limT 
                pickUNC[i] = limT*dt
            else
                pickUNC[i] = firstPot*dt
            end
            #Disable next potential picks until FnS drops below 2
            idxDrop  = findfirst(FnS[tmpTrg[2]:end] .< 2)
            
            if isempty(idxDrop) == false
            ff2=findall(PotTrg .> tmpTrg .&& PotTrg .<= CartesianIndex(1,tmpTrg[2] + idxDrop))
            flagTrg[ff2].=0
            end
        end
end

#Clean the vectors of zeros
#find the index of 0.0
ffpidx=findall(pickIDX .> 0.0)
ffpunc=findall(pickUNC .> 0.0)
fffrq=findall(pickFRQ .> 0.0)
#Remove zero index
PotTrg =PotTrg[ffpidx] 
pickIDX=pickIDX[ffpidx] 
pickUNC=pickUNC[ffpunc]
pickFRQ=pickFRQ[fffrq]
pickTIME=pickIDX*dt

#search for pics :): 2 methods
#method 1: max uncertanty
if isempty(pickUNC) == false
   pick_indx1=Int64((pickIDX[argmax(pickUNC)]))
   pick_time1=pick_indx1*dt
else
   pick_indx1=1
   pick_time1=0.0
   printstyled(" No pick found by max uncertanty method \n",color=:red)
end

#method 2: most repeated sample
if isempty(pickIDX) == false
   pick_indx2=Int64(Int64(mode(pickIDX)))
   pick_time2=pick_indx2*dt
else
   pick_indx2=1
   pick_time2=0.0
   printstyled(" No pick found by most common sample method \n",color=:red)
end

#Organize all the info:
results=[pick_indx1 pick_time1], [pick_indx2 pick_time2],fpstruc(pickIDX,pickUNC,pickFRQ, pickTIME, signal_lp, Fn,  FnS,  FnMAvg)

    return results#pick_time,pick_indx
end




########################################################

"Filter Picker for Julia with Ali's modification"

# Using
using StatsBase:mean 
using DelimitedFiles, GLMakie
using DSP:conv
using StatsBase:mode

"""
filterpicker_ali: This is a reinterpretation of the well known functions in Matlab and Python of the FiltPicker code. By Mariano Arnaiz (mararnai@ipgp.fr). Mars 2024. 

FilterPicker (Lomax et al., 2011) is a general purpose, broad-band, phase detector and picker which is applicable to real-time seismic monitoring and earthquake early-warning. FilterPicker uses an efficient algorithm which operates stably on continuous, real-time, broadband signals, avoids excessive picking during large events, and produces onset timing, realistic timing uncertainty, onset polarity and amplitude information. In practice, it operates on a pre-defined number of frequency bands by generating a set of band-passed time-series with different center frequencies. Characteristic functions are determined for each frequency band and a pick is declared if and when, within a window of predefined time width, the integral of the maximum of the characteristic functions exceeds a pre-defined threshold.

References:

Lomax, A., C. Satriano and M. Vassallo (2012), Automatic picker developments and optimization: FilterPicker - a robust, broadband picker for real-time seismic monitoring and earthquake early-warning, Seism. Res. Lett. , 83, 531-540, doi: 10.1785/gssrl.83.3.531.

Vassallo, M., C. Satriano and A. Lomax, (2012), Automatic picker developments and optimization: A strategy for improving the performances of automatic phase pickers, Seism. Res. Lett. , 83, 541-554, doi: 10.1785/gssrl.83.3.541.

Inputs: 
 - signal: Event data. We will try to pick the P wave and other arrival on this signal.
 - sr: Sampling rate (in Hz) of signal.
 - tflt: in samples, determines how far back in time the previous samples are examined. The filter window will be adjusted upwards to be an integer N power of 2 times the sample interval (deltaTime).  Then numRecursive = N + 1 "filter bands" are created.  For each filter band n = 0,N  the data samples are processed through a simple recursive filter backwards from the current sample, and picking statistics and characteristic function are generated.  Picks are generated based on the maximum of the characteristic funciton values over all filter bands relative to the threshold values threshold1 and threshold2.(default=300)
 - tlong: the long term window (longTermWindow) determines: a) a stabilisation delay time after the beginning of data; before this delay time picks will not be generated. b) the decay constant of a simple recursive filter to accumlate/smooth all picking statistics and characteristic functions for all filter bands.(default=500)
 - tup: determines the maximum time the integral of the (clipped) characteristic function is accumlated after threshold1 is reached (pick event triggered) to check for this integral exceeding threshold2 * tUpEvent (pick declared). (default=20)
 - thr1: sets the threshold to trigger a pick event (potential pick).  This threshold is reached when the (clipped) characteristic function for any filter band exceeds threshold1.(default=10)
 - thr2: sets the threshold to declare a pick (pick will be accepted when tUpEvent reached).  This threshold is reached when the integral of the (clipped) characteristic function for any filter band over the window tUpEvent exceeds threshold2 * tUpEvent (i.e. the average (clipped) characteristic function over tUpEvent is greater than threshold2). (default=10)

Outputs: All te outputs are written to a singles structure array (e.g. picks) that holds all the outputs. To call uppon them use:
 - picks.idx: index (sample) of the probable picks
 - picks.unc: picks uncertanty
 - picks.f: index of the frequency used for each pick
 - picks.t: time [s] after 0.0 (begining of signal) for each pick
 - picks.signal_lp: low passed filtered signal on each frequency band
 - picks.Fn: characteristic function at each frequency band
 - picks.FnS: characteristic function with the max of each frequency band
 - picks.FnMAvg: convolved FnS with the clip set by tup

Notes:
 - I have written most of the resulst to output for later use.
 - No plot is created by the main code.
 - No processing is done to the signal! If you think this is necessary do it before.
 - A fast convolution algorithm is used. See functions at the end of the file for details. 

# Example 1: Run over a signal with sr=48000 with default parameters.
```
julia> picks=filterpicker_ali(signal,48000)
```
# Example 2: Run over a signal with sr=48000 with different parameters.
```
julia> picks=filterpicker_ali(signal,48000; tflt=200, tlong=400, tup=10, thr1=8, thr2=8)
```
"""
function filterpicker_ali(signal,sr; tflt=300, tlong=500, tup=20, thr1=10, thr2=10)

#tflt=300; tlong=500; tup=20; thr1=6; thr2=6; sr=48000


dt=1/sr # dt from sampling rate
numSMP=size(signal,1) #number of samples

#Initial parameters
tflt = 2.0^ceil(log(tflt)/log(2.0))   #Filter window in dt
dc_lng = 1-(1/tlong)  #Decay long constant
tup = max(1,tup)    #tUpEvent in dt
min_sig = eps(0.0) # min value of signal
Navg = min(round(tlong),numSMP) 
y0   = mean(signal[1:Navg])  # set y(0) to mean of y(i) over Tlong
#numBnd  = Int64(ceil(log(tflt)/log(2.0)))+1
numBnd  = Int64(ceil(log2(tflt/(dt))+1))

#Iniriate variables
Fn  = zeros(numBnd,numSMP) 
FnL = zeros(numBnd,numSMP) 
signal_lp = zeros(numBnd,numSMP) 

## Loop over bands
@inbounds for n=1:numBnd 
    w   = (2.0^n*dt)/(2.0*pi)   #Band omega
    cHP = w/(w+dt)          #High pass filter coefficient
    cLP = dt/(w+dt)         #Low pass filter coefficient
    yHP1p = 0.0
    yHP2p = 0.0
    yLPp  = 0.0
    avg  = 0.0
    vrn  = 0.0
    sig  = 0.0
    dy=0.0 # Ali´s filter
    dyp = 0.0 # Ali´s filter
    yHP3 = 0.0 # Ali´s filter
    yHP3p = 0.0 # Ali´s filter
    yHP4 = 0.0 # Ali´s filter
    yHP4p = 0.0 # Ali´s filter
    yLP_new_p = 0.0 # Ali´s filter

    tmpFNL=thr1/2
    
   @inbounds for i=1:numSMP
        if i != 1  
            yHP1 = cHP*(yHP1p + signal[i]-signal[i-1])    #First high-pass
        else
            yHP1 = cHP*(yHP1p + signal[i]-y0)    #First high-pass
        end
        yHP2 = cHP*(yHP2p + yHP1 - yHP1p)   #Second high-pass
        yLP  = yLPp + cLP*(yHP2 - yLPp)     #Low-pass
        dy=yLP-yLPp  # Ali´s filter
        yHP3= cHP*(yHP3p + dy - dyp) # Ali´s filter
        yHP4= cHP*(yHP4p + yHP1 - yHP3p) # Ali´s filter
        yLP_new  = yLP_new_p + cLP*(yHP2 - yHP4p)     #Ali´s Low-pass

        #En   = yLP^2  #Envelope
        En   = yLP_new^2  #Envelope

        yHP1p = yHP1
        yHP2p = yHP2
        yLPp  = yLP
        dyp = dy # Ali´s filter
        yHP3p = yHP3 # Ali´s filter
        yHP4p = yHP4 # Ali´s filter
        yLPp_new  = yLP_new # Ali´s filter

        signal_lp[n,i] = yLP #only for plots
        if sig > min_sig
            if (En-avg)/sig > 5*thr1
                En = 5*thr1*sig + avg
                Fn[n,i] = 5*thr1
            else
                Fn[n,i] = (En-avg)/sig  #Characteristic function
                if Fn[n,i] < 1
                    Fn[n,i] = 0.0
                end
            end
        end
        
        avg = dc_lng*avg + (1-dc_lng)*En     #Mean
        vrn = dc_lng*vrn + (1-dc_lng)*(En-avg)^2  #Variance
        sig = sqrt(vrn) #Standard deviation
        
        tmpFNL = dc_lng*tmpFNL + (1-dc_lng)*Fn[n,i]
        tmpFNL = min(thr1/2,tmpFNL)
        tmpFNL = max(0.5,tmpFNL)
        FnL[n,i]=tmpFNL
    end
end

FnS, iFnS =findmax(Fn,dims=1) #max(Fn,[],1)     #Maximum of all bands
FnS[ FnS .<0.0 ] .= 0.0
FnMAvg = fastconv(FnS,[0 ones(1,tup-1) 0]/tup)
PotTrg = findall(FnS .> thr1) #Vector of potential triggers

## Make an automatic pick    
ff1=findall(PotTrg .< CartesianIndex(1,(tlong+1)) .|| PotTrg .> CartesianIndex(1,numSMP-tup))
deleteat!(PotTrg,ff1)

flagTrg= ones(Int64,size(PotTrg))
pickIDX= zeros(size(PotTrg)) #time of the pick
pickUNC= zeros(size(PotTrg)) #uncertainty ...
pickFRQ= zeros(size(PotTrg)) #frequency ...

@inbounds for i=1:length(flagTrg)
        tmpTrg  = PotTrg[i]
        if FnMAvg[tmpTrg] > thr2
            bandtrg = findfirst(Fn[:,tmpTrg[2]].>thr1)
            firstPot = findfirst(Fn[bandtrg,tmpTrg[2]:-1:1]-FnL[bandtrg,tmpTrg[2]:-1:1].<0)
            if(isempty(firstPot))
                firstPot=0.0
            end
            pickFRQ[i]  = bandtrg
            pickIDX[i]  = tmpTrg[2] - firstPot
            limT        = (2.0^(bandtrg))/20 #(2.0^(bandtrg-1))/20
            if firstPot < limT 
                pickUNC[i] = limT*dt
            else
                pickUNC[i] = firstPot*dt
            end
            #Disable next potential picks until FnS drops below 2
            idxDrop  = findfirst(FnS[tmpTrg[2]:end] .< 2)
            
            if isempty(idxDrop) == false
            ff2=findall(PotTrg .> tmpTrg .&& PotTrg .<= CartesianIndex(1,tmpTrg[2] + idxDrop))
            flagTrg[ff2].=0
            end
        end
end

#Clean the vectors of zeros
#find the index of 0.0
ffpidx=findall(pickIDX .> 0.0)
ffpunc=findall(pickUNC .> 0.0)
fffrq=findall(pickFRQ .> 0.0)
#Remove zero index
PotTrg =PotTrg[ffpidx] 
pickIDX=pickIDX[ffpidx] 
pickUNC=pickUNC[ffpunc]
pickFRQ=pickFRQ[fffrq]
pickTIME=pickIDX*dt

#search for pics :): 2 methods
#method 1: max uncertanty
if isempty(pickUNC) == false
   pick_indx1=min(Int64((pickIDX[argmax(pickUNC)])),Int64(mode(pickIDX)))
   pick_time1=pick_indx1*dt
else
   pick_indx1=1
   pick_time1=0.0
   printstyled(" No pick found by max uncertanty vs mode method \n",color=:red)
end

#method 2: most repeated sample
if isempty(pickIDX) == false
   #pick_indx2=Int64(Int64(mode(pickIDX)))
   pick_indx2=Int64(minimum(pickIDX))
   pick_time2=pick_indx2*dt
else
   pick_indx2=1
   pick_time2=0.0
   printstyled(" No pick found by first sample \n",color=:red)
end

#Organize all the info:
results=[pick_indx1 pick_time1], [pick_indx2 pick_time2],fpstruc(pickIDX,pickUNC,pickFRQ, pickTIME, signal_lp, Fn,  FnS,  FnMAvg)

    return results#pick_time,pick_indx
end



"""
filtpickdas: This is a function to run filterpicker (read ? filterpicker for more info) in channles of a dDAs structure. 
Inputs: 
 - dDAS: iDAS structure with DAS data loaded.
 - x: distance/offset range to pick over.
 - tflt: in samples, determines how far back in time the previous samples are examined. The filter window will be adjusted upwards to be an integer N power of 2 times the sample interval (deltaTime).  Then numRecursive = N + 1 "filter bands" are created.  For each filter band n = 0,N  the data samples are processed through a simple recursive filter backwards from the current sample, and picking statistics and characteristic function are generated.  Picks are generated based on the maximum of the characteristic funciton values over all filter bands relative to the threshold values threshold1 and threshold2.(default=300)
 - tlong: the long term window (longTermWindow) determines: a) a stabilisation delay time after the beginning of data; before this delay time picks will not be generated. b) the decay constant of a simple recursive filter to accumlate/smooth all picking statistics and characteristic functions for all filter bands.(default=500)
 - tup: determines the maximum time the integral of the (clipped) characteristic function is accumlated after threshold1 is reached (pick event triggered) to check for this integral exceeding threshold2 * tUpEvent (pick declared). (default=20)
 - thr1: sets the threshold to trigger a pick event (potential pick).  This threshold is reached when the (clipped) characteristic function for any filter band exceeds threshold1.(default=10)
 - thr2: sets the threshold to declare a pick (pick will be accepted when tUpEvent reached).  This threshold is reached when the integral of the (clipped) characteristic function for any filter band over the window tUpEvent exceeds threshold2 * tUpEvent (i.e. the average (clipped) characteristic function over tUpEvent is greater than threshold2). (default=10)

Outputs: All te outputs are written to a singles structure array (e.g. picks) that holds all the outputs. To call uppon them use:
 - picks table: a table with the picks information (Channel_index Channel_offset Pick_index Pick_time)
Notes:
 - Only the clear picks are written to outpu.
 - No plot is created by the main code.
 - No processing is done to the signal! If you think this is necessary do it before.
 - filterpicker works better if you keep enough data before the arrival!

# Example 1: Run over one channel of the dDAS at x=1300 m with default parameters.
```
julia> picks=filterpicker(dDAS,1300)
```
# Example 2: Run over several channels of the dDAS at x from 1300 to 1360 m.
```
julia> picks=filtpickdas(dDAS,1300:1360; tflt=200, tlong=400, tup=10, thr1=8, thr2=8)
```
"""
function filtpickdas(dDAS, x; tflt=200, tlong=400, tup=10, thr1=8, thr2=8)

    fs= 1/(dDAS.time[2]-dDAS.time[1])
    pocks=zeros(size(x, 1),4)
    
    #if we need to pick in one channel
    if typeof(x) == Int64 || typeof(x) == Float64
        channel = argmin(abs.(dDAS.offset .- x))
        pocks[1]=channel
        schannel = dDAS.offset[channel]
        pocks[2]=schannel
        printstyled(" Picking in Channel $channel at distance $schannel m\n", color=:yellow)
        signal = dDAS.data[:, channel]
        picks=filterpicker(signal,fs; tflt=tflt, tlong=tlong, tup=tup, thr1=thr1, thr2=thr2)
        pocks[3:4]=picks[1]
    end

        #if we need to pick in several channel
    if typeof(x) == UnitRange{Int64} || typeof(x) == StepRange{Int64,Int64} || typeof(x) == StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64} || typeof(x) == StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64}

        for i = 1:size(x, 1)
            channel = argmin(abs.(dDAS.offset .- x[i]))
            pocks[i,1]=channel
            schannel = dDAS.offset[channel]
            pocks[i,2]=schannel
            printstyled(" Picking in Channel $channel at distance $schannel m\n", color=:yellow)
            signal = dDAS.data[:, channel]
            picks=filterpicker(signal,fs; tflt=tflt, tlong=tlong, tup=tup, thr1=thr1, thr2=thr2)
            pocks[i,3:4]=picks[1]
        end
    end

    #we will clean pocks of bad values
    pocks= pocks[findall(pocks[:,4] .!= 0.0),:]

    printstyled(" Output is: Channel_index Channel_offset Pick_index Pick_time\n", color=:green)
    return pocks
end









######3
#Fast convolution
using Base.Cartesian

# fast convolution metond direct version (do not check if threshold is satisfied)
@generated function fastconv(E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote

        retsize = [size(E)...] + [size(k)...] .- 1
        retsize = tuple(retsize...)
        ret = zeros(T, retsize)

        convn!(ret,E,k)
        return ret

    end
end


# in place helper operation to speedup memory allocations
@generated function convn!(out::Array{T}, E::Array{T,N}, k::Array{T,N}) where {T,N}
    quote
        @inbounds begin
            @nloops $N x E begin
                @nloops $N i k begin
                    (@nref $N out d->(x_d + i_d - 1)) += (@nref $N E x) * (@nref $N k i)
                end
            end
        end
        return out
    end
end

#############################################################################################
#############################################################################################
# Structure Arrays section
#############################################################################################
#############################################################################################

## Definition of mutalble structure array attb. This is needed to build variable atrib
# with all those attributes!
mutable struct fpstruc
    idx::Vector{Float64}
    unc::Vector{Float64}
    f::Vector{Float64}
    t::Vector{Float64}
    signal_lp::Matrix{Float64}
    Fn::Matrix{Float64}
    FnS::Matrix{Float64}
    FnMAvg::Matrix{Float64}
end
   
export filterpicker, filterpicker_ali, filtpickdas, fastconv, convn!, fpstruc
