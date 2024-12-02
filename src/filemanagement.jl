"This is the File Management (read, write, etc) part of DAS VADER V1.0"

using HDF5, Dates, JLD2, DelimitedFiles, Seis

export rdas, writejdl, das2sac, das2txt, attb, iDAS

## 1. Using: Loading the package to be used

"""
rdasnew: Read a .h5 DAS file to memory and save the information into variables for latter.

Input:
 - filename: Name of the file to read. Should be written in air quotes (as it is a String). Please include extension.

Outputs:
 - dDAS.data: Strain rate matrix read from the file. We pre-clean the information an assamble the data for easy use.
 - dDAS.time: Time vector in s. This is a relative time to the beginning of the file (time[1] = 0.0). It is written in a range type. If you need the full vector use collect(time).
 - dDAS.htime: Human time vector in Date format. This is the GPS time of the samples as reconstructed from the file.
 - dDAS.offset: Distance vector in m. This is a relative distance to the beginning of the cable (offset[1] = 0.0). It is written in a range type. If you need the full vector use collect(offset).
 - dDAS.atrib: Attributes from the .h5 file. Here all the attributes are stored. Call as:
      atrib.AmpliPower, atrib.BlockRate, atrib.Components, atrib.DataDomain, atrib.DerivationTime, atrib.Extent
      atrib.FiberLength, atrib.GaugeLength, atrib.Hostname, atrib.Origin, atrib.Oversampling, atrib.PipelineTracker
      atrib.PulseRateFreq, atrib.PulseWidth, atrib.SamplingRate, atrib.SamplingRes, atrib.Spacing

Notes:
 - No processing is applied to the data.
 - No figure is created.
 - Reading a DAS .h5 file can easily overload your RAM memory. Proceed with caution.
# Example: Read a DAS .h5 file to memory and save the variables to standard names.
```
julia> dDAS = rdasnew("data.h5")
```
"""
function rdasnew(filename::String)

    # Read the file
    printstyled("\nReading $filename \n", color=:light_blue, bold=true)
    try
        global file = h5open(filename, "r") #Read the file
        printstyled("Outputs: (color coded) \n", color=:light_blue)
        printstyled(" data: strain rate (DAS data) with no time overlaps.\n", color=:green, bold=true)
        printstyled(" time: relative time vector (starts at 00:00).\n", color=:white, bold=true)
        printstyled(" htime: correct human time vector for data.\n", color=:white, bold=true)
        printstyled(" offset: offset vector for data.\n", color=:yellow, bold=true)
        printstyled(" atrib: structure array with original file attributes.\n", color=:magenta)

    catch e
        printstyled("\n $filename is NOT an .h5 file \n", color=:red)
    end

    ## File information: Name of the DAS and datasets in the .h5 file

    # 1. Get Basic Information
    name = keys(file)[1] #Get the Name of the File
    datasets = HDF5.get_datasets(file) #Get the "groups" where the all the data is saved

    # 2. Get the raw strain and time data
    #Raw Strain Rate Data. All of it, in all its raw glory!
    raw_data = read(datasets[1]) # raw strain rate matrix
    raw_utime = read(datasets[2]) # raw unix time vector

    #3. Attributes: Read all attributes from the file
    AmpliPower = attrs(file["$name/Source1/Zone1"])["AmpliPower"]
    BlockRate = attrs(file["$name/Source1/Zone1"])["BlockRate"]
    Components = attrs(file["$name/Source1/Zone1"])["Components"]
    DataDomain = attrs(file["$name/Source1/Zone1"])["DataDomain"]
    DerivationTime = attrs(file["$name/Source1/Zone1"])["DerivationTime"]
    Extent = attrs(file["$name/Source1/Zone1"])["Extent"]
    FiberLength = attrs(file["$name/Source1/Zone1"])["FiberLength"]
    GaugeLength = attrs(file["$name/Source1/Zone1"])["GaugeLength"]
    Hostname = attrs(file["$name/Source1/Zone1"])["Hostname"]
    Origin = attrs(file["$name/Source1/Zone1"])["Origin"]
    Oversampling = attrs(file["$name/Source1/Zone1"])["Oversampling"]
    PipelineTracker = attrs(file["$name/Source1/Zone1"])["PipelineTracker"]
    PulseRateFreq = Float64(attrs(file["$name/Source1/Zone1"])["PulseRateFreq"][1]) / 1000
    PulseWidth = attrs(file["$name/Source1/Zone1"])["PulseWidth"]
    SamplingRate = attrs(file["$name/Source1/Zone1"])["SamplingRate"]
    SamplingRes = attrs(file["$name/Source1/Zone1"])["SamplingRes"]
    Spacing = attrs(file["$name/Source1/Zone1"])["Spacing"]
    BlockOverlap = attrs(file["$name/Source1/Zone1"])["BlockOverlap"]

    #Save all the attributes to the structure
    atrib = attb(AmpliPower, BlockRate, Components, DataDomain, DerivationTime, Extent, FiberLength, GaugeLength, Hostname, Origin, Oversampling, PipelineTracker, PulseRateFreq, PulseWidth, SamplingRate, SamplingRes, Spacing, BlockOverlap)

    # 4. Process the data, time and distance information to build the data

    #4.1 Strain Rate Data
    #raw_data2=permutedims(raw_data,(3,2,1));
    overlap = attrs(file["$name/Source1/Zone1"])["BlockOverlap"][1]
    nb_Block = size(raw_data, 3) #Number of blocks
    Block_Time_size = size(raw_data, 2) #Time size
    Block_Space_size = size(raw_data, 1) #space size
    dim = size(raw_data)

    noverlap = Int64(ceil(dim[2] * (overlap / (100)) / 2))#size(raw_data,2) ÷ 4; #CHECK: Is dim3 the right time dim?
    data = raw_data[:, noverlap+1:end-noverlap, :] #CHECK: are we cleaning in the right dim?
    data = reshape(data, size(data, 1), :)
    data = transpose(data)

    #4.2 Offset or Distance Vectors
    ZI_start = Origin[1] + (Extent[1]) * Spacing[1] #m - Start distance of the zone
    ZI_end = Origin[1] + (Extent[2]) * Spacing[1] #m - End distance of the zone
    offset = ZI_start:Spacing[1]:ZI_end #m - Distance vector of the fiber

    #4.3 Relative time vector
    nb_Block = dim[3] #Number of blocks
    Block_Time_size = dim[2] #Time size
    Block_Space_size = dim[1] #Distance size
    Total_time_size = size(data, 1)#nb_Block*Block_Time_size/2; #### HERE WE CHANGED SOMETHING IMPORTANT!!!!!
    time = 0:(Spacing[2]*1e-3):(Total_time_size-1)*(Spacing[2]*1e-3)

    #4.4 Unix and Human time vector
    utime = range(start=raw_utime[1], step=(Spacing[2] * 1e-3), length=size(time)[1])
    htime = unix2datetime.(utime)
    replace!(data, NaN => 0.0)

    #5: All to the same array
    dDAS = iDAS(data, time, htime, offset, atrib, filename)

    return dDAS
end


################other version of the same
function rdasold(filename::String)

    # Read the file
    printstyled("\nReading $filename \n", color=:light_blue, bold=true)
    try
        global file = h5open(filename, "r") #Read the file
        printstyled("Outputs: (color coded) \n", color=:light_blue)
        printstyled(" strainrate: strain rate (DAS data) with no time overlaps.\n", color=:green, bold=true)
        printstyled(" time: relative time vector (starts at 00:00).\n", color=:white, bold=true)
        printstyled(" htime: correct human time vector for data.\n", color=:white, bold=true)
        printstyled(" offset: offset vector for data.\n", color=:yellow, bold=true)
        printstyled(" atrib: structure array with original file attributes.\n", color=:magenta)

    catch e
        printstyled("\n $filename is NOT an .h5 file \n", color=:red)
    end

    ## File information: Name of the DAS and datasets in the .h5 file

    # 1. Get Basic Information
    name = keys(file)[1] #Get the Name of the File
    datasets = HDF5.get_datasets(file) #Get the "groups" where the all the data is saved

    # 2. Get the raw strain and time data
    #Raw Strain Rate Data. All of it, in all its raw glory!
    raw_data = read(datasets[1]) # raw strain rate matrix
    raw_utime = read(datasets[2]) # raw unix time vector

    #3. Attributes: Read all attributes from the file
    AmpliPower = attrs(file["$name/Source1/Zone1"])["AmpliPower"]
    BlockRate = attrs(file["$name/Source1/Zone1"])["BlockRate"]
    Components = attrs(file["$name/Source1/Zone1"])["Components"]
    DataDomain = attrs(file["$name/Source1/Zone1"])["DataDomain"]
    DerivationTime = attrs(file["$name/Source1/Zone1"])["DerivationTime"]
    Extent = attrs(file["$name/Source1/Zone1"])["Extent"]
    FiberLength = attrs(file["$name/Source1/Zone1"])["FiberLength"]
    GaugeLength = attrs(file["$name/Source1/Zone1"])["GaugeLength"]
    Hostname = attrs(file["$name/Source1/Zone1"])["Hostname"]
    Origin = attrs(file["$name/Source1/Zone1"])["Origin"]
    Oversampling = attrs(file["$name/Source1/Zone1"])["Oversampling"]
    PipelineTracker = attrs(file["$name/Source1/Zone1"])["PipelineTracker"]
    PulseRateFreq = Float64(attrs(file["$name/Source1/Zone1"])["PulseRateFreq"][1]) / 1000
    PulseWidth = attrs(file["$name/Source1/Zone1"])["PulseWidth"]
    SamplingRate = attrs(file["$name/Source1/Zone1"])["SamplingRate"]
    SamplingRes = attrs(file["$name/Source1/Zone1"])["SamplingRes"]
    Spacing = attrs(file["$name/Source1/Zone1"])["Spacing"]
    BlockOverlap = Int32.([50])[:]
    #Save all the attributes to the structure
    atrib = attb(AmpliPower, BlockRate, Components, DataDomain, DerivationTime, Extent, FiberLength, GaugeLength, Hostname, Origin, Oversampling, PipelineTracker, PulseRateFreq, PulseWidth, SamplingRate, SamplingRes, Spacing, BlockOverlap)

    # 4. Process the strainrate, time and distance information to build the data

    #4.1 Strain Rate Data
    #raw_data2=permutedims(raw_data,(3,2,1));
    noverlap = size(raw_data, 2) ÷ 4 #CHECK: Is dim3 the right time dim?
    data = raw_data[:, noverlap+1:end-noverlap, :] #CHECK: are we cleaning in the right dim?
    data = reshape(data, size(data, 1), :)
    data = transpose(data)


    #4.2 Offset or Distance Vectors
    ZI_start = Origin[1] + (Extent[1]) * Spacing[1] #m - Start distance of the zone
    ZI_end = Origin[1] + (Extent[2]) * Spacing[1] #m - End distance of the zone
    offset = ZI_start:Spacing[1]:ZI_end #m - Distance vector of the fiber

    #4.3 Relative time vector
    dim = size(raw_data)
    nb_Block = dim[3] #Number of blocks
    Block_Time_size = dim[2] #Time size
    Block_Space_size = dim[1] #Distance size
    Total_time_size = nb_Block * Block_Time_size / 2
    time = 0:(Spacing[2]*1e-3):(Total_time_size-1)*(Spacing[2]*1e-3)

    #4.4 Unix and Human time vector
    utime = range(start=raw_utime[1], step=(Spacing[2] * 1e-3), length=size(time)[1])
    htime = unix2datetime.(utime)
    replace!(data, NaN => 0.0)
    #5: All to the same array
    dDAS = iDAS(data, time, htime, offset, atrib, filename)

    return dDAS
end

#############################################
#Actual rdas
function rdas(filename::String)

    # Read the file
    printstyled("\nTesting Version of $filename \n", color=:white, bold=true)
    try
        global file = h5open(filename, "r") #Read the fil
    catch e
        printstyled("\n $filename is NOT an .h5 file \n", color=:red)
    end

    ## File information: Name of the DAS and datasets in the .h5 file

    # 1. Get Basic Information
    name = keys(file)[1] #Get the Name of the File

    try
        BlockOverlap = attrs(file["$name/Source1/Zone1"])["BlockOverlap"]
        printstyled("\n $filename is a newer file \n", color=:green)
        global dDAS = rdasnew("$filename")
    catch e
        printstyled("\n No BlockOverlap fild in $filename \n", color=:red)
        printstyled("\n $filename is an older file \n", color=:green)
        global dDAS = rdasold("$filename")
    end

    return dDAS
end


#############################################################################################

"""
writejdl: Write the strain rate data (and its accomapning vectors) to a jdl file.

Input:
 - dDAS: DAS data structure to be written to JDL2 file
 - filename: Name of the file to create. Should be written in air quotes (as it is a String). Please DON'T include extension.

Outputs:
- filename.jdl: A file is created in the working path

Notes:
 - No processing is applied to the data.
 - No figure is created.
 - A jld2 file is a compress version of h5. It can be read in other languages and is much lighter than h5.
# Example: Write dDAS to test.jld2
```
julia> writejdl(dDAS, "test");
```
"""
function writejdl(dDAS, filename::String)
    jldsave("$filename.jld2", true; dDAS)
    printstyled("\n $filename.jld2 created in path \n", color=:green)
    return dDAS
end
#############################################################################################


"""
das2sac: Write the dDAs data to a SAC file for each channel.

Input:
 - dDAS: DAS data structure to be written to JDL2 file
 - x: Offset range or value to transform to sac.

Outputs:
- *.sac: Several .sac files are created in the working path

Notes:
 - No processing is applied to the data.
 - No figure is created.
 - SAC files are created based on Seis.jl library. Big endian is used.
# Example: Write dDAS to sac files
```
julia> das2sac(dDAS; x=100:200);
```
"""
## Show a DAS Channel in time (this is constant X all times)
function das2sac(dDAS; x)

    if typeof(x) == Int64 || typeof(x) == Float64
        channel = argmin(abs.(dDAS.offset .- x))
        schannel = dDAS.offset[channel]
        printstyled(" Saving Channel $channel at distance $schannel m to SAC file\n", color=:yellow)
        signal = dDAS.data[:, channel]
        t = Trace(0.0, dDAS.time[2] - dDAS.time[1], signal)
        t.evt.time = DateTime(year(dDAS.htime[1]), month(dDAS.htime[1]), day(dDAS.htime[1]), hour(dDAS.htime[1]), minute(dDAS.htime[1]), second(dDAS.htime[1]), millisecond(dDAS.htime[1]))
        write_sac(t, "x=$schannel.sac")
    end

    if typeof(x) == UnitRange{Int64} || typeof(x) == StepRange{Int64,Int64} || typeof(x) == StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64} || typeof(x) == StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64}

        for i = 1:size(x, 1)
            channel = argmin(abs.(dDAS.offset .- x[i]))
            schannel = dDAS.offset[channel]
            printstyled(" Saving Channel $channel at distance $schannel m to SAC file\n", color=:yellow)
            signal = dDAS.data[:, channel]
            t = Trace(0.0, dDAS.time[2] - dDAS.time[1], signal)
            t.evt.time = DateTime(year(dDAS.htime[1]), month(dDAS.htime[1]), day(dDAS.htime[1]), hour(dDAS.htime[1]), minute(dDAS.htime[1]), second(dDAS.htime[1]), millisecond(dDAS.htime[1]))
            write_sac(t, "x=$schannel.sac")
        end
    end


end
#############################################################################################


"""
das2txt: Write dDAS data to a ascii file (.txt by deafult).
Note that only the data matrix is written.


Input:
 - dDAS: DAS data structure to be written to JDL2 file
 - filename:  Name of the file to create. Should be written in air quotes (as it is a String). Please DON'T include extension.

Outputs:
- filename.txt: ASCII file is created in the working path.

Notes:
 - No processing is applied to the data.
 - No figure is created.
 - ASCII files might be very large.

# Example: Save one channel and its time vector to a file.
```
julia> das2txt(dDAS,"Export")
```
"""
function das2txt(dDAS, filename::String)
    writedlm("$filename.txt", round.(dDAS.data, digits=5))
end

#############################################################################################


"""
ratrib: Write dDAS data atributes to screen.

Input:
 - dDAS: DAS data structure to be written to JDL2 file

Outputs:
- All to screen.

Notes:
 - No processing is applied to the data.
 - No figure is created.

# Example: Save one channel and its time vector to a file.
```
julia> ratrib(dDAS)
```
"""
function ratrib(dDAS)
    A = dDAS.atrib.AmpliPower
    printstyled("\n AmpliPower: $A \n", color=:green)
    A = dDAS.atrib.BlockOverlap
    printstyled("\n BlockOverlap: $A \n", color=:green)
    A = dDAS.atrib.BlockRate
    printstyled("\n BlockRate: $A \n", color=:green)
    A = dDAS.atrib.Components
    printstyled("\n Components: $A \n", color=:green)
    A = dDAS.atrib.DataDomain
    printstyled("\n DataDomain: $A \n", color=:green)
    A = dDAS.atrib.DerivationTime
    printstyled("\n DerivationTime: $A \n", color=:green)
    A = dDAS.atrib.Extent
    printstyled("\n Extent: $A \n", color=:green)
    A = dDAS.atrib.FiberLength
    printstyled("\n FiberLength: $A \n", color=:green)
    A = dDAS.atrib.GaugeLength
    printstyled("\n GaugeLength: $A \n", color=:green)
    A = dDAS.atrib.Hostname
    printstyled("\n Hostname: $A \n", color=:green)
    A = dDAS.atrib.Origin
    printstyled("\n Origin: $A \n", color=:green)
    A = dDAS.atrib.Oversampling
    printstyled("\n Oversampling: $A \n", color=:green)
    A = dDAS.atrib.PipelineTracker
    printstyled("\n PipelineTracker: $A \n", color=:green)
    A = dDAS.atrib.PulseRateFreq
    printstyled("\n PulseRateFreq: $A \n", color=:green)
    A = dDAS.atrib.PulseWidth
    printstyled("\n PulseWidth: $A \n", color=:green)
    A = dDAS.atrib.SamplingRate
    printstyled("\n SamplingRate: $A \n", color=:green)
    A = dDAS.atrib.SamplingRes
    printstyled("\n SamplingRes: $A \n", color=:green)
    A = dDAS.atrib.Spacing
    printstyled("\n Spacing: $A \n", color=:green)


















end


#############################################################################################
#############################################################################################
# Structure Arrays section
#############################################################################################
#############################################################################################

## Definition of mutalble structure array attb. This is needed to build variable atrib
# with all those attributes!
mutable struct attb
    AmpliPower::Vector{Int8}
    BlockRate::Vector{Int32}
    Components::Vector{Int32}
    DataDomain::Vector{Int32}
    DerivationTime::Vector{Float64}
    Extent::Vector{Int32}
    FiberLength::Vector{Int32}
    GaugeLength::Vector{Float64}
    Hostname::String
    Origin::Vector{Float64}
    Oversampling::Vector{Int8}
    PipelineTracker::String
    PulseRateFreq::Float64
    PulseWidth::Vector{Int8}
    SamplingRate::Vector{Int32}
    SamplingRes::Vector{Int32}
    Spacing::Vector{Float64}
    BlockOverlap::Vector{Int32}
end

mutable struct attb
    AmpliPower::Vector{Int8}
    BlockRate::Vector{Int32}
    Components::Vector{Int32}
    DataDomain::Vector{Int32}
    DerivationTime::Vector{Float64}
    Extent::Vector{Int32}
    FiberLength::Vector{Int32}
    GaugeLength::Vector{Float64}
    Hostname::String
    Origin::Vector{Float64}
    Oversampling::Vector{Int8}
    PipelineTracker::String
    PulseRateFreq::Float64
    PulseWidth::Vector{Int8}
    SamplingRate::Vector{Int32}
    SamplingRes::Vector{Int32}
    Spacing::Vector{Float64}
    BlockOverlap::Vector{Int32}
end


## Definition of mutalble structure array DASdata. Where we save DAS data to files
mutable struct iDAS
    data::Array{Float32,2}
    time::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64}
    htime::Array{DateTime,1}
    offset::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64}
    atrib::attb
    name::String
end
