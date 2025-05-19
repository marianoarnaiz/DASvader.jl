using GLMakie, CairoMakie, Dates

#################################################################
# Main FUNCTIONS
#

function viewdas(dDAS; cm=:RdBu_11, climit=10000, picks=[])
    GLMakie.activate!() # Activate GL plotting
    fig = Figure(; size=(1000, 600))

    # Add the plotting elements
    ihm = iheatmap(fig[1, 1], dDAS.time, dDAS.offset, dDAS.data, colormap=cm, colorrange=(-climit, climit))
    Colorbar(fig[1, 2], ihm.plot, label="StrainRate [nm/m/s]")
    ihm.axis.xlabel = "Relative Time [s]"
    ihm.axis.ylabel = "Offset [m]"
    ihm.axis.title = dDAS.name
    DataInspector(fig)

    if !isempty(picks)
        scatter!(ihm.axis, picks[:, 4], picks[:, 2], marker='|', color=:black, markersize=20)
        println("Plotting picks")
    end

    # Add the buttons
    add_buttons!(fig)

    # Display the figure
    display(fig)
    return fig
end


function viewchannel(dDAS; x, type=identity)
    GLMakie.activate!() # Activate GL plotting
    chind = argmin(abs.(dDAS.offset .- x)) # Get the channel index
    println("\nPlotting channel $x at index $chind\n")

    ich = dDAS.offset[chind]

    fv, smv = get_spectrum(dDAS.data[:, chind], 1 / (dDAS.time[2] - dDAS.time[1]))

    #GET SPECTOGRAM
    spec = gspectrogram(dDAS.data[:, chind], 1 / delta)



    fig = Figure(; size=(600, 1000))
    ax1 = Axis(fig[1, 1], xlabel="Relative Time [s]", ylabel="ηε", title="Channel $chind at offset $ich")
    ilines!(ax1, dDAS.time, dDAS.data[:, chind])
    ax2 = Axis(fig[2, 1][1, 1], xlabel="Relative Time [s]", ylabel="Frequency [Hz]", title="Spectrogram")
    cf = contourf!(ax2, spec.time .+ dDAS.time[1], spec.freq, sqrt.(spec.amp'), levels=25, colormap=:plasma)
    Colorbar(fig[2, 1][1, 2], cf, label="dB")
    linkxaxes!(ax1, ax2)
    ax3 = Axis(fig[3, 1], xlabel="Frequency [Hz]", ylabel="Power", title="Frequency Spectrum", xscale=log10, yscale=log10)
    ilines!(ax3, fv, abs.(smv) .^ 2, color=:black)


    # Add the buttons
    add_buttons!(fig)

    # Display the figure
    display(fig)
    return fig
end

function recordsection(dDAS; scale=1.0)

    GLMakie.activate!() #activate GL plotting
    dx = dDAS.offset[2] - dDAS.offset[1]
    n = length(dDAS.time)
    scaleddata = dDAS.data ./ maximum(abs.(dDAS.data))
    scaleddata = scale * dx * scaleddata
    nc = size(dDAS.data, 2)
    # Plot signals with positive/negative areas filled:

    fig = Figure(; size=(1000, 1000))
    ax = Axis(fig[1, 1], title=dDAS.name, xlabel="Time [s]", ylabel="Offset [m]", xminorgridvisible=true, yminorgridvisible=true, limits=(nothing, nothing))
    for i in 1:nc
        #band!(ax, time, map(x-> x < 0 ? x + i*dx : i*dx, strainrate[i,:]), i*dx, color=("red", 0.3))
        #band!(time, map(x-> x >= 0 ? x + i*dx : i*dx, strainrate[i,:]), i*dx, color=("blue", 0.3))
        ilines!(ax, dDAS.time, scaleddata[:, i] .+ i * dx .+ dDAS.offset[1], linewidth=0.5, color=:black)
    end
    tightlimits!(ax)

    # Add the buttons
    add_buttons!(fig)

    # Display the figure
    display(fig)

    #fig
    return fig
end



#################################################################
# AUXILIARY FUNCTIONS
#

# Function to save current figure as GLMakie snapshot
function savescreenshot(fig, name)
    GLMakie.save("$name.png", fig.scene)
    println("Spanshot saved as $name.png \n")
end

# Function to save figure via CairoMakie as PDF
function savefig(fig, name)
    #using CairoMakie
    save("$name", fig, backend=CairoMakie)
    println("Figure saved as $name \n")
    return
end

# Function to add buttons to the figure
function add_buttons!(fig)
    # Create grid layout for buttons
    buttongrid = GridLayout(tellwidth=false)
    fig[end+1, 1] = buttongrid

    # Save Snapshot button
    screenshotbutton = Button(fig, label="Save Snapshot")
    on(screenshotbutton.clicks) do _
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")
        savescreenshot(fig, "snapshot_$timestamp")
    end

    # Save PDF button
    savebutton = Button(fig, label="Save Figure")
    on(savebutton.clicks) do _
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")
        savefig(fig, "Figure_$timestamp.pdf")
    end

    # Close button
    closebutton = Button(fig, label="Close")
    screen = display(fig)
    on(closebutton.clicks) do _
        GLMakie.close(screen)
        println("Figure window closed. \n")
    end

    # Add buttons to the grid layout
    buttongrid[1, 1:3] = [screenshotbutton, savebutton, closebutton]
end



export viewdas, recordsection, viewchannel, savescreenshot, savefig, add_buttons!
