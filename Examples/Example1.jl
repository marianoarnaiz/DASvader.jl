"DASVader V0.1 Example 1. Here we try to introduce the most basic concepts and tools of DASVader"

#Step 0: Load the Package
using DASVader

#Step 1: Read a raw Febus HDF5 file to memory.
#We will use the rdas command to read the name of the file.
# rdas will output the data of the file to a variable called "dDAS".
#Notive that you must use "" to indicate the name of the file.
#You should be able to read any file comming from the FEBUS A1 DAS.
dDAS = rdas("SR_DS_2023-10-30_12-01-40_UTC_Microevent.h5")

#Step 2: View the data on the file.
# We will use the command viewdas to plot the color matrix of the data in dDAS.
#You will need to input the name of the variable to plot (in this case dDAS).
#There are also other inputs required:
# - cm is the colormap to use. THere are many of them: a few options are: ":grays", ":viridis" and ":RdBu_9"
# if you need more options go to this website: https://docs.makie.org/dev/explanations/colors
# - climit is the end of the color bar to use. Everything above and bellow this value will be saturated. Try 
# several ranges to see your data more clearly.
fig = viewdas(dDAS; cm=:RdBu_9, climit=10000)
# Why not save that figure
savefig(fig, "MicroEvent.pdf")

#Step 3: Cut the Event with slicedas. This will require you to input the time range (as tmin and tmax) as
# well as the offset range (xmin and xmax) to cut the data.

cDAS = slicedas(dDAS; tmin=4.65, tmax=4.85, xmin=1250, xmax=1400)

# and now we can make another figure
fig2 = viewdas(cDAS; cm=:RdBu_9, climit=10000)
savefig(fig2, "MicroEvent_Zoom.pdf")

#Step 4: See the data a record section!
#This is very straight forward, just input the variable of the data
#and a scaling factor that is related to the offset
rs= recordsection(cDAS, scale=5.0)
savefig(rs, "Record_Section.pdf")

#Step 5. Study the waveform of one channel with a bit more detail!
# For this we use viewchannel. Here you need the data to use and the 
#offset of the channel. 
channel1, spectrum1, spectogram1, f_Channel1 = viewchannel(dDAS; x=1326.5)

#Step 6: Plot the channel 1 with a time vector and add the next and previous channel!
ilines(dDAS.time, channel1)
ilines!(dDAS.time, dDAS.data[:,829])
ilines!(dDAS.time, dDAS.data[:,831])

# Step 7: save the selected channels to SAC 
#To save the interesting channels to a bunch of SAC files can be done if you know the 
#range of offsets you with to export
das2sac(dDAS; x=1310:1340)

#Step 8: Write the data in the structure to a JLD2 file (compressed H5)
writejdl(cDAS, "Test")
#@load "Test.jld2"

