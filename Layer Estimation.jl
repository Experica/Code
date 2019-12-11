using NeuroAnalysis,Statistics,FileIO,Plots,Interact,Dierckx,PyCall

# Combine all layer tests of one recording site for layer estimation
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"
settimeunit(1)

subject = "AE9";recordsession = "";recordsite = "u003"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
sitedir = joinpath(resultroot,subject,siteid)
layer = Dict("WM"=>[0,0],"Out"=>[3500,3500])

# Plot all CSD
testids = ["$(siteid)_00$i" for i in 0:3]
testn=length(testids)
testtitles = ["Eyeₚ","Eyeₙₚ","Eyes","EyeₚS"]
csds = load.(joinpath.(sitedir,testids,"csd.jld2"),"csd","depth","fs")
depths=csds[1][2];fs=csds[1][3]

csdr = [0.025 0.1]
csdi = epoch2samplerange(csdr,fs)
csdx = ((1:length(csdi))./fs .+ csdr[1])*1000
csdss = map(i->i[1][:,csdi],csds)
clim=maximum(abs.(vcat(csdss...)))

p=plot(layout=(1,testn),link=:all,legend=false)
for i in 1:testn
    heatmap!(p,subplot=i,csdx,depths,csdss[i],color=:RdBu,clims=(-clim,clim),title=testtitles[i])
    if !isnothing(layer)
        hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(csdx[1]+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
end
foreach(i->savefig(joinpath(sitedir,"layer_csd$i")),[".png",".svg"])

# Plot all Power Spectrum
pss = load.(joinpath.(sitedir,testids,"powerspectrum.jld2"),"rcps","depth","freq")
depths=pss[1][2];freq=pss[1][3]

pss = map(i->i[1],pss)
clim=maximum(vcat(pss...))

p=plot(layout=(1,testn),link=:all,legend=false)
for i in 1:testn
    heatmap!(p,subplot=i,freq,depths,pss[i],color=:fire,clims=(0,clim),title=testtitles[i])
    if !isnothing(layer)
        hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(freq[1]+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
end
foreach(i->savefig(joinpath(sitedir,"layer_power_rc$i")),[".png",".svg"])

# Plot all Depth PSTH
depthpsths = load.(joinpath.(sitedir,testids,"depthpsth.jld2"),"depthpsth","depth","x")
depths=depthpsths[1][2];x=depthpsths[1][3];bw = x[2]-x[1]

psthr = [0.025 0.1]
psthi = epoch2samplerange(psthr,1/bw)
psthx = ((1:length(psthi)).*bw .+ psthr[1])*1000
psthss = map(i->i[1][:,psthi],depthpsths)
clim=maximum(vcat(psthss...))

p=plot(layout=(1,testn),link=:all,legend=false)
for i in 1:testn
    heatmap!(p,subplot=i,psthx,depths,psthss[i],color=:Reds,clims=(0,clim),title=testtitles[i])
    if !isnothing(layer)
        hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(psthx[1]+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
end
foreach(i->savefig(joinpath(sitedir,"layer_depthpsth$i")),[".png",".svg"])

# Plot all unit position
spikes = load.(joinpath.(sitedir,testids,"spike.jld2"),"spike")

p=plot(layout=(1,testn),link=:all,legend=false,grid=false,xlims=(10,60))
for i in 1:testn
    scatter!(p,subplot=i,spikes[i]["unitposition"][:,1],spikes[i]["unitposition"][:,2],color=:darkgreen,alpha=0.5,markerstrokewidth=0,markersize=3,title=testtitles[i])
    if !isnothing(layer)
        hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
end
foreach(i->savefig(joinpath(sitedir,"layer_unitposition$i")),[".png",".svg"])












# earliest response should be due to LGN M,P input to 4Ca,4Cb
ln = ["4Cb","4Ca","4B","4A","2/3","Out"]
lcsd = dropdims(mean(mncsd[:,epoch2samplerange([0.045 0.055],fs)],dims=2),dims=2)
ldepths = depths[1]:depths[end]
lcsd = Spline1D(depths,lcsd)(ldepths);threshold = 1.2std(lcsd)

scipysignal = pyimport("scipy.signal")
di,dv=scipysignal.find_peaks(-lcsd,prominence=threshold,height=threshold)
peaks =ldepths[di.+1];bases = hcat(ldepths[dv["left_bases"].+1],ldepths[dv["right_bases"].+1])

plot(lcsd,ldepths,label="CSD Profile")
vline!([-threshold,threshold],label="Threshold")
hline!(peaks,label="CSD Sink Peak")
hline!(bases[:,1],label="CSD Sink Low Border")
hline!(bases[:,2],label="CSD Sink High Border")

layer = Dict(ln[i]=>bases[i,:] for i in 1:size(bases,1))
plotanalog(mncsd,fs=fs,color=:RdBu,layer=layer)


# Layers from Depth PSTH
depthpsth,depths,x = load(joinpath(resultdir,"depthpsth.jld2"),"depthpsth","depth","x")
plotpsth(depthpsth,x,depths,layer=layer)

bw = x[2]-x[1]
lpsth = dropdims(mean(depthpsth[:,epoch2samplerange([0.045 0.055],1/bw)],dims=2),dims=2)
ldepths = depths[1]:depths[end]
lpsth = Spline1D(depths,lpsth)(ldepths);threshold = 1.2std(lpsth)

scipysignal = pyimport("scipy.signal")
di,dv=scipysignal.find_peaks(lpsth,prominence=threshold,height=threshold)
peaks =ldepths[di.+1];bases = hcat(ldepths[dv["left_bases"].+1],ldepths[dv["right_bases"].+1])

plot(lpsth,ldepths,label="PSTH Profile")
vline!([-threshold,threshold],label="Threshold")
hline!(peaks,label="PSTH Peak")
hline!(bases[:,1],label="PSTH Low Border")
hline!(bases[:,2],label="PSTH High Border")












layer["Out"]=[3260,2900]
layer["2/3"]=[2645,1800]
layer["4A"]=[2550,1800]
layer["4B"]=[2460,1800]
layer["4Ca"]=[2250,1800]
layer["4Cb"]=[2100,1300]
layer["5"]=[1800,1800]
layer["6"]=[1500,1800]
layer["WM"]=[0,0]
# Finalize Layers
save(joinpath(sitedir,"layer.jld2"),"layer",checklayer(layer))



# Layer Verification with Tuning Properties
layer = load(joinpath(sitedir,"layer.jld2"),"layer")
testids = ["$(siteid)_$(lpad(i,3,'0'))" for i in [8,12,13,14]]
testn=length(testids)
testtitles = ["Lum","L","M","S"]
ds = load.(joinpath.(sitedir,testids,"factorresponse.jld2"),"factorstats","fms","fses","fa","responsive","modulative")
spikes = load.(joinpath.(sitedir,testids,"spike.jld2"),"spike")

f = :Dir
p=plot(layout=(1,testn),link=:all,legend=false,grid=false,xlims=(10,60))
for i in 1:testn
    if f==:Ori
        color = map((i,j)->j ? HSV(i.oo,1-i.ocv,1) : HSVA(0,0,0.5,0.2),ds[i][1][:Ori],ds[i][end])
    elseif f==:Dir
        color = map((i,j)->j ? HSV(i.od,1-i.dcv,1) : HSVA(0,0,0.5,0.2),ds[i][1][:Ori],ds[i][end])
    end
    scatter!(p,subplot=i,spikes[i]["unitposition"][:,1],spikes[i]["unitposition"][:,2],color=color,markerstrokewidth=0,markersize=3,title=testtitles[i])
    if !isnothing(layer)
        hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
end
foreach(i->savefig(joinpath(sitedir,"layer_unitposition_$(f)_Tuning$i")),[".png",".svg"])




pyplot()
