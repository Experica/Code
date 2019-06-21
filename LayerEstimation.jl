using NeuroAnalysis,Statistics,FileIO,Plots,Dierckx,PyCall

# Combine all layer estimations for one recording site
env = Dict{Any,Any}(
    :dataroot => "../Data",
    :dataexportroot => "../DataExport",
    :resultroot => "../Result",
    :imageroot => "../NaturalStimuli")
settimeunit(1)
subject = "AE9";recordsession = "";recordsite = "u003";
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
sitedir = joinpath(env[:resultroot],subject,siteid)
layer = nothing


testid = "$(siteid)_000"
# Layers from CSD
resultdir = joinpath(sitedir,testid)
mncsd,depths,fs = load(joinpath(resultdir,"csd.jld2"),"csd","depth","fs")
plotanalog(mncsd,fs=fs,color=:RdBu,layer=layer)

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




layer["Out"]=[3160,2900]
layer["2/3"]=[2675,1800]
layer["4A"]=[2565,1800]
layer["4B"]=[2390,1800]
layer["4Ca"]=[2150,1800]
layer["4Cb"]=[1340,1300]
layer["5"]=[1500,1800]
layer["6"]=[1270,1800]
layer["WM"]=[0,0]
# Plot Layers with all CSD
pcsd,depths,fs = load(joinpath(sitedir,"$(siteid)_000","csd.jld2"),"csd","depth","fs")
npcsd,depths,fs = load(joinpath(sitedir,"$(siteid)_001","csd.jld2"),"csd","depth","fs")
bcsd,depths,fs = load(joinpath(sitedir,"$(siteid)_002","csd.jld2"),"csd","depth","fs")
pscsd,depths,fs = load(joinpath(sitedir,"$(siteid)_003","csd.jld2"),"csd","depth","fs")

csdeidx = epoch2samplerange([0.025 0.1],fs)
csdx = ((1:length(csdeidx))./fs .+ 0.025)*1000
epcsd = pcsd[:,csdeidx]
enpcsd = npcsd[:,csdeidx]
ebcsd = bcsd[:,csdeidx]
epscsd = pscsd[:,csdeidx]

lcsdx = csdx[1]
csdlim=maximum(abs.(vcat(epcsd,enpcsd,ebcsd,epscsd)))
p=plot(layout=(1,4),link=:all,legend=false)
heatmap!(p,subplot=1,csdx,depths,epcsd,color=:RdBu,clims=(-csdlim,csdlim),title="Prefered Eye")
hline!(p,subplot=1,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lcsdx+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
heatmap!(p,subplot=2,csdx,depths,enpcsd,color=:RdBu,clims=(-csdlim,csdlim),title="Non-Prefered Eye")
hline!(p,subplot=2,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lcsdx+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
heatmap!(p,subplot=3,csdx,depths,ebcsd,color=:RdBu,clims=(-csdlim,csdlim),title="Both Eye")
hline!(p,subplot=3,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lcsdx+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
heatmap!(p,subplot=4,csdx,depths,epscsd,color=:RdBu,clims=(-csdlim,csdlim),title="Prefered Eye Siso")
hline!(p,subplot=4,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lcsdx+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)

foreach(i->savefig(joinpath(sitedir,"layer_csd$i")),[".png",".svg"])


# Plot layers with all DepthPSTH
pdepthpsth,depths,x = load(joinpath(sitedir,"$(siteid)_000","depthpsth.jld2"),"depthpsth","depth","x")
npdepthpsth,depths,x = load(joinpath(sitedir,"$(siteid)_001","depthpsth.jld2"),"depthpsth","depth","x")
bdepthpsth,depths,x = load(joinpath(sitedir,"$(siteid)_002","depthpsth.jld2"),"depthpsth","depth","x")
psdepthpsth,depths,x = load(joinpath(sitedir,"$(siteid)_003","depthpsth.jld2"),"depthpsth","depth","x")

dpeidx = epoch2samplerange([0.025 0.1],1/bw)
dpx = ((1:length(dpeidx)).*bw .+ 0.025)*1000
epdepthpsth = pdepthpsth[:,dpeidx]
enpdepthpsth = npdepthpsth[:,dpeidx]
ebdepthpsth = bdepthpsth[:,dpeidx]
epsdepthpsth = psdepthpsth[:,dpeidx]

ldpx = dpx[1]
dpclim=maximum(vcat(epdepthpsth,enpdepthpsth,ebdepthpsth,epsdepthpsth))
p=plot(layout=(1,4),link=:all,legend=false)
heatmap!(p,subplot=1,dpx,depths,epdepthpsth,color=:Reds,clims=(0,dpclim),title="Prefered Eye")
hline!(p,subplot=1,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(ldpx+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
heatmap!(p,subplot=2,dpx,depths,enpdepthpsth,color=:Reds,clims=(0,dpclim),title="Non-Prefered Eye")
hline!(p,subplot=2,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(ldpx+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
heatmap!(p,subplot=3,dpx,depths,ebdepthpsth,color=:Reds,clims=(0,dpclim),title="Both Eye")
hline!(p,subplot=3,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(ldpx+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
heatmap!(p,subplot=4,dpx,depths,epsdepthpsth,color=:Reds,clims=(0,dpclim),title="Prefered Eye Siso")
hline!(p,subplot=4,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(ldpx+5,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)

foreach(i->savefig(joinpath(sitedir,"layer_depthpsth$i")),[".png",".svg"])


# Plot layers with all unitposition
pspike = load(joinpath(sitedir,"$(siteid)_000","spike.jld2"),"spike")
npspike = load(joinpath(sitedir,"$(siteid)_001","spike.jld2"),"spike")
bspike = load(joinpath(sitedir,"$(siteid)_002","spike.jld2"),"spike")
psspike = load(joinpath(sitedir,"$(siteid)_003","spike.jld2"),"spike")

p=plot(layout=(1,4),link=:all,xlims=(10,60),grid=false)
scatter!(p,subplot=1,pspike["unitposition"][:,1],pspike["unitposition"][:,2],color=:darkgreen,alpha=0.5,markerstrokewidth=0,markersize=3)
hline!(p,subplot=1,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
scatter!(p,subplot=2,npspike["unitposition"][:,1],npspike["unitposition"][:,2],color=:darkgreen,alpha=0.5,markerstrokewidth=0,markersize=3)
hline!(p,subplot=2,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
scatter!(p,subplot=3,bspike["unitposition"][:,1],bspike["unitposition"][:,2],color=:darkgreen,alpha=0.5,markerstrokewidth=0,markersize=3)
hline!(p,subplot=3,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
scatter!(p,subplot=4,psspike["unitposition"][:,1],psspike["unitposition"][:,2],color=:darkgreen,alpha=0.5,markerstrokewidth=0,markersize=3)
hline!(p,subplot=4,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(15,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)

foreach(i->savefig(joinpath(sitedir,"layer_unitposition$i")),[".png",".svg"])


# Finalize Layers
save(joinpath(sitedir,"layer.jld2"),"layer",checklayer(layer))
