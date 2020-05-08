using NeuroAnalysis,Statistics,FileIO,Plots,Interact#,Dierckx,PyCall

# Combine all layer tests of one recording site
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL3"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
resultsitedir = joinpath(resultroot,subject,siteid)

layer = Dict("WM"=>[0,0],"Out"=>[3500,3500])

# testids = ["$(siteid)_00$i" for i in 0:3]
# testn=length(testids)
# testtitles = ["Eyeₚ","Eyeₙₚ","Eyes","EyeₚS"]
# csds = load.(joinpath.(sitedir,testids,"csd.jld2"),"csd","depth","fs")
# depths=csds[1][2];fs=csds[1][3]

testids = ["$(siteid)_Flash2Color_$i" for i in 1:4]
testn=length(testids)
## All CSD
csds = load.(joinpath.(resultsitedir,testids,"csd.jld2"))
testlogs = ["$(i["log"])_$(i["color"])" for i in csds]
fs=csds[1]["fs"]
csddepths=csds[1]["depth"]
csdtimes = csds[1]["time"]

csdb = [0 15]
csdr = [15 100]
csdbi = epoch2sample(csdb,fs)
csdri = epoch2sample(csdr,fs)
csdx = csdtimes[csdri]
csdss = map(i->imfilter(stfilter(i["csd"],temporaltype=:sub,ti=csdbi,hbordervalue=0),Kernel.gaussian((0.9,1)))[:,csdri],csds)
csdclim=maximum(abs.(vcat(csdss...)))

p=plot(layout=(1,testn),link=:all,legend=false,grid=false,size=(testn*200,700))
for i in 1:testn
    yticks = i==1 ? true : false
    xlabel = i==1 ? "Time (ms)" : ""
    ylabel = i==1 ? "Depth (um)" : ""
    heatmap!(p,subplot=i,csdx,csddepths,csdss[i],color=:RdBu,clims=(-csdclim,csdclim),title=testlogs[i],titlefontsize=8,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
    if !isnothing(layer)
        hline!(p,subplot=i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray30,legend=false,
        annotations=[(csdx[1]+1,layer[k][1],text(k,4,:gray20,:bottom,:left)) for k in keys(layer)])
    end
end
p
foreach(i->savefig(joinpath(resultsitedir,"Layer_dCSD$i")),[".png",".svg"])

## All Depth PSTH
depthpsths = load.(joinpath.(resultsitedir,testids,"depthpsth.jld2"))
psthtimes=depthpsths[1]["time"];bw = psthtimes[2]-psthtimes[1];psthdepths = depthpsths[1]["depth"]

psthb = [0 15]
psthr = [15 100]
psthbi = epoch2sample(psthb,1/(bw*SecondPerUnit))
psthri = epoch2sample(psthr,1/(bw*SecondPerUnit))
psthx = psthtimes[psthri]
psthss = map(i->imfilter(stfilter(i["depthpsth"],temporaltype=:sub,ti=psthbi),Kernel.gaussian((0.9,1)))[:,psthri],depthpsths)
psthclims=extrema(vcat(psthss...))

p=plot(layout=(1,testn),link=:all,legend=false,grid=false,size=(testn*200,700))
for i in 1:testn
    yticks = i==1 ? true : false
    xlabel = i==1 ? "Time (ms)" : ""
    ylabel = i==1 ? "Depth (um)" : ""
    heatmap!(p,subplot=i,psthx,psthdepths,psthss[i],color=:Reds,clims=psthclims,title=testlogs[i],titlefontsize=8,yticks=yticks,xlabel=xlabel,ylabel=ylabel,ylims=extrema(depths),xlims=extrema(psthx))
    if i==testn
        n = depthpsths[i]["n"]
        pn = maximum(psthx) .- n./maximum(n) .* maximum(psthx) .* 0.2
        plot!(p,subplot=testn,pn,depths,label="Number of Units",color=:seagreen,lw=0.5)
    end
    if !isnothing(layer)
        hline!(p,subplot=i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray30,legend=false,
        annotations=[(psthx[1]+1,layer[k][1],text(k,4,:gray20,:bottom,:left)) for k in keys(layer)])
    end
end
p
foreach(i->savefig(joinpath(resultsitedir,"Layer_dPSTH$i")),[".png",".svg"])


## All Power Spectrum
pss = load.(joinpath.(resultsitedir,testids,"powerspectrum.jld2"))
psdepths=pss[1]["depth"];freq=pss[1]["freq"]

psss = map(i->i["rcps"],pss)
psclims=extrema(vcat(psss...))

p=plot(layout=(1,testn),link=:all,legend=false,grid=false,size=(testn*200,700))
for i in 1:testn
    yticks = i==1 ? true : false
    xlabel = i==1 ? "Frequency (Hz)" : ""
    ylabel = i==1 ? "Depth (um)" : ""
    heatmap!(p,subplot=i,freq,psdepths,psss[i],color=:PuRd,clims=psclims,title=testlogs[i],titlefontsize=8,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
    if !isnothing(layer)
        hline!(p,subplot=i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray30,legend=false,
        annotations=[(freq[1]+1,layer[k][1],text(k,4,:gray20,:bottom,:left)) for k in keys(layer)])
    end
end
p
foreach(i->savefig(joinpath(resultsitedir,"Layer_PowerSpectrum_RelativeChange$i")),[".png",".svg"])

## Plot all unit position
spikes = load.(joinpath.(resultsitedir,testids,"spike.jld2"),"spike")

p=plot(layout=(1,testn),link=:all,legend=false,grid=false,xlims=(10,60),size=(testn*200,700))
for i in 1:testn
    yticks = i==1 ? true : false
    xlabel = i==1 ? "Position_X (um)" : ""
    ylabel = i==1 ? "Position_Y (um)" : ""
    scatter!(p,subplot=i,spikes[i]["unitposition"][:,1],spikes[i]["unitposition"][:,2],color=map(i->i ? :darkgreen : :gray30,spikes[i]["unitgood"]),yticks=yticks,xlabel=xlabel,ylabel=ylabel,
    alpha=0.5,markerstrokewidth=0,markersize=3,series_annotations=text.(spikes[i]["unitid"],2,:gray10,:center),title=testlogs[i],titlefontsize=8)
    if !isnothing(layer)
        hline!(p,subplot=i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray30,legend=false,
        annotations=[(11,layer[k][1],text(k,4,:gray20,:bottom,:left)) for k in keys(layer)])
    end
end
p
foreach(i->savefig(joinpath(resultsitedir,"Layer_UnitPosition$i")),[".png",".svg"])











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












layer["Out"]=[3170,3900]
layer["1"]=[2160,1800]
layer["2/3"]=[2670,1800]
layer["2"]=[2550,1800]
layer["3"]=[2550,1800]
layer["4A/B"]=[2575,1800]
layer["4A"]=[2575,1800]
layer["4B"]=[2430,1800]
layer["4C"]=[1800,1800]
layer["4Ca"]=[2170,1800]
layer["4Cb"]=[1870,1300]
layer["5/6"]=[1510,1800]
layer["5"]=[1300,1800]
layer["6"]=[1500,1800]
layer["WM"]=[0,0]
# Finalize Layers
save(joinpath(resultsitedir,"layer.jld2"),"layer",checklayer(layer))



## Tuning Properties in layers
layer = load(joinpath(resultsitedir,"layer.jld2"),"layer")
# testids = ["$(siteid)_$(lpad(i,3,'0'))" for i in [8,12,13,14]]

testids = ["$(siteid)_OriSF_$i" for i in 1:5]
testn=length(testids)
ds = load.(joinpath.(resultsitedir,testids,"factorresponse.jld2"))
testtitles = ["$(i["color"])" for i in ds]
spikes = load.(joinpath.(resultsitedir,testids,"spike.jld2"),"spike")

f = :Ori
p=plot(layout=(1,testn),link=:all,legend=false,grid=false,xlims=(10,60))
for i in 1:testn
    vi = spikes[i]["unitgood"].&ds[i]["responsive"].&ds[i]["modulative"]
    if f in [:Ori,:Ori_Final]
        color = map((i,j)->j ? HSV(i.oo,1-i.ocv/1.5,1) : RGBA(0.5,0.5,0.5,0.1),ds[i]["factorstats"][:Ori_Final],vi)
    elseif f==:Dir
        color = map((i,j)->j ? HSV(i.od,1-i.dcv/1.5,1) : RGBA(0.5,0.5,0.5,0.1),ds[i]["factorstats"][:Ori_Final],vi)
    end
    scatter!(p,subplot=i,spikes[i]["unitposition"][:,1],spikes[i]["unitposition"][:,2],color=color,markerstrokewidth=0,markersize=2,title=testtitles[i],titlefontsize=8)
    if !isnothing(layer)
        hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(17,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
end
p
foreach(i->savefig(joinpath(resultsitedir,"Layer_UnitPosition_$(f)_Tuning$i")),[".png",".svg"])
