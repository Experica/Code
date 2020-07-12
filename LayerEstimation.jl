using NeuroAnalysis,Statistics,FileIO,Plots,Images,Interact

# Combine all layer tests of one recording site
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL5"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
resultsitedir = joinpath(resultroot,subject,siteid)
figfmt = [".png",".svg"]

layer = load(joinpath(resultsitedir,"layer.jld2"),"layer")
layer = Dict("WM"=>[0,0],"Out"=>[3800,3800])
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
csdbi = epoch2sampleindex(csdb,fs)
csdri = epoch2sampleindex(csdr,fs)
csdx = csdtimes[csdri]
csdss = map(i->imfilter(stfilter(i["csd"],temporaltype=:sub,ti=csdbi,hbordervalue=0),Kernel.gaussian((1,1)))[:,csdri],csds)
csdss1 = map(i->imfilter(stfilter(i["ccsd"][1],temporaltype=:sub,ti=csdbi,hbordervalue=0),Kernel.gaussian((1,1)))[:,csdri],csds)
csdss2 = map(i->imfilter(stfilter(i["ccsd"][2],temporaltype=:sub,ti=csdbi,hbordervalue=0),Kernel.gaussian((1,1)))[:,csdri],csds)
absmax(x) = mapreduce(i->maximum(abs.(i)),max,x)
csdclim=max(absmax(csdss),absmax(csdss1),absmax(csdss2))

## All Depth PSTH
depthpsths = load.(joinpath.(resultsitedir,testids,"depthpsth.jld2"))
psthtimes=depthpsths[1]["time"];bw = psthtimes[2]-psthtimes[1];psthdepths = depthpsths[1]["depth"];depthnunit=depthpsths[end]["n"]

psthb = [0 15]
psthr = [15 100]
psthbi = epoch2sampleindex(psthb,1/(bw*SecondPerUnit))
psthri = epoch2sampleindex(psthr,1/(bw*SecondPerUnit))
psthx = psthtimes[psthri]
psthss = map(i->imfilter(stfilter(i["depthpsth"],temporaltype=:sub,ti=psthbi),Kernel.gaussian((1,1)))[:,psthri],depthpsths)
psthss1 = map(i->imfilter(stfilter(i["cdepthpsth"][1],temporaltype=:sub,ti=psthbi),Kernel.gaussian((1,1)))[:,psthri],depthpsths)
psthss2 = map(i->imfilter(stfilter(i["cdepthpsth"][2],temporaltype=:sub,ti=psthbi),Kernel.gaussian((1,1)))[:,psthri],depthpsths)
psthclim=max(absmax(psthss),absmax(psthss1),absmax(psthss2))

## All Power Spectrum
pcs = load.(joinpath.(resultsitedir,testids,"powercontrast.jld2"))
pcdepths=pcs[1]["depth"];freq=pcs[1]["freq"]

pcss = map(i->imfilter(i["pcs"],Kernel.gaussian((1,1))),pcs)
pcss1 = map(i->imfilter(i["cpcs"][1],Kernel.gaussian((1,1))),pcs)
pcss2 = map(i->imfilter(i["cpcs"][2],Kernel.gaussian((1,1))),pcs)
pcclim=max(absmax(pcss),absmax(pcss1),absmax(pcss2))

## Set Layers
plotlayer=(o...;w=150,h=700,csd=csdss1,psth=psthss1,pc=pcss1)->begin
    p=plot(layout=(1,3testn),link=:y,legend=false,grid=false,size=(3testn*w,h))
    for i in 1:testn
        yticks = i==1 ? true : false
        xlabel = i==1 ? "Time (ms)" : ""
        ylabel = i==1 ? "Depth (um)" : ""
        heatmap!(p,subplot=i,csdx,csddepths,csd[i],color=:RdBu,clims=(-csdclim,csdclim),title=testlogs[i],titlefontsize=8,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
        if !isnothing(layer)
            hline!(p,subplot=i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
            annotations=[(csdx[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
        end
    end
    for i in 1:testn
        yticks = false
        xlabel = i==1 ? "Time (ms)" : ""
        ylabel = ""
        heatmap!(p,subplot=testn+i,psthx,psthdepths,psth[i],color=:coolwarm,clims=(-psthclim,psthclim),title=testlogs[i],titlefontsize=8,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
        if i==testn
            pn = maximum(psthx) .- depthnunit./maximum(depthnunit) .* maximum(psthx) .* 0.2
            plot!(p,subplot=testn+i,pn,psthdepths,label="Number of Units",color=:seagreen,lw=0.5)
        end
        if !isnothing(layer)
            hline!(p,subplot=testn+i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
            annotations=[(psthx[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
        end
    end
    for i in 1:testn
        yticks = false
        xlabel = i==1 ? "Frequency (Hz)" : ""
        ylabel = ""
        heatmap!(p,subplot=2testn+i,freq,pcdepths,pc[i],color=:vik,clims=(-pcclim,pcclim),title=testlogs[i],titlefontsize=8,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
        if !isnothing(layer)
            hline!(p,subplot=2testn+i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
            annotations=[(freq[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
        end
    end
    p
end

lw = Dict(k=>widget(0:3800,label=k,value=layer[k][1]) for k in keys(layer))
foreach(k->on(v->layer[k][1]=v,lw[k]),keys(lw))
lp = map(plotlayer,values(lw)...)
vbox(values(lw)...,lp)

plotlayer(csd=csdss,psth=psthss,pc=pcss)
foreach(ext->savefig(joinpath(resultsitedir,"Layer_dCSD_dPSTH_PowerContrast$ext")),figfmt)
plotlayer(csd=csdss1,psth=psthss1,pc=pcss1)
foreach(ext->savefig(joinpath(resultsitedir,"Cond_1_Layer_dCSD_dPSTH_PowerContrast$ext")),figfmt)
plotlayer(csd=csdss2,psth=psthss2,pc=pcss2)
foreach(ext->savefig(joinpath(resultsitedir,"Cond_2_Layer_dCSD_dPSTH_PowerContrast$ext")),figfmt)

## Finalize Layers
save(joinpath(resultsitedir,"layer.jld2"),"layer",checklayer!(layer))






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




# # earliest response should be due to LGN M,P input to 4Ca,4Cb
# ln = ["4Cb","4Ca","4B","4A","2/3","Out"]
# lcsd = dropdims(mean(mncsd[:,epoch2samplerange([0.045 0.055],fs)],dims=2),dims=2)
# ldepths = depths[1]:depths[end]
# lcsd = Spline1D(depths,lcsd)(ldepths);threshold = 1.2std(lcsd)
#
# scipysignal = pyimport("scipy.signal")
# di,dv=scipysignal.find_peaks(-lcsd,prominence=threshold,height=threshold)
# peaks =ldepths[di.+1];bases = hcat(ldepths[dv["left_bases"].+1],ldepths[dv["right_bases"].+1])
#
# plot(lcsd,ldepths,label="CSD Profile")
# vline!([-threshold,threshold],label="Threshold")
# hline!(peaks,label="CSD Sink Peak")
# hline!(bases[:,1],label="CSD Sink Low Border")
# hline!(bases[:,2],label="CSD Sink High Border")
#
# layer = Dict(ln[i]=>bases[i,:] for i in 1:size(bases,1))
# plotanalog(mncsd,fs=fs,color=:RdBu,layer=layer)
#
#
# # Layers from Depth PSTH
# depthpsth,depths,x = load(joinpath(resultdir,"depthpsth.jld2"),"depthpsth","depth","x")
# plotpsth(depthpsth,x,depths,layer=layer)
#
# bw = x[2]-x[1]
# lpsth = dropdims(mean(depthpsth[:,epoch2samplerange([0.045 0.055],1/bw)],dims=2),dims=2)
# ldepths = depths[1]:depths[end]
# lpsth = Spline1D(depths,lpsth)(ldepths);threshold = 1.2std(lpsth)
#
# scipysignal = pyimport("scipy.signal")
# di,dv=scipysignal.find_peaks(lpsth,prominence=threshold,height=threshold)
# peaks =ldepths[di.+1];bases = hcat(ldepths[dv["left_bases"].+1],ldepths[dv["right_bases"].+1])
#
# plot(lpsth,ldepths,label="PSTH Profile")
# vline!([-threshold,threshold],label="Threshold")
# hline!(peaks,label="PSTH Peak")
# hline!(bases[:,1],label="PSTH Low Border")
# hline!(bases[:,2],label="PSTH High Border")




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
