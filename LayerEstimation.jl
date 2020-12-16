using NeuroAnalysis,Statistics,FileIO,Plots,Images,Interact

# Combine all layer tests of one recording site
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL3"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)
figfmt = [".png",".svg"]

layer = load(joinpath(siteresultdir,"layer.jld2"),"layer")
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


df=DataFrame(power=[405.5,239.0,92.8,410.5,115.2,493.6],hue=repeat(["180","0"],outer=3),location=repeat(["P5","P4","P3"],inner=2))

df |> @vlplot(:bar,x={:hue,axis={labelAngle=0}},y={:power,axis={grid=false}},column=:location,color={:hue,scale={range=["#FF5A81","#00A57E"]}})


# P5
mm=405.5
lm=239.0
# P4
mm=92.8
lm=410.5
# P3
mm = 115.2
lm = 493.6

lfp=load(joinpath(siteresultdir,"$(siteid)_Flash2Color_2","lfp.jld2"))
freq=lfp["freq"]
depths = -lfp["depth"].+layer["Out"][1]
pc = lfp["cmpc"]
mpc = pc["Color=[0.0, 0.6462, 0.4939, 1.0]"]
lpc = pc["Color=[1.0, 0.3538, 0.5061, 1.0]"]
li = 0 .<=depths.<=300
fi = 30 .<=freq.<=100

mm = sum(mpc[li,fi])
lm = sum(lpc[li,fi])


lfp=load(joinpath(siteresultdir,"$(siteid)_Flash2Color_3","lfp.jld2"))
freq=lfp["freq"]
depths = -lfp["depth"].+layer["Out"][1]
pc = lfp["cmpc"]
lmpc = pc["Color=[0.3672, 0.6168, 0.0, 1.0]"]
spc = pc["Color=[0.6328, 0.3832, 1.0, 1.0]"]
li = 100 .<=depths.<=400
fi = 30 .<=freq.<=100

lmm = mean(lmpc[li,fi])
sm = mean(spc[li,fi])

# testids = ["$(siteid)_00$i" for i in 0:3]
# testn=length(testids)
# testtitles = ["Eyeₚ","Eyeₙₚ","Eyes","EyeₚS"]
# csds = load.(joinpath.(sitedir,testids,"csd.jld2"),"csd","depth","fs")
# depths=csds[1][2];fs=csds[1][3]

testids = ["$(siteid)_Flash2Color_$i" for i in 1:4]
testn=length(testids)
## All CSD and Power Contrast from LFP
tlfps = load.(joinpath.(siteresultdir,testids,"lfp.jld2"))
testlogs = ["$(i["log"])_$(i["color"])" for i in tlfps]
fs=tlfps[1]["fs"]
freq=tlfps[1]["freq"]
lfpdepths=tlfps[1]["depth"]
lfptimes = tlfps[1]["time"]

csdb = [0 15]
csdr = [15 100]
csdbi = epoch2sampleindex(csdb,fs)
csdri = epoch2sampleindex(csdr,fs)
csdx = lfptimes[csdri]
csds = mapreduce((l,i)->["$(l)_$k"=>imfilter(stfilter(v,temporaltype=:sub,ti=csdbi,hbordervalue=0),Kernel.gaussian((1,1)))[:,csdri] for (k,v) in i["cmcsd"]],
append!,testlogs,tlfps)
absmax(x) = mapreduce(i->maximum(abs.(i.second)),max,x)
csdclim=absmax(csds)

pcs = mapreduce((l,i)->["$(l)_$k"=>imfilter(v,Kernel.gaussian((1,1))) for (k,v) in i["cmpc"]],
append!,testlogs,tlfps)
pcclim=absmax(pcs)
## All PSTH
tpsths = load.(joinpath.(siteresultdir,testids,"depthpsth.jld2"))
ffpsth = first(values(tpsths[1]["cmdepthpsth"]))
psthtimes=ffpsth.x;bw = psthtimes[2]-psthtimes[1];psthdepths = ffpsth.y;psthdn=ffpsth.n

psthb = [0 15]
psthr = [15 100]
psthbi = epoch2sampleindex(psthb,1/(bw*SecondPerUnit))
psthri = epoch2sampleindex(psthr,1/(bw*SecondPerUnit))
psthx = psthtimes[psthri]
psths = mapreduce((l,i)->["$(l)_$k"=>imfilter(stfilter(v.psth,temporaltype=:sub,ti=psthbi),Kernel.gaussian((1,1)))[:,psthri] for (k,v) in i["cmdepthpsth"]],
append!,testlogs,tpsths)
psthclim=absmax(psths)

## Set Layers
plotlayer=(o...;w=230,h=510)->begin
    cn = length(csds)
    p=plot(layout=(3,cn),link=:y,legend=false,grid=false,size=(cn*w,3h))
    for i in 1:cn
        yticks = i==1 ? true : false
        xlabel = i==1 ? "Time (ms)" : ""
        ylabel = i==1 ? "Depth (um)" : ""
        heatmap!(p,subplot=i,csdx,lfpdepths,csds[i].second,color=:RdBu,clims=(-csdclim,csdclim),title=csds[i].first,titlefontsize=6,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
        if !isnothing(layer)
            hline!(p,subplot=i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
            annotations=[(csdx[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
        end
    end
    for i in 1:cn
        yticks = i==1 ? true : false
        xlabel = i==1 ? "Time (ms)" : ""
        ylabel = i==1 ? "Depth (um)" : ""
        heatmap!(p,subplot=cn+i,psthx,psthdepths,psths[i].second,color=:coolwarm,clims=(-psthclim,psthclim),title=psths[i].first,titlefontsize=6,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
        if i==1
            pn = maximum(psthx) .- psthdn./maximum(psthdn) .* maximum(psthx) .* 0.2
            plot!(p,subplot=cn+i,pn,psthdepths,label="Number of Units",color=:seagreen,lw=0.5)
        end
        if !isnothing(layer)
            hline!(p,subplot=cn+i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
            annotations=[(psthx[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
        end
    end
    for i in 1:cn
        yticks = i==1 ? true : false
        xlabel = i==1 ? "Frequency (Hz)" : ""
        ylabel = i==1 ? "Depth (um)" : ""
        heatmap!(p,subplot=2cn+i,freq,lfpdepths,pcs[i].second,color=:vik,clims=(-pcclim,pcclim),title=pcs[i].first,titlefontsize=6,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
        if !isnothing(layer)
            hline!(p,subplot=2cn+i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
            annotations=[(freq[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
        end
    end
    p
end

lw = Dict(k=>widget(0:3800,label=k,value=layer[k][1]) for k in keys(layer))
foreach(k->on(v->layer[k][1]=v,lw[k]),keys(lw))
lp = map(plotlayer,values(lw)...)
vbox(values(lw)...,lp)

plotlayer()
foreach(ext->savefig(joinpath(siteresultdir,"Layer_dCSD_dPSTH_PowerContrast$ext")),figfmt)

## Finalize Layer
save(joinpath(siteresultdir,"layer.jld2"),"layer",checklayer!(layer),"siteid",siteid)

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
