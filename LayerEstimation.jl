using NeuroAnalysis,Statistics,FileIO,JLD2,Plots,Images,Interact,ProgressMeter,XLSX

## Combine all layer tests of one RecordSite
dataroot = "X:/"
dataexportroot = "Y:/"
resultroot = "Z:/"

subject = "AG2";recordsession = "V1";recordsite = "ODR1"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)
figfmt = [".png"]


# df=DataFrame(power=[405.5,239.0,92.8,410.5,115.2,493.6],hue=repeat(["180","0"],outer=3),location=repeat(["P5","P4","P3"],inner=2))
#
# df |> @vlplot(:bar,x={:hue,axis={labelAngle=0}},y={:power,axis={grid=false}},column=:location,color={:hue,scale={range=["#FF5A81","#00A57E"]}})
#
#
# # P5
# mm=405.5
# lm=239.0
# # P4
# mm=92.8
# lm=410.5
# # P3
# mm = 115.2
# lm = 493.6


# lfp=load(joinpath(siteresultdir,"$(siteid)_Flash2Color_2","lfp.jld2"))
# freq=lfp["freq"]
# depths = -lfp["depth"].+layer["Out"][1]
# pc = lfp["cmpc"]
# mpc = pc["Color=[0.0, 0.6462, 0.4939, 1.0]"]
# lpc = pc["Color=[1.0, 0.3538, 0.5061, 1.0]"]
# li = 0 .<=depths.<=300
# fi = 30 .<=freq.<=100
#
# mm = sum(mpc[li,fi])
# lm = sum(lpc[li,fi])
#
#
#
# lfp=load(joinpath(siteresultdir,"$(siteid)_Flash2Color_3","lfp.jld2"))
# freq=lfp["freq"]
# depths = -lfp["depth"].+layer["Out"][1]
# pc = lfp["cmpc"]
# lmpc = pc["Color=[0.3672, 0.6168, 0.0, 1.0]"]
# spc = pc["Color=[0.6328, 0.3832, 1.0, 1.0]"]
# li = 100 .<=depths.<=400
# fi = 30 .<=freq.<=100
#
# lmm = mean(lmpc[li,fi])
# sm = mean(spc[li,fi])

# testids = ["$(siteid)_00$i" for i in 0:3]
# testn=length(testids)
# testtitles = ["Eyeₚ","Eyeₙₚ","Eyes","EyeₚS"]
# csds = load.(joinpath.(sitedir,testids,"csd.jld2"),"csd","depth","fs")
# depths=csds[1][2];fs=csds[1][3]


## Flash2Color
flashdepth = (siteid,siteresultdir) -> begin

ii = '0'
test = "Flash2Color"
testids = ["$(siteid)_$(test)_$i" for i in 0:3]
testn=length(testids)
# AP dRMS
aps = load.(joinpath.(siteresultdir,testids,"ap$ii.jld2"))
titles = repeat(["$(i["eye"])_$(i["color"])" for i in aps],inner=2)
colors = mapreduce(i->[RGBA(parse.(Float64,split(match(r".*Color=\[(.*)\]",k).captures[1],", "))...) for k in keys(i["cdrms"])],append!,aps)
drms = mapreduce(i->[1e6v for v in values(i["cdrms"])],append!,aps)
time = mapreduce(i->[i["time"],i["time"]],append!,aps)
depths = mapreduce(i->[i["depths"],i["depths"]],append!,aps)
rmslim=absmax(drms)

plotflashdepthresponse(drms,rmslim,time,depths,colors,titles)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_dRMS$ext")),figfmt)

# Unit dPSTH
dpsths = load.(joinpath.(siteresultdir,testids,"depthpsth$ii.jld2"))
dpsth = mapreduce(i->[v.psth for v in values(i["cdpsth"])],append!,dpsths)
time = mapreduce(i->[v.x for v in values(i["cdpsth"])],append!,dpsths)
depths = mapreduce(i->[v.y for v in values(i["cdpsth"])],append!,dpsths)
psthlim=absmax(dpsth)

plotflashdepthresponse(dpsth,psthlim,time,depths,colors,titles)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_dPSTH$ext")),figfmt)

# LFP and CSD
lfps = load.(joinpath.(siteresultdir,testids,"lfp$ii.jld2"))
lfp = mapreduce(i->[1e6v for v in values(i["cmlfp"])],append!,lfps)
dcsd = mapreduce(i->[v for v in values(i["cdcsd"])],append!,lfps)
time = mapreduce(i->[i["time"],i["time"]],append!,lfps)
depths = mapreduce(i->[i["depths"],i["depths"]],append!,lfps)
lfplim=absmax(lfp)
csdlim=absmax(dcsd)

plotflashdepthresponse(lfp,lfplim,time,depths,colors,titles)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_LFP$ext")),figfmt)

plotflashdepthresponse(dcsd,csdlim,time,depths,colors,titles,color=:RdBu)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_dCSD$ext")),figfmt)

# 5x1 gaussian kernal(80μm)
fdcsd = mapreduce(i->[imfilter(v,Kernel.gaussian((1,0)),Fill(0)) for v in values(i["cdcsd"])],append!,lfps)
plotflashdepthresponse(fdcsd,absmax(fdcsd),time,depths,colors,titles,color=:RdBu)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_fdCSD$ext")),figfmt)
end


absmax(x) = mapreduce(i->maximum(abs.(i)),max,x)
plotflashdepthresponse=(resp,lim,time,depths,colors,titles;w=140,h=600,color=:coolwarm)->begin
    n = length(resp)
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:500:depths[i][end]) : false
        xticks = (0:50:time[i][end])
        lmargin = i==1 ? 4mm : -4mm
        xlabel = i==1 ? "Time (ms)" : ""
        ylabel = i==1 ? "Depth (μm)" : ""

        heatmap!(p[i],time[i],depths[i],resp[i];color=color,clims=(-lim,lim),
        title=titles[i],titlefontsize=10,yticks,xticks,tickor=:out,xlabel,ylabel,
        annotation=[(5,50,Plots.text("■",15,colors[i],:left,:bottom))],
        left_margin=lmargin,bottom_margin=3mm)
    end
    p
end


## OriSF
orisfdepth = (siteid,siteresultdir) -> begin

ii = '0'
test = "OriSF"
testids = ["$(siteid)_$(test)_$i" for i in 0:3]
testn=length(testids)
# AP dRMS
aps = load.(joinpath.(siteresultdir,testids,"ap$ii.jld2"))
titles = repeat(["$(i["eye"])_$(i["color"])" for i in aps],inner=2)
colors = mapreduce(i->[RGBA(parse.(Float64,split(match(r".*Color=\[(.*)\]",k).captures[1],", "))...) for k in keys(i["cdrms"])],append!,aps)
drms = mapreduce(i->[1e6v for v in values(i["cdrms"])],append!,aps)
time = mapreduce(i->[i["time"],i["time"]],append!,aps)
depths = mapreduce(i->[i["depths"],i["depths"]],append!,aps)
rmslim=absmax(drms)

cmdrms = map(i->mapreduce(v->1e6v,.+, values(i["cdrms"])),aps)

p=plot(layout=(1,4),link=:y,legend=false,grid=false)
for i in 1:4
    yticks = i==1 ? (0:500:depths[i][end]) : false
    xticks = (0:50:time[i][end])
    lmargin = i==1 ? 4mm : -4mm
    xlabel = i==1 ? "Time (ms)" : ""
    ylabel = i==1 ? "Depth (μm)" : ""

    lim = maximum(abs.(cmdrms[i]))
    clims = (-lim,lim)
    heatmap!(p[i],time[i],depths[i],cmdrms[i];yticks,xticks,color=:coolwarm,clims)
    anno = [(time[i][1]+3,mean(layer[k]),text(k,6,:gray10,:center,:left)) for k in keys(layer)]
    hline!(p[i],[l[2] for l in values(layer)],linecolor=:gray25,legend=false,lw=0.5,annotations=anno)
end

p


plotflashdepthresponse(drms,rmslim,time,depths,colors,titles)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_dRMS$ext")),figfmt)

# Unit dPSTH
dpsths = load.(joinpath.(siteresultdir,testids,"depthpsth$ii.jld2"))
dpsth = mapreduce(i->[v.psth for v in values(i["cdpsth"])],append!,dpsths)
time = mapreduce(i->[v.x for v in values(i["cdpsth"])],append!,dpsths)
depths = mapreduce(i->[v.y for v in values(i["cdpsth"])],append!,dpsths)
psthlim=absmax(dpsth)

plotflashdepthresponse(dpsth,psthlim,time,depths,colors,titles)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_dPSTH$ext")),figfmt)

# LFP and CSD
lfps = load.(joinpath.(siteresultdir,testids,"lfp$ii.jld2"))
lfp = mapreduce(i->[1e6v for v in values(i["cmlfp"])],append!,lfps)
dcsd = mapreduce(i->[v for v in values(i["cdcsd"])],append!,lfps)
time = mapreduce(i->[i["time"],i["time"]],append!,lfps)
depths = mapreduce(i->[i["depths"],i["depths"]],append!,lfps)
lfplim=absmax(lfp)
csdlim=absmax(dcsd)

plotflashdepthresponse(lfp,lfplim,time,depths,colors,titles)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_LFP$ext")),figfmt)

plotflashdepthresponse(dcsd,csdlim,time,depths,colors,titles,color=:RdBu)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_dCSD$ext")),figfmt)

# 5x1 gaussian kernal(80μm)
fdcsd = mapreduce(i->[imfilter(v,Kernel.gaussian((1,0)),Fill(0)) for v in values(i["cdcsd"])],append!,lfps)
plotflashdepthresponse(fdcsd,absmax(fdcsd),time,depths,colors,titles,color=:RdBu)
foreach(ext->savefig(joinpath(siteresultdir,"$(test)_fdCSD$ext")),figfmt)
end










## Batch Penetration Sites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1")...)
@showprogress "Batch FlashDepth ... " for r in eachrow(penetration)
    flashdepth(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid))
    orisfdepth(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid))
end


dcsd = mapreduce(i->[mapwindow(mean,v,(11,11),border=Fill(0)) for v in values(i["cdcsd"])],append!,lfps)
dcsd = mapreduce(i->[imfilter(v,Kernel.gaussian((1,1)),Fill(0)) for v in values(i["cdcsd"])],append!,lfps)

dcsd = mapreduce(i->[csd(v[1:2:end,:],h=hy) for v in values(i["cmlfp"])],append!,lfps)
dcsd = mapreduce(i->[imfilter(csd(v[1:2:end,:],h=hy),Kernel.gaussian((1,0)),Fill(0)) for v in values(i["cmlfp"])],append!,lfps)
depths = mapreduce(i->[i["depths"][1:2:end],i["depths"][1:2:end]],append!,lfps)







## Set Layers
# plotlayer=(o...;w=230,h=510)->begin
#     cn = length(csds)
#     p=plot(layout=(3,cn),link=:y,legend=false,grid=false,size=(cn*w,3h))
#     for i in 1:cn
#         yticks = i==1 ? true : false
#         xlabel = i==1 ? "Time (ms)" : ""
#         ylabel = i==1 ? "Depth (um)" : ""
#         heatmap!(p,subplot=i,csdx,lfpdepths,csds[i].second,color=:RdBu,clims=(-csdclim,csdclim),title=csds[i].first,titlefontsize=6,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
#         if !isnothing(layer)
#             hline!(p,subplot=i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
#             annotations=[(csdx[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
#         end
#     end
#     for i in 1:cn
#         yticks = i==1 ? true : false
#         xlabel = i==1 ? "Time (ms)" : ""
#         ylabel = i==1 ? "Depth (um)" : ""
#         heatmap!(p,subplot=cn+i,psthx,psthdepths,psths[i].second,color=:coolwarm,clims=(-psthclim,psthclim),title=psths[i].first,titlefontsize=6,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
#         if i==1
#             pn = maximum(psthx) .- psthdn./maximum(psthdn) .* maximum(psthx) .* 0.2
#             plot!(p,subplot=cn+i,pn,psthdepths,label="Number of Units",color=:seagreen,lw=0.5)
#         end
#         if !isnothing(layer)
#             hline!(p,subplot=cn+i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
#             annotations=[(psthx[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
#         end
#     end
#     for i in 1:cn
#         yticks = i==1 ? true : false
#         xlabel = i==1 ? "Frequency (Hz)" : ""
#         ylabel = i==1 ? "Depth (um)" : ""
#         heatmap!(p,subplot=2cn+i,freq,lfpdepths,pcs[i].second,color=:vik,clims=(-pcclim,pcclim),title=pcs[i].first,titlefontsize=6,yticks=yticks,xlabel=xlabel,ylabel=ylabel)
#         if !isnothing(layer)
#             hline!(p,subplot=2cn+i,[l[1] for l in values(layer)],linestyle=:dash,linecolor=:gray25,legend=false,
#             annotations=[(freq[1]+1,layer[k][1],text(k,4,:gray10,:bottom,:left)) for k in keys(layer)])
#         end
#     end
#     p
# end
#
# lw = Dict(k=>widget(0:3800,label=k,value=layer[k][1]) for k in keys(layer))
# foreach(k->on(v->layer[k][1]=v,lw[k]),keys(lw))
# lp = map(plotlayer,values(lw)...)
# vbox(values(lw)...,lp)
#
# plotlayer()
# foreach(ext->savefig(joinpath(siteresultdir,"Layer_dCSD_dPSTH_PowerContrast$ext")),figfmt)


# earliest response should be due to LGN M,P input to 4Ca,4Cb
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
