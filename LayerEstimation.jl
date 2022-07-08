using NeuroAnalysis,Statistics,StatsBase,FileIO,JLD2,StatsPlots,Images,ProgressMeter,DataFrames,XLSX,Dierckx

## Combine all layer tests of one RecordSite
dataroot = "X:/"
dataexportroot = "Y:/"
resultroot = "Z:/"
figfmt = [".png"]


# Neuropixels 1.0 layout
xs = repeat([16,48,0,32],outer=96)
ys = repeat(range(0,step=20,length=192),inner=2)

deleteat!(xs,192)
deleteat!(ys,192)

scatter(xs,ys,m=:rect,ms=8,msw=0,color=:gray35,size=(300,3840),leg=false,frame=:none,ratio=:equal)
scatter!([32],[1900],m=:rect,ms=8,msw=0,color=HSL(0,0.9,0.6),leg=false)


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

cmcode = Dict("DKL_HueL0"=>ColorMaps["lidkl_mcchue_l0"],"HSL_HueYm"=>ColorMaps["hsl_mshue_l0.4"])
absmax(x) = mapreduce(i->maximum(abs.(i)),max,x)
abspct(x,p=99.9) = percentile(vec(abs.(x)),p)
plotlayerunitfeature=(udf,depths;w=140,h=600,color=:tab10,layer=nothing)->begin
    n = 4;cs = palette(color).colors.colors[1:n]
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:200:depths[end]) : false
        leftmargin = i==1 ? 4mm : -4mm
        bottommargin = 3mm
        ylabel = i==1 ? "Depth (μm)" : ""

        if i == 1
            x = udf.Density
            xmin,xmax = extrema(x)
            plot!(p[i],x,depths;xlabel="Density",ylabel,leg=false,color=cs[i],
            lw=1.5,yticks,xticks=[ceil(xmin,sigdigits=1),floor(xmax,sigdigits=1)],
            leftmargin,bottommargin,tickor=:out)
        elseif i == 2
            x = [udf.upspread;reverse(udf.downspread)]
            y = [depths;reverse(depths)]
            plot!(p[i],x,y;st=:shape,xlabel="Spread",lw=0,yticks,alpha=0.8,xticks=[0],color=cs[i],
            leftmargin,bottommargin,tickor=:out)
        elseif i == 3
            x = [0;udf.peaktroughratio;0]
            y = [depths[1];depths;depths[end]]
            plot!(p[i],x,y;st=:shape,xlabel="PTRatio",lw=0,yticks,alpha=0.8,xticks=[0],color=cs[i],
            leftmargin,bottommargin,tickor=:out)
        elseif i == 4
            x = [0;udf.duration;0]
            y = [depths[1];depths;depths[end]]
            plot!(p[i],x,y;st=:shape,xlabel="Duration",lw=0,yticks,alpha=0.8,xticks=[0],color=cs[i],
            leftmargin,bottommargin,tickor=:out)
        end

        xmin,xmax = extrema(x)
        if !isnothing(layer)
            ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
            hline!(p[i],[l[2] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ann)
        end
    end
    p
end
plotcondresponse=(resp,times,depths;colors=[],minmaxcolors=[],colormaps=[],titles="",titlefontsize=7,w=140,h=600,color=:coolwarm,layer=nothing)->begin
    n = length(resp);lim = absmax(resp)
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        isnothing(layer) || (lim=abspct(resp[i]))
        yticks = i==1 ? (0:200:depths[i][end]) : false
        xticks = (0:50:times[i][end])
        leftmargin = i==1 ? 4mm : -4mm
        xlabel = i==1 ? "Time (ms)" : ""
        ylabel = i==1 ? "Depth (μm)" : ""
        title = isempty(titles) ? titles : titles[i]
        xmin,xmax = extrema(times[i])
        ann = isempty(colors) ? [] : [(5,50,Plots.text("■",15,colors[i],:left,:bottom))]
        ann = isempty(minmaxcolors) ? ann : [(5,50,Plots.text("▮",15,minmaxcolors[i][2],:left,:bottom)),(17,50,Plots.text("▮",15,minmaxcolors[i][1],:left,:bottom))]
        ann = isempty(colormaps) ? ann : [(5+7.5(j-1),50,Plots.text("▮",8,colormaps[i][j],:left,:bottom)) for j in eachindex(colormaps[i])]

        heatmap!(p[i],times[i],depths[i],resp[i];color,clims=(-lim,lim),
        title,titlefontsize,yticks,xticks,tickor=:out,xlabel,ylabel,
        ann,leftmargin,bottommargin=3mm)

        if !isnothing(layer)
            ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
            hline!(p[i],[l[2] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ann)
            # yticks = i==1 ? [layer["1"][2],layer["WM"][2]] : false
            # yformatter = x->Int(layer["1"][2]-x)
            # hline!(p[i],[l[2] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ann,yticks,yformatter)
        end
    end
    p
end
plottrialresponse=(resp,trials,depths;colors=[],minmaxcolors=[],colormaps=[],titles="",titlefontsize=7,w=140,h=600,color=:heat,layer=nothing)->begin
    n = length(resp);lim = absmax(resp)
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        isnothing(layer) || (lim=abspct(resp[i]))
        yticks = i==1 ? (0:200:depths[i][end]) : false
        xticks = range(start=1,stop=round(Int,0.8*trials[i][end]),length=2)
        leftmargin = i==1 ? 4mm : -4mm
        xlabel = i==1 ? "Trial" : ""
        ylabel = i==1 ? "Depth (μm)" : ""
        title = isempty(titles) ? titles : titles[i]
        xmin,xmax = extrema(trials[i])
        ann = isempty(colors) ? [] : [(30,50,Plots.text("▮",15,colors[1+2(i-1)],:left,:bottom)),(60,50,Plots.text("▮",15,colors[2+2(i-1)],:left,:bottom))]
        ann = isempty(minmaxcolors) ? ann : [(30,50,Plots.text("▮",15,minmaxcolors[i][2],:left,:bottom)),(90,50,Plots.text("▮",15,minmaxcolors[i][1],:left,:bottom))]
        ann = isempty(colormaps) ? ann : [(0.1*xmax+0.045*xmax*(j-1),50,Plots.text("▮",8,colormaps[i][j],:left,:bottom)) for j in eachindex(colormaps[i])]

        heatmap!(p[i],trials[i],depths[i],resp[i];color,clims=(0,lim),
        title,titlefontsize,yticks,xticks,tickor=:out,xlabel,ylabel,
        ann,leftmargin,bottommargin=3mm)

        n = dropdims(mean(resp[i],dims=2),dims=2)
        pn = clampscale(n) .* (xmax-xmin) .* 0.3 .+ xmin
        plot!(p[i],pn,depths[i],color=:gray40,leg=false,xlim=(xmin,xmax),ylim=extrema(depths[i]))
        if !isnothing(layer)
            ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
            hline!(p[i],[l[2] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ann)
            # yticks = i==1 ? [layer["1"][2],layer["WM"][2]] : false
            # yformatter = x->Int(layer["1"][2]-x)
            # hline!(p[i],[l[2] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ann,yticks,yformatter)
        end
    end
    p
end



## Manual Interactive Layer Assignment
subject = "AG2";recordsession = "V1";recordsite = "ODR2"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)

layer = load(joinpath(siteresultdir,"layer.jld2"),"layer")
# layer = Dict()

layer["1"] = [3100,3160]
layer["2"] = [3000,3000]
layer["3"] = [2650,2830]
layer["4A"] = [2700,2620]
layer["4B"] = [2350,2560]
layer["4Cα"] = [2250,2395]
layer["4Cβ"] = [2050,2290]
layer["5"] = [1800,2180]
layer["6"] = [1600,1800]
layer["WM"] = [0,1625]

# Finalize Layer
layer = checklayer!(layer)
jldsave(joinpath(siteresultdir,"layer.jld2");layer,siteid)

flashdepth(siteid,siteresultdir;layer)
colordepth(siteid,siteresultdir;layer)
orisfdepth(siteid,siteresultdir;layer)



## Flash2Color
flashdepth = (siteid,siteresultdir;ii='0',layer=nothing) -> begin

test = "Flash2Color"
testids = ["$(siteid)_$(test)_$i" for i in 0:3]
testn=length(testids)
if layer==:batch
    layerpath = joinpath(siteresultdir,"layer.jld2")
    layer = ispath(layerpath) ? load(layerpath,"layer") : nothing
end
# Unit
if !isnothing(layer)
    df = load(joinpath(siteresultdir,"unitdepthfeature.jld2"),"df")
    udf = (;(Symbol(k)=>df["depthfeature"][:,k.==df["feature"][:]] for k in df["feature"])...)
    plotlayerunitfeature(udf,df["depth"];layer)
    foreach(ext->savefig(joinpath(siteresultdir,"Layer_UnitFeature$ext")),figfmt)
end
# AP
aps = load.(joinpath.(siteresultdir,testids,"ap$ii.jld2"))
titles = repeat(["$(i["exenv"]["eye"])_$(i["exenv"]["color"])" for i in aps],inner=2)
colors = mapreduce(i->[RGBA(parse.(Float64,split(match(r".*Color=\[(.*)\]",k).captures[1],", "))...) for k in keys(i["crms"])],append!,aps)
crms = mapreduce(i->[v for v in values(i["crms"])],append!,aps)
wcrms = mapreduce(i->[v for v in values(i["wcrms"])],append!,aps)
times = mapreduce(i->[i["times"],i["times"]],append!,aps)
depths = mapreduce(i->[i["depths"],i["depths"]],append!,aps)

plotcondresponse(crms,times,depths;colors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dRMS$ext")),figfmt)
plotcondresponse(wcrms,times,depths;colors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wdRMS$ext")),figfmt)

prms = map(i->i["prms"],aps)
wprms = map(i->i["wprms"],aps)
trials = map(i->1:size(i,2),prms)
plottrialresponse(prms,trials,depths;colors,titles=titles[1:2:end],layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_RMS$ext")),figfmt)
plottrialresponse(wprms,trials,depths;colors,titles=titles[1:2:end],layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wRMS$ext")),figfmt)

pc = map(i->i["pc"],aps)
wpc = map(i->i["wpc"],aps)
plottrialresponse(pc,trials,depths;colors,titles=titles[1:2:end],layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_Coherence$ext")),figfmt)
plottrialresponse(wpc,trials,depths;colors,titles=titles[1:2:end],layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wCoherence$ext")),figfmt)

pdc = map(i->abs.(i["pc"].-i["pbc"]),aps)
wpdc = map(i->abs.(i["wpc"].-i["wpbc"]),aps)
plottrialresponse(pdc,trials,depths;colors,titles=titles[1:2:end],layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCoherence$ext")),figfmt)
plottrialresponse(wpdc,trials,depths;colors,titles=titles[1:2:end],layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wdCoherence$ext")),figfmt)
jldsave(joinpath(siteresultdir,"$(test)_ap.jld2");crms,pc,pdc,times,trials,depths,colors,titles,siteid)

# Unit PSTH
cpsths = load.(joinpath.(siteresultdir,testids,"depthpsth$ii.jld2"))
cpsth = mapreduce(i->[v.psth for v in values(i["cpsth"])],append!,cpsths)
xs = mapreduce(i->[v.x for v in values(i["cpsth"])],append!,cpsths)
ys = mapreduce(i->[v.y for v in values(i["cpsth"])],append!,cpsths)

plotcondresponse(cpsth,xs,ys;colors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dPSTH$ext")),figfmt)

# LFP and CSD
lfs = load.(joinpath.(siteresultdir,testids,"lf$ii.jld2"))
clfp = mapreduce(i->[1e6v for v in values(i["clfp"])],append!,lfs)
ccsd = mapreduce(i->[v for v in values(i["ccsd"])],append!,lfs)
times = mapreduce(i->[i["times"],i["times"]],append!,lfs)
depths = mapreduce(i->[i["depths"],i["depths"]],append!,lfs)
cm = cgrad(:jet,rev=true)

# heatmap(rand(20,20),color=c,size=(400,600))

plotcondresponse(clfp,times,depths;colors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_LFP$ext")),figfmt)
plotcondresponse(ccsd,times,depths;colors,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSD$ext")),figfmt)

# σ=1(20μm): 5(80μm) diameter gaussian kernal
fccsd = mapreduce(i->[imfilter(v,Kernel.gaussian((1,0)),Fill(0)) for v in values(i["ccsd"])],append!,lfs)
plotcondresponse(fccsd,times,depths;colors,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf1$ext")),figfmt)

# σ=2(40μm): 9(160μm) diameter gaussian kernal
fccsd = mapreduce(i->[imfilter(v,Kernel.gaussian((2,0)),Fill(0)) for v in values(i["ccsd"])],append!,lfs)
plotcondresponse(fccsd,times,depths;colors,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf2$ext")),figfmt)
jldsave(joinpath(siteresultdir,"$(test)_lf.jld2");fccsd,times,depths,colors,titles,siteid)

# σ=3(60μm): 13(240μm) diameter gaussian kernal
fccsd = mapreduce(i->[imfilter(v,Kernel.gaussian((3,0)),Fill(0)) for v in values(i["ccsd"])],append!,lfs)
plotcondresponse(fccsd,times,depths;colors,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf3$ext")),figfmt)
end


## Color
colordepth = (siteid,siteresultdir;ii='0',layer=nothing) -> begin

test = "Color"
testids = ["$(siteid)_$(test)_$i" for i in 0:1]
testn=length(testids)
if layer==:batch
    layerpath = joinpath(siteresultdir,"layer.jld2")
    layer = ispath(layerpath) ? load(layerpath,"layer") : nothing
end
# AP
aps = load.(joinpath.(siteresultdir,testids,"ap$ii.jld2"))
titles = ["$(i["exenv"]["eye"])_$(i["exenv"]["color"])" for i in aps]
colormaps = map(i->cmcode[i["exenv"]["color"]].colors[range(1,step=45,length=8)],aps)
cmrms = map(i->reduce(.+, values(i["crms"])),aps)
wcmrms = map(i->reduce(.+, values(i["wcrms"])),aps)
times = map(i->i["times"],aps)
depths = map(i->i["depths"],aps)

plotcondresponse(cmrms,times,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dRMS$ext")),figfmt)
plotcondresponse(wcmrms,times,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wdRMS$ext")),figfmt)

prms = map(i->i["prms"],aps)
wprms = map(i->i["wprms"],aps)
trials = map(i->1:size(i,2),prms)
plottrialresponse(prms,trials,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_RMS$ext")),figfmt)
plottrialresponse(wprms,trials,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wRMS$ext")),figfmt)

pc = map(i->i["pc"],aps)
wpc = map(i->i["wpc"],aps)
plottrialresponse(pc,trials,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_Coherence$ext")),figfmt)
plottrialresponse(wpc,trials,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wCoherence$ext")),figfmt)

pdc = map(i->abs.(i["pc"].-i["pbc"]),aps)
wpdc = map(i->abs.(i["wpc"].-i["wpbc"]),aps)
plottrialresponse(pdc,trials,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCoherence$ext")),figfmt)
plottrialresponse(wpdc,trials,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wdCoherence$ext")),figfmt)

# Unit PSTH
cpsths = load.(joinpath.(siteresultdir,testids,"depthpsth$ii.jld2"))
cmpsth = map(i->mapreduce(v->v.psth ,.+, values(i["cpsth"])),cpsths)
xs = map(i->first(values(i["cpsth"])).x,cpsths)
ys = map(i->first(values(i["cpsth"])).y,cpsths)

plotcondresponse(cmpsth,xs,ys;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dPSTH$ext")),figfmt)

# LFP and CSD
lfs = load.(joinpath.(siteresultdir,testids,"lf$ii.jld2"))
cmlfp = map(i->mapreduce(j->1e6j,.+, values(i["clfp"])),lfs)
cmcsd = map(i->reduce(.+, values(i["ccsd"])),lfs)
times = map(i->i["times"],lfs)
depths = map(i->i["depths"],lfs)
cm = cgrad(:jet,rev=true)

plotcondresponse(cmlfp,times,depths;colormaps,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_LFP$ext")),figfmt)
plotcondresponse(cmcsd,times,depths;colormaps,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSD$ext")),figfmt)

# σ=1(20μm): 5(80μm) diameter gaussian kernal
fcmcsd = map(i->mapreduce(j->imfilter(j,Kernel.gaussian((1,0)),Fill(0)),.+,values(i["ccsd"])),lfs)
plotcondresponse(fcmcsd,times,depths;colormaps,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf1$ext")),figfmt)

# σ=2(40μm): 9(160μm) diameter gaussian kernal
fcmcsd = map(i->mapreduce(j->imfilter(j,Kernel.gaussian((2,0)),Fill(0)),.+,values(i["ccsd"])),lfs)
plotcondresponse(fcmcsd,times,depths;colormaps,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf2$ext")),figfmt)

# σ=3(60μm): 13(240μm) diameter gaussian kernal
fcmcsd = map(i->mapreduce(j->imfilter(j,Kernel.gaussian((3,0)),Fill(0)),.+,values(i["ccsd"])),lfs)
plotcondresponse(fcmcsd,times,depths;colormaps,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf3$ext")),figfmt)
end


## OriSF
orisfdepth = (siteid,siteresultdir;ii='0',layer=nothing) -> begin

test = "OriSF"
testids = ["$(siteid)_$(test)_$i" for i in 0:3]
testn=length(testids)
if layer==:batch
    layerpath = joinpath(siteresultdir,"layer.jld2")
    layer = ispath(layerpath) ? load(layerpath,"layer") : nothing
end
# AP
aps = load.(joinpath.(siteresultdir,testids,"ap$ii.jld2"))
titles = ["$(i["exenv"]["eye"])_$(i["exenv"]["color"])" for i in aps]
minmaxcolors = [i["exenv"]["minmaxcolor"] for i in aps]
cmrms = map(i->reduce(.+, values(i["crms"])),aps)
wcmrms = map(i->reduce(.+, values(i["wcrms"])),aps)
times = map(i->i["times"],aps)
depths = map(i->i["depths"],aps)

plotcondresponse(cmrms,times,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dRMS$ext")),figfmt)
plotcondresponse(wcmrms,times,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wdRMS$ext")),figfmt)

prms = map(i->i["prms"],aps)
wprms = map(i->i["wprms"],aps)
trials = map(i->1:size(i,2),prms)
plottrialresponse(prms,trials,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_RMS$ext")),figfmt)
plottrialresponse(wprms,trials,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wRMS$ext")),figfmt)

pc = map(i->i["pc"],aps)
wpc = map(i->i["wpc"],aps)
plottrialresponse(pc,trials,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_Coherence$ext")),figfmt)
plottrialresponse(wpc,trials,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wCoherence$ext")),figfmt)

pdc = map(i->abs.(i["pc"].-i["pbc"]),aps)
wpdc = map(i->abs.(i["wpc"].-i["wpbc"]),aps)
plottrialresponse(pdc,trials,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCoherence$ext")),figfmt)
plottrialresponse(wpdc,trials,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_wdCoherence$ext")),figfmt)

# Unit PSTH
cpsths = load.(joinpath.(siteresultdir,testids,"depthpsth$ii.jld2"))
cmpsth = map(i->mapreduce(v->v.psth ,.+, values(i["cpsth"])),cpsths)
xs = map(i->first(values(i["cpsth"])).x,cpsths)
ys = map(i->first(values(i["cpsth"])).y,cpsths)

plotcondresponse(cmpsth,xs,ys;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dPSTH$ext")),figfmt)

# LFP and CSD
lfs = load.(joinpath.(siteresultdir,testids,"lf$ii.jld2"))
cmlfp = map(i->mapreduce(j->1e6j,.+, values(i["clfp"])),lfs)
cmcsd = map(i->reduce(.+, values(i["ccsd"])),lfs)
times = map(i->i["times"],lfs)
depths = map(i->i["depths"],lfs)
cm = cgrad(:jet,rev=true)

plotcondresponse(cmlfp,times,depths;minmaxcolors,titles,layer)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_LFP$ext")),figfmt)
plotcondresponse(cmcsd,times,depths;minmaxcolors,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSD$ext")),figfmt)

# σ=1(20μm): 5(80μm) diameter gaussian kernal
fcmcsd = map(i->mapreduce(j->imfilter(j,Kernel.gaussian((1,0)),Fill(0)),.+,values(i["ccsd"])),lfs)
plotcondresponse(fcmcsd,times,depths;minmaxcolors,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf1$ext")),figfmt)

# σ=2(40μm): 9(160μm) diameter gaussian kernal
fcmcsd = map(i->mapreduce(j->imfilter(j,Kernel.gaussian((2,0)),Fill(0)),.+,values(i["ccsd"])),lfs)
plotcondresponse(fcmcsd,times,depths;minmaxcolors,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf2$ext")),figfmt)

# σ=3(60μm): 13(240μm) diameter gaussian kernal
fcmcsd = map(i->mapreduce(j->imfilter(j,Kernel.gaussian((3,0)),Fill(0)),.+,values(i["ccsd"])),lfs)
plotcondresponse(fcmcsd,times,depths;minmaxcolors,titles,layer,color=cm)
foreach(ext->savefig(joinpath(siteresultdir,"$(isnothing(layer) ? "" : "Layer_")$(test)_dCSDf3$ext")),figfmt)
end



## Batch Penetration Sites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1")...)
@showprogress "Batch Layer ... " for r in eachrow(penetration)
    # flashdepth(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=nothing)
    # colordepth(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=nothing)
    # orisfdepth(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=nothing)
    flashdepth(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=:batch)
    # colordepth(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=:batch)
    # orisfdepth(r.siteid,joinpath(resultroot,r.Subject_ID,r.siteid),layer=:batch)
end



dcsd = mapreduce(i->[csd(v[1:2:end,:],h=hy) for v in values(i["cmlfp"])],append!,lfps)
dcsd = mapreduce(i->[imfilter(csd(v[1:2:end,:],h=hy),Kernel.gaussian((1,0)),Fill(0)) for v in values(i["cmlfp"])],append!,lfps)
depths = mapreduce(i->[i["depths"][1:2:end],i["depths"][1:2:end]],append!,lfps)







ln = ["1", "2", "3", "4A", "4B", "4Cα", "4Cβ", "5", "6", "WM"]


lwd = []

p2b(p) = cumsum([p[1];-1*p[2:end]])

ofun(p;λ=10,data=data,lwd=lwd) = bfun(p2b(p),data) + λ*mfun(p[2:end],lwd)

mfun(lw,lwd) = mapreduce((w,d)->d(w),+,lw,lwd)

function bfun(b,data)
    body
end





## Collect Layers of All RecordSites
function collectlayer!(indir;layer=DataFrame(),datafile="layer.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            siteid,l = load(joinpath(root,datafile),"siteid","layer")
            df = DataFrame("siteid"=>siteid,(k=>[l[k]] for k in sort(collect(keys(l))))...)
            append!(layer,df,cols=:union)
        end
    end
    return layer
end

alllayer = collectlayer!(resultroot)
save(joinpath(resultroot,"alllayer.jld2"),"alllayer",alllayer)

## Layer Statistics and Transformation
using KernelDensity
alllayer = load(joinpath(resultroot,"alllayer.jld2"),"alllayer")
transform!(alllayer,"1"=>ByRow(i->(x->last(i).-x))=>"e2cfun")
transform!(alllayer,vcat.("e2cfun",names(alllayer,Not([:siteid,:e2cfun]))) .=> ByRow((f,x)->f(x)) => last)

layerwidth = combine(alllayer,names(alllayer,Not([:siteid,:e2cfun,:WM])) .=> ByRow(i->reduce(-,i)) => identity)
layerboundary = combine(alllayer,names(alllayer,Not([:siteid,:e2cfun])).=>ByRow(last)=>identity)
layerwidthtemplate = combine(layerwidth,names(layerwidth) .=> round ∘ mean => identity)
layerboundarytemplate = combine(layerwidthtemplate,names(layerwidthtemplate)=>ByRow((i...)->cumsum([0,i...]))=>[names(layerwidthtemplate);"WM"])
ln = names(layerwidth)
lbt = collect(first(layerboundarytemplate))
nlbt = lbt./last(lbt)
layertemplate = Dict("Out"=>[-Inf,0],(ln[i]=>[nlbt[i],nlbt[i+1]] for i in eachindex(ln))...,"WM"=>[1,Inf])
save(joinpath(resultroot,"layertemplate.jld2"),"layertemplate",layertemplate)


function gettcfun(x,y,isnorm=false)
    i->begin
        ii = Spline1D([first(x)-100;x;last(x)+100],[first(y)-100;y;last(y)+100],k=1,bc="extrapolate")(i)
        isnorm ? ii./last(y) : ii
    end
end
alllayer.tcfun = [gettcfun(collect(r),lbt,true) for r in eachrow(layerboundary)]


plotlayer = (indir)-> begin

@df layerwidth groupedbar(cols(ncol(layerwidth):-1:1),size=(850,650),bar_position=:stack,yflip=true,palette=:tab10,lw=0,bar_width=1,
    xticks=1:nrow(layerwidth),xlim=(0,nrow(layerwidth)+1),xtickfontsize=6,xformatter=i->alllayer.siteid[Int(i)],xrotation=30,tickor=:out,
    leg=(0.23,0.17),legendfontsize=6,ylim=(0,2250),ylabel="Depth (μm)",xlabel="Penetration",yticks=0:100:2500)
foreach(ext->savefig(joinpath(indir,"Layers$ext")),figfmt)

p=plot(size=(250,650),palette=:tab10,yflip=true,xlim=(-0.005,0.027),ylim=(0,2250),yticks=0:100:2500,
    xticks=0:0.02:0.5,ylabel="Depth (μm)",xlabel="Probability",grid=false,tickor=:out,leftmargin=5mm)
for l in reverse(ln)
    ld = kde(layerwidth[!,l])
    plot!(p,ld.density,ld.x .+ layerboundarytemplate[1,l],label=l,lw=2)
end
ann = [(-0.0043,mean(lbt[i:i+1]),text(ln[i],7,:gray10,:left,:vcenter)) for i in 1:length(ln)]
hline!(p,lbt;linecolor=:gray25,legend=false,lw=1,ann)
foreach(ext->savefig(joinpath(indir,"LayerTemplate$ext")),figfmt)

end

plotlayer(resultroot)


deleteat!(layerwidth,15)

## Aligned Depth Feature
function collectunitdepthfeature!(indir;udf=DataFrame(),datafile="unitdepthfeature.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            d = load(joinpath(root,datafile),"df")
            df = DataFrame("siteid"=>d["siteid"],"depth"=>[d["depth"]],(d["feature"][i]=>[d["depthfeature"][:,i]] for i in 1:length(d["feature"]))...)
            append!(udf,df,cols=:union)
        end
    end
    return udf
end

function collectflashfeature!(indir;aplf=DataFrame(),test="Flash2Color")
    apfile = "$(test)_ap.jld2"
    lffile = "$(test)_lf.jld2"
    for (root,dirs,files) in walkdir(indir)
        if apfile in files
            ap = load(joinpath(root,apfile))
            lf = load(joinpath(root,lffile))
            df = DataFrame("siteid"=>ap["siteid"],"colors"=>[ap["colors"]],"titles"=>[ap["titles"]],
                "aptimes"=>[ap["times"]],"apdepths"=>[ap["depths"]],"aptrials"=>[ap["trials"]],
                "crms"=>[ap["crms"]],"pc"=>[ap["pc"]],"pdc"=>[ap["pdc"]],
                "fccsd"=>[lf["fccsd"]],"lftimes"=>[lf["times"]],"lfdepths"=>[lf["depths"]])
            append!(aplf,df,cols=:union)
        end
    end
    return aplf
end

allflashf = collectflashfeature!(resultroot)
leftjoin!(allflashf,alllayer[:,[:siteid,:e2cfun,:tcfun]],on=:siteid)
transform!(allflashf,[:tcfun,:e2cfun,:apdepths]=>ByRow((g,f,x)->g.(f.(x)))=>:apdepths,
            [:tcfun,:e2cfun,:lfdepths]=>ByRow((g,f,x)->g.(f.(x)))=>:lfdepths)
transform!(allflashf,[:fccsd,:lfdepths]=>ByRow((c,d)->camn.(c,d))=>:nfccsd)
camn(c,d) = c/maximum(abs.(c[0 .<= d .< 1,:]))



alludf = collectunitdepthfeature!(resultroot)
udf = select(alludf,[:siteid,:depth,:Density,:uppvinv,:downpvinv],[:upspread,:downspread]=>ByRow(.-)=>:spread,
    :uppvinv=>ByRow(i->inv.(i))=>:uppv,:downpvinv=>ByRow(i->inv.(i))=>:downpv)
leftjoin!(udf,alllayer[:,[:siteid,:e2cfun,:tcfun]],on=:siteid)
transform!(udf,[:tcfun,:e2cfun,:depth]=>ByRow((g,f,x)->g(f(x)))=>:depth)
transform!(udf,[:depth,:Density]=>ByRow((y,d)->mean(d[0 .<= y .< 1]))=>:cmd)

function getdffun(x,y)
    order = sortperm(x)
    i->Spline1D(x[order],y[order],k=3,bc="zero")(i)
end

plotfeature = (udf,f;dir=nothing,color=:red,xlabel="",xlim=())->begin

y = range(0,1.2,length=2000)
af = hcat(map((i,j)->getdffun(i,j)(y),udf.depth,udf[!,f])...)
xmin,xmax = isempty(xlim) ? extrema(af) : xlim
afm = mean(af,dims=2)
afse = std(af,dims=2)/sqrt(size(af,2))

p=plot(af,y;grid=false,yflip=true,color=:gray70,lw=0.5,size=(350,500),ylim=(0,1.2),xlim,ylabel="Normalized Depth",tickor=:out)
plot!(p,[afm.-afse;reverse(afm.+afse)],[y;reverse(y)]; st=:shape,lw=0,alpha=0.2,color,xlabel)
plot!(p,afm,y;label=false,lw=2,color)
ann = [(xmin+0.02(xmax-xmin),mean(nlbt[i:i+1]),text(ln[i],7,:gray10,:left,:vcenter)) for i in 1:length(ln)]
hline!(p,nlbt;linecolor=:gray25,legend=false,lw=1,ann)

isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"AlignedLayer_$f$ext")),figfmt)
end

plotfeature(udf,"Density",xlabel="Density (unit/mm³)",dir=nothing)
plotfeature(udf,"spread",xlabel="Spread (μm)",dir=nothing)
plotfeature(udf,"uppvinv",xlabel="1/Propagation Speed (ms/μm)",xlim=(-1e5,2e5),dir=nothing)
plotfeature(udf,"downpvinv",xlabel="1/Propagation Speed (ms/μm)",xlim=(-1e5,2e5),dir=nothing)
plotfeature(udf,"uppv",xlabel="Propagation Speed (μm/ms)",xlim=(-1e5,2e5),dir=nothing)
plotfeature(udf,"downpv",xlabel="Propagation Speed (μm/ms)",xlim=(-1e5,2e5),dir=nothing)

cmdp = percentile.([udf.cmd],[0,25,50,75,100])
t = filter(r->cmdp[3] <= r.cmd < cmdp[4],udf)

plotfeature(t,"Density",xlabel="Density (unit/mm³)",dir=nothing)
plotfeature(t,"spread",xlabel="Spread (μm)",dir=nothing)


pyplot()

plotcondresponse(allflashf.fccsd[1],allflashf.lftimes[1],allflashf.lfdepths[1];colors=allflashf.colors[1],titles=allflashf.titles[1],layer=layertemplate)


allflashf.lfdepths

t = allflashf.nfccsd[1][4]
heatmap(t,yflip=true)
heatmap(tt,yflip=true)
heatmap(mnfccsd,yflip=true)
heatmap(allflashf.lftimes[1],-allflashf.lfdepths[1],allflashf.nfccsd[1][4],yflip=true)



function getdffun2(x,y,z)
    (i,j)->evalgrid(Spline2D(x,y,z,kx=3,ky=3),i,j)
end

tt=getdffun2(reverse(allflashf.lfdepths[1][4]),collect(allflashf.lftimes[1][4]),reverse(allflashf.nfccsd[1][4],dims=1))(y,allflashf.lftimes[1][4])






mnfccsd=mapreduce((d,t,c)->getdffun2(reverse(d[4]),collect(t[4]),reverse(c[4],dims=1))(y,t[4]),.+,allflashf.lfdepths,allflashf.lftimes,allflashf.nfccsd)/nrow(allflashf)





t=ColorMaps["lms_mccsiso"].colors


a=HSVA(t[1])
a=HSVA(a.h,a.s,1,1)
b=HSVA(t[end])
b=HSVA(b.h,b.s,1,1)
range(RGBA(a),t[180],RGBA(b),length=360)




plotfeature2 = (allflashf,f;dir=nothing,color=:red,xlabel="",xlim=())->begin

y = range(0,1,length=2000)
af = hcat(map((i,j)->getdffun(i,j)(y),udf.depth,udf[!,f])...)
xmin,xmax = isempty(xlim) ? extrema(af) : xlim
afm = mean(af,dims=2)
afse = std(af,dims=2)/sqrt(size(af,2))

p=plot(af,y;grid=false,yflip=true,color=:gray70,lw=0.5,size=(350,500),ylim=(0,1.2),xlim,ylabel="Normalized Depth",tickor=:out)
plot!(p,[afm.-afse;reverse(afm.+afse)],[y;reverse(y)]; st=:shape,lw=0,alpha=0.2,color,xlabel)
plot!(p,afm,y;label=false,lw=2,color)
ann = [(xmin+0.02(xmax-xmin),mean(nlbt[i:i+1]),text(ln[i],7,:gray10,:left,:vcenter)) for i in 1:length(ln)]
hline!(p,nlbt;linecolor=:gray25,legend=false,lw=1,ann)

isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"AlignedLayer_$f$ext")),figfmt)
end
