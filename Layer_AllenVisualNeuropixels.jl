using NeuroAnalysis,Statistics,StatsBase,FileIO,JLD2,StatsPlots,Images,ProgressMeter,DataFrames,XLSX,Dierckx

resultroot = "Z:/"

cohcm = cgrad([HSL(170,1,0.94),HSL(190,1,0.56),HSL(233,1,0.38)])
powcm = cgrad([HSL(60,1,0.92),HSL(40,1,0.62),HSL(0,1,0.32)])
csdcm = cgrad([HSL(358,0.7,0.3),HSL(358,0.8,0.5),HSL(58,0.9,0.6),HSL(125,0.4,0.7),HSL(192,0.9,0.6),HSL(238,0.8,0.5),HSL(238,0.7,0.3)],
              [0,0.12,0.37,0.5,0.63,0.88,1])
absmax(x) = mapreduce(i->maximum(abs.(i)),max,x)
absmin(x) = mapreduce(i->minimum(abs.(i)),min,x)
abspct(x,p=99.9) = percentile(vec(abs.(x)),p)
abspct2(x,p=99.9) = percentile.((vec(abs.(x)),),(100-p,p))

plotunitfeatureallen=(resps,depths;ylims=(0,3820),xl="",titles="",titlefontsize=7,w=140,h=600,color=:red,layers=[],nss=nothing)->begin
    n = length(resps)
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:200:ylims[2]) : false
        leftmargin = i==1 ? 8Plots.mm : -4Plots.mm
        ylabel = i==1 ? "Depth (μm)" : ""
        xlabel = i==1 ? xl : ""
        title = isempty(titles) ? titles : titles[i]

        xmin,xmax = extrema(resps[i])
        annxl = xmin+0.02(xmax-xmin)
        annxr = xmin+0.98(xmax-xmin)

        if !isnothing(nss)
            ns = nss[i]
            ann = [(annxr,mean(ns[k]),text(k,6,:gray10,:right,:vcenter)) for k in keys(ns)]
            hspan!(p[i],mapreduce(l->[l...],append!,values(ns));linecolor=:gray25,color=coloralpha(colorant"gray65",0.1),legend=false,lw=0.5,ann)
        end

        plot!(p[i],resps[i],depths[i];color,ylims,
        title,titlefontsize,yticks,tickor=:out,xlabel,ylabel,lw=2,
        xticks = range(start=ceil(xmin,sigdigits=1),stop=floor(0.75xmax,sigdigits=1),length=2),
        leftmargin,bottommargin=6Plots.mm,topmargin=12Plots.mm)

        if !isempty(layers)
            layer = layers[i]
            if !isempty(layer)
                ann = [(annxl,mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
                hline!(p[i],[l[end] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ls=:dash,ann)
            end
        end
    end
    p
end

plotcondresponseallen=(resps,times,depths;colors=[],minmaxcolors=[],xl="Time (ms)",xt=[1000,2000],rl="",rlims=(0,1),rscale=1,
                        rcolor=HSL(120,0.5,0.25),titles="",titlefontsize=7,w=140,h=600,color=:coolwarm,layers=[],nss=nothing)->begin
    n = length(resps)
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:200:depths[i][end]) : false
        xticks = xl == "Time (ms)" ? (0:50:times[i][end]) : xt
        leftmargin = i==1 ? 8Plots.mm : -4Plots.mm
        xlabel = i==1 ? xl : ""
        rlabel = i==1 ? rl : ""
        ylabel = i==1 ? "Depth (μm)" : ""
        title = isempty(titles) ? titles : titles[i]

        if isnothing(nss)
            lim = absmax(resps)
            clims = (-lim,lim)
        else
            ns = nss[i]
            visk = filter(i->startswith(String(i),"VIS"),keys(ns))|>first
            if xl == "Time (ms)"
                lim=abspct(resps[i][ns[visk][1].<=depths[i].<=ns[visk][2],:])
                clims = (-lim,lim)
            else
                clims=abspct2(resps[i][ns[visk][1].<=depths[i].<=ns[visk][2],:])
            end
        end

        xmin,xmax = extrema(times[i])
        ymin,ymax = extrema(depths[i])
        annxl = xmin+0.02(xmax-xmin)
        annxr = xmin+0.98(xmax-xmin)
        anny = ymin+0.01(ymax-ymin)
        adx = 30
        ann = isempty(colors) ? [] : [(annxl,anny,Plots.text("■",15,colors[i],:left,:bottom))]
        ann = isempty(minmaxcolors) ? ann : [(annxl,anny,Plots.text("▮",15,minmaxcolors[i][2],:left,:bottom)),(annxl+adx,anny,Plots.text("▮",15,minmaxcolors[i][1],:left,:bottom))]

        heatmap!(p[i],times[i],depths[i],resps[i];color,clims,ylims=extrema(depths[i]),
        title,titlefontsize,yticks,xticks,tickor=:out,xlabel,ylabel,
        ann,leftmargin,bottommargin=6Plots.mm,topmargin=12Plots.mm)

        if xl != "Time (ms)"
            rm = dropdims(mean(resps[i]*rscale,dims=2),dims=2)
            rmse = dropdims(std(resps[i]*rscale,dims=2),dims=2)/sqrt(size(resps[i],2))
            rticks=[rlims[1],round(0.5(rlims[2]-rlims[1]),sigdigits=1)]
            pp=mytwiny(p[i])
            plot!(pp,[rm.-rmse;reverse(rm.+rmse)],[depths[i];reverse(depths[i])]; st=:shape,lw=0,alpha=0.25,color=rcolor,leftmargin,xlabel=rlabel)
            plot!(pp,rm,depths[i];color=rcolor,leg=false,tickor=:out,xlims=rlims,xticks=rticks,ylims=extrema(depths[i]),lw=1,leftmargin)
        end

        if !isnothing(nss)
            ns=nss[i]
            ann = [(annxr,mean(ns[k]),text(k,6,:gray10,:right,:vcenter)) for k in keys(ns)]
            hline!(p[i],mapreduce(l->[l...],append!,values(ns));linecolor=:gray25,legend=false,lw=0.5,ann)
        end

        if !isempty(layers)
            layer = layers[i]
            if !isempty(layer)
                ann = [(annxl,mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
                hline!(p[i],[l[end] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ls=:dash,ann)
            end
        end
    end
    p
end

plottrialresponseallen=(resps,trials,depths;colors=[],minmaxcolors=[],xl="Trial",rl="",rlims=(0,1),rscale=1,
                        rcolor=HSL(120,0.5,0.25),titles="",titlefontsize=7,w=140,h=600,color=:heat,layers=[],nss=nothing)->begin
    n = length(resps)
    clims = (absmin(resps),absmax(resps))
    p=plot(layout=(1,n),link=:y,legend=false,grid=false,size=(n*w,h))
    for i in 1:n
        yticks = i==1 ? (0:200:depths[i][end]) : false
        xticks = range(start=1,stop=round(Int,0.75*trials[i][end]),length=2)
        leftmargin = i==1 ? 6Plots.mm : -4Plots.mm
        xlabel = i==1 ? xl : ""
        rlabel = i==1 ? rl : ""
        ylabel = i==1 ? "Depth (μm)" : ""
        title = isempty(titles) ? titles : titles[i]

        if !isnothing(nss)
            ns = nss[i]
            visk = filter(i->startswith(String(i),"VIS"),keys(ns))|>first
            clims=abspct2(resps[i][ns[visk][1].<=depths[i].<=ns[visk][2],:])
        end

        xmin,xmax = extrema(trials[i])
        ymin,ymax = extrema(depths[i])
        annxl = xmin+0.02(xmax-xmin)
        annxr = xmin+0.98(xmax-xmin)
        anny = ymin+0.01(ymax-ymin)
        adx = 30
        ann = isempty(colors) ? [] : [(annxl,anny,Plots.text("▮",15,colors[2+2(i-1)],:left,:bottom)),(annxl+adx,anny,Plots.text("▮",15,colors[1+2(i-1)],:left,:bottom))]
        ann = isempty(minmaxcolors) ? ann : [(annxl,anny,Plots.text("▮",15,minmaxcolors[i][2],:left,:bottom)),(annxl+adx,anny,Plots.text("▮",15,minmaxcolors[i][1],:left,:bottom))]

        heatmap!(p[i],trials[i],depths[i],resps[i];color,clims,ylims=extrema(depths[i]),
        title,titlefontsize,yticks,xticks,tickor=:out,xlabel,ylabel,
        ann,leftmargin,bottommargin=3Plots.mm)

        rm = dropdims(mean(resps[i]*rscale,dims=2),dims=2)
        rmse = dropdims(std(resps[i]*rscale,dims=2),dims=2)/sqrt(size(resps[i],2))
        rticks=[rlims[1],round(0.5(rlims[2]-rlims[1]),sigdigits=1)]
        pp=mytwiny(p[i])
        plot!(pp,[rm.-rmse;reverse(rm.+rmse)],[depths[i];reverse(depths[i])]; st=:shape,lw=0,alpha=0.25,color=rcolor,leftmargin,xlabel=rlabel)
        plot!(pp,rm,depths[i];color=rcolor,leg=false,tickor=:out,xlims=rlims,xticks=rticks,ylims=extrema(depths[i]),lw=1,leftmargin)

        if !isnothing(nss)
            ns = nss[i]
            ann = [(annxr,mean(ns[k]),text(k,6,:gray10,:right,:vcenter)) for k in keys(ns)]
            hline!(p[i],mapreduce(l->[l...],append!,values(ns));linecolor=:gray25,legend=false,lw=0.5,ann)
        end

        if !isempty(layers)
            layer = layers[i]
            if !isempty(layer)
                ann = [(annxl,mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
                hline!(p[i],[l[end] for l in values(layer)];linecolor=:gray25,legend=false,lw=0.5,ls=:dash,ann)
            end
        end
    end
    p
end

function mytwiny(sp, letter=:y)
    plt = sp.plt
    # orig_sp = first(plt.subplots)
    orig_sp = sp
    for letter in filter(!=(letter), Plots.axes_letters(orig_sp, letter))
        ax = orig_sp[Plots.get_attr_symbol(letter, :axis)]
        ax[:grid] = false  # disable the grid (overlaps with twin axis)
    end
    if orig_sp[:framestyle] === :box
        # incompatible with shared axes (see github.com/JuliaPlots/Plots.jl/issues/2894)
        orig_sp[:framestyle] = :axes
    end
    plot!(
        # plt;
        sp;
        inset = (sp[:subplot_index], bbox(0, 0, 1, 1)),
        left_margin = orig_sp[:left_margin],
        top_margin = orig_sp[:top_margin],
        right_margin = orig_sp[:right_margin],
        bottom_margin = orig_sp[:bottom_margin],
    )
    twin_sp = last(plt.subplots)
    letters = Plots.axes_letters(twin_sp, letter)
    tax, oax = map(l -> twin_sp[Plots.get_attr_symbol(l, :axis)], letters)
    tax[:grid] = false
    tax[:showaxis] = false
    tax[:ticks] = :none
    oax[:grid] = false
    oax[:mirror] = true
    twin_sp[:background_color_inside] = RGBA{Float64}(0, 0, 0, 0)
    Plots.link_axes!(sp[Plots.get_attr_symbol(letter, :axis)], tax)
    twin_sp
end



## Manual Layer Assignment
session_id = 847657808
probe_id = 848037572
sessionresultdir = joinpath(resultroot,"$session_id")
sessionproberesultdir = joinpath(sessionresultdir,"$probe_id")

layer = load(joinpath(sessionproberesultdir,"layer.jld2"),"layer")
# layer = Dict()

layer["1"] = [1100,2285]
layer["2/3"] = [1000,2200]
layer["4"] = [1000,1970]
layer["5"] = [1800,1780]
layer["6"] = [1600,1595]
layer["WM"] = [1300,1390]

# Finalize Layer
layer = checklayer!(layer)
jldsave(joinpath(sessionproberesultdir,"layer.jld2");layer,session_id,probe_id)

# Plot layer and responses
allenflash(session_id,sessionresultdir;spudf,layer=:batch)



function allenflash(session_id,sessionresultdir;spudf=nothing,layer=nothing,figfmt=[".png"])

    aps=[];lfs=[];pls=[]
    for (root,dirs,files) in walkdir(sessionresultdir)
        if "ap.jld2" in files
            push!(aps,load(joinpath(root,"ap.jld2")))
            push!(lfs,load(joinpath(root,"lf.jld2")))
            if layer == :batch
                layerfile = joinpath(root,"layer.jld2")
                l = isfile(layerfile) ? load(layerfile,"layer") : Dict()
                push!(pls,l)
            end
        end
    end


    ## Unit
    if !isnothing(spudf)
        pidx = [findfirst(r->r.probe_id==i["probe_id"] && r.session_id==session_id,eachrow(spudf)) for i in aps]
        depths = spudf.depth[pidx]
        titles = spudf.probe_id[pidx]
        nss = spudf.neurostructure[pidx]

        fs = ["Density","Spread","Duration","PTR"]
        fc = palette(:tab10,length(fs)).colors.colors

        for i in eachindex(fs)
            plotunitfeatureallen(spudf[pidx,fs[i]],depths;color=fc[i],xl=fs[i],titles,layers=pls,nss)
            foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_Unit$(fs[i])$ext")),figfmt)
        end
    end


    ## AP
    titles = repeat(["$(i["probe_id"])" for i in aps],inner=2)
    plss = repeat(pls,inner=2)
    colors = mapreduce(i->[contains(k,"-1") ? Gray(0) : Gray(1) for k in keys(i["crms"])],append!,aps)
    minmaxcolors=[colors[i:i+1] for i in 1:2:length(colors)]
    nss = mapreduce(i->[i["neurostructure"],i["neurostructure"]],append!,aps)
    depths = mapreduce(i->[i["depths"],i["depths"]],append!,aps)

    crms = mapreduce(i->[v for v in values(i["crms"])],append!,aps)
    times = mapreduce(i->[i["times"],i["times"]],append!,aps)
    plotcondresponseallen(crms,times,depths;colors,titles,layers=plss,nss)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_dRMS$ext")),figfmt)

    tps = map(i->i["tps"],aps)
    trials = map(i->1:size(i,2),tps)
    plottrialresponseallen(tps,trials,depths;minmaxcolors,titles=titles[1:2:end],layers=pls,nss=nss[1:2:end],color=powcm,rl="Power (μV²)",rlims=(0,2e5))
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_tPower$ext")),figfmt)
    # σ=1(20μm): 5(80μm) diameter gaussian kernal to filter depth line power 
    ftps = map(i -> imfilter(i,Kernel.gaussian((1,0)),Fill(0)), tps)
    plottrialresponseallen(ftps, trials, depths; minmaxcolors, titles=titles[1:2:end], layers=pls, nss=nss[1:2:end], color=powcm,rl="Power (μV²)",rlims=(0,2e5))
    foreach(ext -> savefig(joinpath(sessionresultdir, "$(isnothing(layer) ? "" : "Layer_")$(session_id)_tPowerf1$ext")), figfmt)
    tlc = map(i -> i["tlc"], aps)
    plottrialresponseallen(tlc, trials, depths; minmaxcolors, titles=titles[1:2:end], layers=pls, nss=nss[1:2:end], color=cohcm,rl="Coherence")
    foreach(ext -> savefig(joinpath(sessionresultdir, "$(isnothing(layer) ? "" : "Layer_")$(session_id)_tCoherence$ext")), figfmt)

    cfps = mapreduce(i -> [v for v in values(i["cfps"])], append!, aps)
    psfreqs = mapreduce(i->[i["psfreqs"],i["psfreqs"]],append!,aps)
    plotcondresponseallen(cfps, psfreqs, depths;xl="Frequency (Hz)", colors, titles, layers=plss, nss,color=powcm,rl="Power (μV²)",rlims=(0,2e5))
    foreach(ext -> savefig(joinpath(sessionresultdir, "$(isnothing(layer) ? "" : "Layer_")$(session_id)_fPower$ext")), figfmt)
    # σ=1(20μm): 5(80μm) diameter gaussian kernal to filter depth line power
    fcfps = map(i -> imfilter(i,Kernel.gaussian((1,0)),Fill(0)), cfps)
    # w=21(147Hz) frequency window to filter line noises
    fcfps = map(i -> mapwindow(minimum,i,(1,21)), fcfps)
    plotcondresponseallen(fcfps, psfreqs, depths;xl="Frequency (Hz)", colors, titles, layers=plss, nss,color=powcm,rl="Power (μV²)",rlims=(0,2e5))
    foreach(ext -> savefig(joinpath(sessionresultdir, "$(isnothing(layer) ? "" : "Layer_")$(session_id)_fPowerf1w21$ext")), figfmt)
    
    cflc = mapreduce(i -> [v for v in values(i["cflc"])], append!, aps)
    lcfreqs = mapreduce(i->[i["lcfreqs"],i["lcfreqs"]],append!,aps)
    plotcondresponseallen(cflc, lcfreqs, depths;xl="Frequency (Hz)", colors, titles, layers=plss, nss,color=cohcm,rl="Coherence")
    foreach(ext -> savefig(joinpath(sessionresultdir, "$(isnothing(layer) ? "" : "Layer_")$(session_id)_fCoherence$ext")), figfmt)
    # w=21(147Hz) frequency window to filter line noises
    fcflc = map(i -> mapwindow(minimum,i,(1,21)), cflc)
    plotcondresponseallen(fcflc, lcfreqs, depths;xl="Frequency (Hz)", colors, titles, layers=plss, nss,color=cohcm,rl="Coherence")
    foreach(ext -> savefig(joinpath(sessionresultdir, "$(isnothing(layer) ? "" : "Layer_")$(session_id)_fCoherencew21$ext")), figfmt)
    

    # condition combined AP
    cctimes = times[1:2:end]
    ccdepths = depths[1:2:end]
    cctitles = titles[1:2:end]
    ccnss = nss[1:2:end]

    ccrms = mapreduce(i->[reduce((i,j)->(i.+j)/2,values(i["crms"]))],append!,aps)
    plotcondresponseallen(ccrms,cctimes,ccdepths;minmaxcolors,titles=cctitles,nss=ccnss,layers=pls)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_ccdRMS$ext")),figfmt)

    ccfps = mapreduce(i->[reduce((i,j)->(i.+j)/2,values(i["cfps"]))],append!,aps)
    ccpsfreqs = psfreqs[1:2:end]
    # σ=1(20μm): 5(80μm) diameter gaussian kernal to filter depth line power
    fccfps = map(i -> imfilter(i,Kernel.gaussian((1,0)),Fill(0)), ccfps)
    # w=21(147Hz) frequency window to filter line noises
    fccfps = map(i -> mapwindow(minimum,i,(1,21)), fccfps)
    plotcondresponseallen(fccfps, ccpsfreqs, ccdepths;xl="Frequency (Hz)", minmaxcolors, titles=cctitles, layers=pls, nss=ccnss,color=powcm,rl="Power (μV²)",rlims=(0,2e5))
    foreach(ext -> savefig(joinpath(sessionresultdir, "$(isnothing(layer) ? "" : "Layer_")$(session_id)_ccfPowerf1w21$ext")), figfmt)
    
    ccflc = mapreduce(i->[reduce((i,j)->(i.+j)/2,values(i["cflc"]))],append!,aps)
    cclcfreqs = lcfreqs[1:2:end]
    # w=21(147Hz) frequency window to filter line noises
    fccflc = map(i -> mapwindow(minimum,i,(1,21)), ccflc)
    plotcondresponseallen(fccflc,cclcfreqs,ccdepths;xl="Frequency (Hz)",minmaxcolors,titles=cctitles,nss=ccnss,layers=pls,color=cohcm,rl="Coherence")
    foreach(ext -> savefig(joinpath(sessionresultdir, "$(isnothing(layer) ? "" : "Layer_")$(session_id)_ccfCoherencew21$ext")), figfmt)
    
    
    ## LFP and CSD
    clfp = mapreduce(i -> [v for v in values(i["clfp"])], append!, lfs)
    cdcsd = mapreduce(i -> [v for v in values(i["ccsd"])], append!, lfs)
    times = mapreduce(i -> [i["times"], i["times"]], append!, lfs)
    depths = mapreduce(i -> [i["depths"], i["depths"]], append!, lfs)
    hy = depths[1][2] - depths[1][1]
    baseindex = lfs[1]["baseindex"]

    plotcondresponseallen(clfp,times,depths;colors,titles,layers=plss,nss)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_LFP$ext")),figfmt)
    plotcondresponseallen(cdcsd,times,depths;colors,titles,layers=plss,nss,color=csdcm)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_dCSD$ext")),figfmt)

    # σ=1.5(30μm): 9(160μm) diameter gaussian kernal
    fcdcsd = map(i -> imfilter(i, Kernel.gaussian((1.5, 0)), Fill(0)), cdcsd)
    plotcondresponseallen(fcdcsd,times,depths;colors,titles,layers=plss,nss,color=csdcm)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_dCSDf1.5$ext")),figfmt)

    # CSD
    ccsd = map(v->csd(v,h=hy),clfp)
    fccsd = map(i->imfilter(i,Kernel.gaussian((1.5,0)),Fill(0)),ccsd)
    plotcondresponseallen(fccsd,times,depths;colors,titles,layers=plss,nss,color=csdcm)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_CSDf1.5$ext")),figfmt)


    # condition combined LFP and CSD
    cctimes = times[1:2:end]
    ccdepths = depths[1:2:end]
    cctitles = titles[1:2:end]
    ccnss = nss[1:2:end]
    cclfp = mapreduce(i->[reduce((i,j)->(i.+j)/2,values(i["clfp"]))],append!,lfs)
    plotcondresponseallen(cclfp,cctimes,ccdepths;minmaxcolors,titles=cctitles,nss=ccnss,layers=pls)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_ccLFP$ext")),figfmt)

    # ccCSD
    cccsd = map(v->csd(v,h=hy),cclfp)
    fcccsd = map(i->imfilter(i,Kernel.gaussian((1.5,0)),Fill(0)),cccsd)
    plotcondresponseallen(fcccsd,cctimes,ccdepths;minmaxcolors,titles=cctitles,nss=ccnss,layers=pls,color=csdcm)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_ccCSDf1.5$ext")),figfmt)

    # ccdCSD
    fccdcsd = map(v->imfilter(stfilter(v,temporaltype=:sub,ti=baseindex),Kernel.gaussian((1.5,0)),Fill(0)),cccsd)
    plotcondresponseallen(fccdcsd,cctimes,ccdepths;minmaxcolors,titles=cctitles,nss=ccnss,layers=pls,color=csdcm)
    foreach(ext->savefig(joinpath(sessionresultdir,"$(isnothing(layer) ? "" : "Layer_")$(session_id)_ccdCSDf1.5$ext")),figfmt)
end


## batch each session
@showprogress "Batch Allen Session Plots ... " for s in levels(sps.session_id)
    # allenflash(s,joinpath(resultroot,"$s"))
    allenflash(s,joinpath(resultroot,"$s");spudf,layer=:batch,figfmt=[".png",".svg"])
end