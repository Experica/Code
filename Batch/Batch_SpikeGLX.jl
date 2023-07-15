using StatsBase,StatsPlots,Images,DSP,ePPR

function process_flash_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files),spikesorter=param[:spikesorter])

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    datadir = joinpath(param[:dataroot],subject,siteid)
    siteresultdir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(siteresultdir,testid)
    isdir(resultdir) || mkpath(resultdir)

    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"]
    unitgood=spike["unitgood"];unitposition=spike["unitposition"];unitsync=spike["unitsync"]
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    jldsave(joinpath(resultdir,"spike.jld2");spike,siteid)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff.-condon)
    exenv=Dict()
    exenv["ID"] = ex["ID"]
    exenv["conddur"] = conddur
    if haskey(ex,"Eye")
        exenv["eye"] = ex["Eye"]
    else
        exenv["eye"] = ex["Log"]
    end
    exenv["color"] = "$(exparam["ColorSpace"])_$(exparam["Color"])"

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)

    # for ii in dataset["imecindex"]
    #     # Prepare AP
    #     # d = "ap$ii"
    #     # meta = dataset[d]["meta"]
    #     # binfile = meta["fileName"]
    #     # nsavedch = meta["nSavedChans"]
    #     # nsample = meta["nFileSamp"]
    #     # fs = meta["fs"]
    #     # nch = nsavedch-1 # exclude sync channel
    #     # hx,hy,hz = meta["probespacing"]
    #     # exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
    #     # pnrow,pncol = size(exchmask)
    #     # depths = hy*(0:(pnrow-1))
    #     # mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)

    #     # Prepare AP
    #     d = "ap$ii"
    #     meta = dataset[d]["meta"]
    #     binfile = spike["datapath"]
    #     chmapraw = spike["chmapraw"]
    #     nch = spike["nch"]
    #     nsample = spike["nsample"]
    #     hx,hy,hz = meta["probespacing"]
    #     pversion = meta["probeversion"]
    #     t0 = spike["t0"]
    #     winvraw = spike["winvraw"]
    #     exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
    #     pnrow,pncol = size(exchmask)
    #     depths = hy*(0:(pnrow-1))
    #     mmbinfile = mmap(binfile,Matrix{Int16},(nch,nsample),0)
    #     synccondon = ref2sync(condon,dataset,ii)
    #     exenv["hy"] = hy
    #     exenv["t0"] = t0
    #     exenv["synccondon"] = synccondon
    #     exenv["exchmask"] = exchmask
    #     exenv["chmapraw"] = chmapraw
    #     exenv["nch"] = nch
    #     exenv["nsample"] = nsample
    #     exenv["winvraw"] = winvraw
    #     exenv["fs"] = spike["fs"]
    #     exenv["binfile"] = binfile

    #     # # epoch AP
    #     # epoch = [-40 150]
    #     # epochs = synccondon.+t0.+epoch
    #     # # All AP epochs(mapped to concat file), unwhiten, gain corrected(voltage), bandpass filtered,
    #     # # with all channels mapped in the shape of probe where excluded channels are replaced with rand local average
    #     # ys = fill2mask(epochsamplenp(mmbinfile,spike["fs"],epochs,1:nch;meta,bandpass=[300,3000],whiten=winvraw),exchmask,chmap=chmapraw,randreplace=true)
    #     # # 3 fold downsampling(10kHz)
    #     # ys = resample(ys,1//3,dims=3)
    #     # fs = spike["fs"]/3
    #     # baseindex = epoch2sampleindex([0 50],fs) # [-40, 10] ms

    #     # # power spectrum of same depth
    #     # @views pys = ys[.!exchmask,:,:] # exclude and flat channels
    #     # chpos = vcat(chpositionnp(pversion)[.!exchmask]...) # exclude and flat channel positions 
    #     # chgi = [findall(chpos[:,2].==up) for up in unique(chpos[:,2])] # group channels with same depth
    #     # @views ppss,psfreqs = powerspectrum(pys[:,baseindex[end]+1:end,:],fs;freqrange=[300,3000])
    #     # pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
    #     # @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    #     # tps = dropdims(mean(pss,dims=2),dims=2) # freq average
    #     # @views cfps = Dict(condstring(r)=>
    #     #     dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average       
    #     # for r in eachrow(cond))

    #     # # power contrast of same depth
    #     # prmss = mapwindow(x->rms(x)^2,pys,(1,101,1),border="symmetric") # 10ms rms window
    #     # rmss = Array{Float64}(undef,length(chgi),size(prmss)[2:end]...)
    #     # @views foreach(i->rmss[i,:,:] = mean(prmss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    #     # @views crms = Dict(condstring(r)=>
    #     #     stfilter(dropdims(mean(rmss[:,:,r.i],dims=3),dims=3),temporaltype=:rcb,ti=baseindex) # trial average
    #     # for r in eachrow(cond))
    #     # times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(pys,2))

    #     # # local coherence
    #     # @views lcs,lcfreqs = localcoherence(pys[:,baseindex[end]+1:end,:],chpos,fs;freqrange=[300,3000],lr=55,sigma=25,chgroupdim=2)
    #     # tlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
    #     # @views cflc = Dict(condstring(r)=>
    #     #     dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average       
    #     # for r in eachrow(cond))

    #     # # @views pc,freqs = coherencespectrum(ys[:,baseindex[end]+1:end,:],fs,freqrange=[300,3000],ismean=true) # freq average

    #     # # dsn = 2:5 # downsample n, 40,60,80,100μm for hy=20μm
    #     # # dsr = [3,2,1,1] # downsample r, 120,120,80,100μm for dsn
    #     # # @views pcd = Dict(dsn[d]=>hcat(map(i->bandmean(pc[1:dsn[d]:end,1:dsn[d]:end,i],r=dsr[d],s=1),1:size(pc,3))...) for d in 1:4)
    #     # # @views pc = hcat(map(i->bandmean(pc[:,:,i],r=5,s=1),1:size(pc,3))...) # local average with gaussian weights of -100:100μm, σ=20μm
        
    #     # if plot
    #     #     plotanalog(tps;hy,color=:heat,n=mean(tps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tPower$ext")),figfmt)
    #     #     plotanalog(tlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(tlc,dims=2),layer)
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tCoherence$ext")),figfmt)
    #     #     fpslims = extrema(mapreduce(extrema,union,values(cfps)))
    #     #     rmslim = mapreduce(x->maximum(abs.(x)),max,values(crms))
    #     #     flclims = extrema(mapreduce(extrema,union,values(cflc)))
    #     #     for k in keys(crms)
    #     #         plotanalog(cfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(cfps[k],dims=2),layer)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fPower$ext")),figfmt)
    #     #         plotanalog(crms[k];x=times,hy,clims=(-rmslim,rmslim),layer)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dRMS$ext")),figfmt)
    #     #         plotanalog(cflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(cflc[k],dims=2),layer)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fCoherence$ext")),figfmt)
    #     #     end
    #     # end
    #     # jldsave(joinpath(resultdir,"ap$(ii).jld2");tps,cfps,crms,tlc,cflc,times,psfreqs,lcfreqs,depths,baseindex,fs,siteid,exenv)


    #     # # baseline power spectrum of same depth
    #     # @views ppss,psfreqs = powerspectrum(pys[:,baseindex,:],fs;freqrange=[300,3000])
    #     # pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
    #     # @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    #     # tps = dropdims(mean(pss,dims=2),dims=2) # freq average
    #     # @views cfps = Dict(condstring(r)=>
    #     #     dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average       
    #     # for r in eachrow(cond))

    #     # # baseline local coherence
    #     # @views lcs,lcfreqs = localcoherence(pys[:,baseindex,:],chpos,fs;freqrange=[300,3000],lr=55,sigma=25,chgroupdim=2)
    #     # tlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
    #     # @views cflc = Dict(condstring(r)=>
    #     #     dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average       
    #     # for r in eachrow(cond))
        
    #     # if plot
    #     #     plotanalog(tps;hy,color=:heat,n=mean(tps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tbPower$ext")),figfmt)
    #     #     plotanalog(tlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(tlc,dims=2),layer)
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tbCoherence$ext")),figfmt)
    #     #     fpslims = extrema(mapreduce(extrema,union,values(cfps)))
    #     #     flclims = extrema(mapreduce(extrema,union,values(cflc)))
    #     #     for k in keys(cfps)
    #     #         plotanalog(cfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(cfps[k],dims=2),layer)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fbPower$ext")),figfmt)
    #     #         plotanalog(cflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(cflc[k],dims=2),layer)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fbCoherence$ext")),figfmt)
    #     #     end
    #     # end
    #     # jldsave(joinpath(resultdir,"ap$(ii)+.jld2");tps,cfps,tlc,cflc,psfreqs,lcfreqs,depths,baseindex,fs,siteid,exenv)

    #     # # Whiten epoch AP
    #     # epoch = [-40 150]
    #     # epochs = synccondon.+t0.+epoch
    #     # # All AP epochs(mapped to concat file), remain whiten, bandpass filtered,
    #     # # with all channels mapped in the shape of probe where excluded channels are replaced with rand local average
    #     # ys = fill2mask(epochsamplenp(mmbinfile,spike["fs"],epochs,1:nch;bandpass=[300,3000]),exchmask,chmap=chmapraw,randreplace=true)
    #     # # 3 fold downsampling(10kHz)
    #     # ys = resample(ys,1//3,dims=3)
    #     # fs = spike["fs"]/3
    #     # baseindex = epoch2sampleindex([0 50],fs) # [-40, 10] ms

    #     # # power spectrum of same depth
    #     # @views pys = ys[.!exchmask,:,:] # exclude and flat channels
    #     # chpos = vcat(chpositionnp(pversion)[.!exchmask]...) # exclude and flat channel positions 
    #     # chgi = [findall(chpos[:,2].==up) for up in unique(chpos[:,2])] # group channels with same depth
    #     # @views ppss,psfreqs = powerspectrum(pys[:,baseindex[end]+1:end,:],fs;freqrange=[300,3000])
    #     # pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
    #     # @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    #     # wtps = dropdims(mean(pss,dims=2),dims=2) # freq average
    #     # @views wcfps = Dict(condstring(r)=>
    #     #     dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average       
    #     # for r in eachrow(cond))
        
    #     # # power contrast of same depth
    #     # prmss = mapwindow(x->rms(x)^2,pys,(1,101,1),border="symmetric") # 10ms rms window
    #     # rmss = Array{Float64}(undef,length(chgi),size(prmss)[2:end]...)
    #     # @views foreach(i->rmss[i,:,:] = mean(prmss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    #     # @views wcrms = Dict(condstring(r)=>
    #     #     stfilter(dropdims(mean(rmss[:,:,r.i],dims=3),dims=3),temporaltype=:rcb,ti=baseindex) # trial average
    #     # for r in eachrow(cond))
    #     # times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(pys,2))

    #     # # local coherence
    #     # @views lcs,lcfreqs = localcoherence(pys[:,baseindex[end]+1:end,:],chpos,fs;freqrange=[300,3000],lr=55,sigma=25,chgroupdim=2)
    #     # wtlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
    #     # @views wcflc = Dict(condstring(r)=>
    #     #     dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average       
    #     # for r in eachrow(cond))

    #     # if plot
    #     #     plotanalog(wtps;hy,color=:heat,n=mean(wtps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_wtPower$ext")),figfmt)
    #     #     plotanalog(wtlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(wtlc,dims=2),layer)
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_wtCoherence$ext")),figfmt)
    #     #     fpslims = extrema(mapreduce(extrema,union,values(wcfps)))
    #     #     rmslim = mapreduce(x->maximum(abs.(x)),max,values(wcrms))
    #     #     flclims = extrema(mapreduce(extrema,union,values(wcflc)))
    #     #     for k in keys(wcrms)
    #     #         plotanalog(wcfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(wcfps[k],dims=2),layer)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_wfPower$ext")),figfmt)
    #     #         plotanalog(wcrms[k];x=times,hy,clims=(-rmslim,rmslim),layer)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_wdRMS$ext")),figfmt)
    #     #         plotanalog(wcflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(wcflc[k],dims=2),layer)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_wfCoherence$ext")),figfmt)
    #     #     end
    #     # end
    #     # jldsave(joinpath(resultdir,"wap$(ii).jld2");wtps,wcfps,wcrms,wtlc,wcflc,times,psfreqs,lcfreqs,depths,baseindex,fs,siteid,exenv)
        
        
    #     # Depth Unit PSTH
    #     epoch = [-40 150]
    #     epochs = synccondon.+epoch
    #     bw = 2 # ms
    #     psthbins = epoch[1]:bw:epoch[2]
    #     baseindex = epoch2sampleindex([0 50],1/(bw*SecondPerUnit)) # [-40, 10] ms
        
    #     ui = unitsync.==ii
    #     @views unitepochpsth = map(st->psthspiketrains(epochspiketrain(st,epochs,isminzero=true,shift=-epoch[1]).y,psthbins,israte=true,ismean=false),unitspike[ui])
    #     x = unitepochpsth[1].x
    #     @views cpsth = Dict(condstring(r)=>begin
    #                         ucp = map(p->vmeanse(p.mat[r.i,:]).m,unitepochpsth)
    #                         p = spacepsth(ucp,unitposition[ui,:],dims=2,spacerange=depths,bw=2hy,step=hy)
    #                         (;psth=mapwindow(mean,p.psth,(1,5),border="symmetric"),p.y) # 10ms mean window
    #                     end  for r in eachrow(cond))
    #     y = first(values(cpsth)).y
    #     cpsth = Dict(k=>cpsth[k].psth for k in keys(cpsth))
    #     cdpsth = Dict(k=>stfilter(cpsth[k],temporaltype=:sub,ti=baseindex) for k in keys(cpsth))

    #     if plot
    #         plim = mapreduce(i->maximum(abs.(i)),max,values(cdpsth))
    #         for k in keys(cdpsth)
    #             plotanalog(cdpsth[k];x,y,clims=(-plim,plim))
    #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dPSTH$ext")),figfmt)
    #         end
    #     end
    #     jldsave(joinpath(resultdir,"psth$(ii).jld2");cpsth,cdpsth,x,y,baseindex,siteid,exenv)

    #     continue
    #     # Prepare LF
    #     d = "lf$ii"
    #     meta = dataset[d]["meta"]
    #     binfile = meta["fileName"]
    #     nsavedch = meta["nSavedChans"]
    #     nsample = meta["nFileSamp"]
    #     nch = nsavedch-1 # exclude sync channel
    #     hx,hy,hz = meta["probespacing"]
    #     pversion = meta["probeversion"]
    #     t0 = 0 # not concat drift correct binary yet
    #     exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
    #     pnrow,pncol = size(exchmask)
    #     depths = hy*(0:(pnrow-1))
    #     mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)
    #     synccondon = ref2sync(condon,dataset,ii)
    #     exenv["hy"] = hy
    #     exenv["t0"] = t0
    #     exenv["synccondon"] = synccondon
    #     exenv["exchmask"] = exchmask
    #     exenv["nch"] = nch
    #     exenv["nsample"] = nsample
    #     exenv["fs"] = meta["fs"]
    #     exenv["binfile"] = binfile

    #     # epoch LF
    #     epoch = [-40 150]
    #     epochs = synccondon.+t0.+epoch
    #     # All LF epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered,
    #     # with all channels mapped in the shape of probe where excluded channels are replaced with local average
    #     ys = fill2mask(epochsamplenp(mmbinfile,meta["fs"],epochs,1:nch;meta,bandpass=[1,100]),exchmask)
    #     # 2.5 fold downsampling(1kHz)
    #     ys = resample(ys,1/2.5,dims=3)
    #     fs = meta["fs"]/2.5
    #     baseindex = epoch2sampleindex([0 50],fs) # [-40, 10] ms

    #     # @views ys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    #     # @views clfp = Dict(condstring(r)=>
    #     #     dropdims(mean(ys[:,:,vcat((r.i .+ c*nrow(ctc) for c in 0:pncol-1)...)],dims=3),dims=3)
    #     #     for r in eachrow(cond))

    #     # if plot
    #     #     for c in 1:pncol
    #     #         cmcys = Dict(condstring(r)=>dropdims(mean(ys[:,c,:,r.i],dims=3),dims=3) for r in eachrow(cond))
    #     #         cmccsd = Dict(condstring(r)=>dropdims(mean(csd(ys[:,c,:,r.i],h=hy),dims=3),dims=3) for r in eachrow(cond))
    #     #         for k in keys(cmcys)
    #     #             plotanalog(cmcys[k];fs,cunit=:uv,hy)
    #     #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_Column_$(c)_$(k)_MeanLFP$ext")),figfmt)
        
    #     #             plotanalog(imfilter(cmccsd[k],Kernel.gaussian((1,1)));fs,hy)
    #     #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_Column_$(c)_$(k)_MeanCSD$ext")),figfmt)
    #     #         end
    #     #     end
    #     # end

    #     # LFP of same depth
    #     @views pys = ys[.!exchmask,:,:] # exclude and flat channels
    #     chpos = vcat(chpositionnp(pversion)[.!exchmask]...) # exclude and flat channel positions 
    #     chgi = [findall(chpos[:,2].==up) for up in unique(chpos[:,2])] # group channels with same depth
    #     lfps = Array{Float64}(undef,length(chgi),size(pys)[2:end]...)
    #     @views foreach(i->lfps[i,:,:] = mean(pys[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    #     @views clfp = Dict(condstring(r)=>
    #         dropdims(mean(lfps[:,:,r.i],dims=3),dims=3) # trial average
    #     for r in eachrow(cond))
        
    #     # CSD of same depth
    #     ccsd = Dict(k=>stfilter(csd(v,h=hy),temporaltype=:sub,ti=baseindex) for (k,v) in clfp)
    #     times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(pys,2))

    #     if plot
    #         lfplim = mapreduce(x->maximum(abs.(x)),max,values(clfp))
    #         csdlim = mapreduce(x->maximum(abs.(x)),max,values(ccsd))
    #         for k in keys(clfp)
    #             plotanalog(clfp[k];x=times,hy,clims=(-lfplim,lfplim))
    #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_LFP$ext")),figfmt)
    #             plotanalog(ccsd[k];x=times,color=:RdBu,hy,clims=(-csdlim,csdlim))
    #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dCSD$ext")),figfmt)
    #         end
    #     end
    #     jldsave(joinpath(resultdir,"lf$(ii).jld2");clfp,ccsd,times,depths,fs,baseindex,siteid,exenv)


    #     # power spectrum of same depth
    #     nw = 2
    #     @views ppss,psfreqs = powerspectrum(pys[:,baseindex[end]+1:end,:],fs;freqrange=[1,100],nw)
    #     pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
    #     @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    #     tps = dropdims(mean(pss,dims=2),dims=2) # freq average
    #     @views cfps = Dict(condstring(r)=>
    #         dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average       
    #     for r in eachrow(cond))

    #     # local coherence
    #     @views lcs,lcfreqs = localcoherence(pys[:,baseindex[end]+1:end,:],chpos,fs;freqrange=[1,100],lr=55,sigma=25,chgroupdim=2,nw)
    #     tlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
    #     @views cflc = Dict(condstring(r)=>
    #         dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average       
    #     for r in eachrow(cond))

    #     if plot
    #         plotanalog(tps;hy,color=:heat,n=mean(tps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
    #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tPower.lf$ext")),figfmt)
    #         plotanalog(tlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(tlc,dims=2),layer)
    #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_tCoherence.lf$ext")),figfmt)
    #         fpslims = extrema(mapreduce(extrema,union,values(cfps)))
    #         flclims = extrema(mapreduce(extrema,union,values(cflc)))
    #         for k in keys(cfps)
    #             plotanalog(cfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(cfps[k],dims=2),layer)
    #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fPower.lf$ext")),figfmt)
    #             plotanalog(cflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(cflc[k],dims=2),layer)
    #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_fCoherence.lf$ext")),figfmt)
    #         end
    #     end
    #     jldsave(joinpath(resultdir,"lf$(ii)+.jld2");tps,cfps,psfreqs,tlc,cflc,lcfreqs,depths,fs,baseindex,siteid,exenv)
    # end

    # # Unit Spike Trian of Epochs
    # if plot
    #     epochext = max(preicidur,suficidur)
    #     epochs = [condon.-epochext condoff.+epochext]
    #     for i in eachindex(unitspike)
    #         ys,ns,ws,is = epochspiketrain(unitspike[i],ref2sync(epochs,dataset,unitsync[i]),isminzero=true,shift=epochext)
    #         title = "IMEC$(unitsync[i])_$(unitgood[i] ? "S" : "M")U$(unitid[i])_SpikeTrian"
    #         plotspiketrain(ys,timeline=[0,conddur],title=title)
    #         foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
    #     end
    # end
    
    # Single Unit Binary Spike Train of Condition Tests
    onsetshift = 300 # avoid transient onset response
    epoch = [onsetshift minconddur]
    epochs = condon.+epoch
    epochbins = epoch[1]:timetounit(1):epoch[2] # 1ms bin histogram convert spike times to binary sequence
    bst = map((st,si)->permutedims(psthspiketrains(epochspiketrain(st,ref2sync(epochs,dataset,si),isminzero=true,shift=-onsetshift).y,epochbins,israte=false,ismean=false).mat),
                unitspike[unitgood],unitsync[unitgood]) # nBin x nEpoch for each single unit
              
    # Single Unit Correlogram and Projection
    lag=100
    c = correlogramprojection(bst;lag,correction=(shufflejitter=true,l=25),maxprojlag=10,minbaselag=50,unitid=unitid[unitgood],condis=cond.i)

    if !isempty(c.ccgs)
        if plot
            for i in eachindex(c.ccgs)
                title = "Correlogram of SU$(c.ccgis[i][1]) and SU$(c.ccgis[i][2])"
                hline(c.ccgths[i];lw=0.6,color=:gray)
                bar!(c.x,c.ccgs[i];bar_width=1,legend=false,color=:gray15,linecolor=:match,title,xlabel="Lag (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,-50,-10,0,10,50,lag])
                foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
            end
        end
    end
    jldsave(joinpath(resultdir,"projection.jld2");c...,siteid,exenv)
end

function process_hartley_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files),spikesorter=param[:spikesorter])

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    datadir = joinpath(param[:dataroot],subject,siteid)
    siteresultdir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(siteresultdir,testid)
    isdir(resultdir) || mkpath(resultdir)

    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"]
    unitgood=spike["unitgood"];unitposition=spike["unitposition"];unitsync=spike["unitsync"]
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    jldsave(joinpath(resultdir,"spike.jld2");spike,siteid)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    exenv=Dict()
    exenv["ID"] = ex["ID"]
    exenv["conddur"] = conddur
    if haskey(ex,"Eye")
        exenv["eye"] = ex["Eye"]
    else
        exenv["eye"] = ex["Log"]
    end
    exenv["color"] = "$(exparam["ColorSpace"])_$(exparam["Color"])"

    # Prepare Conditions
    condtable = DataFrame(ex["Cond"])
    ctc = condtestcond(ex["CondTestCond"])
    replace!(ctc.SpatialFreq,NaN=>0) # convert old blank(SF=NaN) to new blank(SF=0)
    cond = condin(ctc)
    factors = finalfactor(ctc)
    blank = haskey(param,:blank) ? param[:blank] : :SpatialFreq=>0
    ci = ctc[!,blank.first].!=blank.second
    cctc = ctc[ci,factors]
    ccond = condin(cctc)
    ccondon = condon[ci]
    ccondoff = condoff[ci]
    ccondidx = condidx[ci]
    bi = .!ci
    bctc = ctc[bi,factors]
    bcondon = condon[bi]
    bcondoff = condoff[bi]
    bcondidx = condidx[bi]
    isblank = !isempty(bcondon)

    # Prepare Imageset
    bgcolor = RGBA(envparam["BGColor"]...)
    maxcolor = RGBA(envparam["MaxColor"]...)
    mincolor = RGBA(envparam["MinColor"]...)
    masktype = envparam["MaskType"]
    maskradius = envparam["MaskRadius"]
    masksigma = envparam["MaskSigma"]
    exenv["bgcolor"] = bgcolor
    exenv["maxcolor"] = maxcolor
    exenv["mincolor"] = mincolor
    exenv["masktype"] = masktype
    exenv["maskradius"] = maskradius
    exenv["masksigma"] = masksigma

    diameter = envparam["Diameter"]
    ppd = haskey(param,:ppd) ? param[:ppd] : 45 # pixel/deg, 4.5 pixel for SF=10 RF
    # diameter = 5
    # ii = round(Int,4.5ppd):round(Int,8.5ppd)
    # jj = range(round(Int,5.5ppd),length=length(ii))
    # d=round(length(ii)/ppd,digits=1)

    sizedeg = (diameter,diameter)
    imagesetname = splitext(splitdir(ex["CondPath"])[2])[1] * "_size$(sizedeg)_ppd$ppd" # hartley subspace, degree size and ppd define a unique image set
    if !haskey(param,imagesetname)
        # generate images of all conditions
        imageset = Dict{Any,Any}(:image => map(i->Gray.(grating(θ=deg2rad(i.Ori),sf=i.SpatialFreq,phase=i.SpatialPhase,size=sizedeg,ppd=ppd)),eachrow(condtable)))
        # imageset = Dict{Any,Any}(:image => map(i->Gray.(grating(θ=deg2rad(i.Ori),sf=i.SpatialFreq,phase=i.SpatialPhase,size=sizedeg,ppd=ppd)[ii,jj]),eachrow(condtable)))
        imageset[:sizepx] = size(imageset[:image][1])
        param[imagesetname] = imageset
    end
    # sizedeg=(d,d)

    # Prepare Stimuli
    imageset = param[imagesetname]
    bgcolor = oftype(imageset[:image][1][1],bgcolor)
    stimuliname = "bgcolor$(bgcolor)_masktype$(masktype)_maskradius$(maskradius)_masksigma$(masksigma)" # bgcolor and mask define a unique masking on the imageset
    if !haskey(imageset,stimuliname)
        imagestimuli = Dict{Any,Any}(:stimuli => map(i->alphablend.(alphamask(i,radius=maskradius,sigma=masksigma,masktype=masktype).y,bgcolor),imageset[:image]))
        imagestimuli[:unmaskindex] = alphamask(imageset[:image][1],radius=maskradius,sigma=masksigma,masktype=masktype).i
        imageset[stimuliname] = imagestimuli
    end
    imagestimuli = imageset[stimuliname]


    # imgfile = joinpath(param[:dataroot],"$(testid)_image.mat")
    # imgs = readmat(imgfile)["imageArray"]
    # mat2julia!(imgs)
    # px=613;py=135;rp=32
    #
    # rawimageset = map(i->dropdims(mean(i[py-rp:py+rp,px-rp:px+rp,:],dims=3),dims=3),imgs)
    # imageset = map(i->gaussian_pyramid(i, 1, 4, 1.5)[2],rawimageset)
    #
    # imagesize = size(imageset[1])
    # x = Array{Float64}(undef,length(ccondidx),prod(imagesize))
    # foreach(i->x[i,:]=vec(imageset[ccondidx[i]]),1:size(x,1))


    if :STA in param[:model]
        sizepx = imageset[:sizepx]
        xi = imagestimuli[:unmaskindex]
        xin = length(xi)
        ugood=Dict();uy=Dict();usta = Dict();uŷ=Dict();ugof=Dict()
        epochs = [condon condoff]
        delays = -40:5:180 # ms

        if isblank
            uci = unique(ccondidx)
            ucii = map(i->findall(condidx.==i),uci)
            buci = unique(bcondidx)
            bucii = map(i->findall(condidx.==i),buci)
            buciis = vcat(bucii...)
            # mean stimuli for all blank trials, hartley subspace stimuli usually have only one blank condition (gray)
            bx = vec(mean(mapreduce(i->gray.(imagestimuli[:stimuli][i][xi]),hcat,mapreduce((i,n)->fill(i,n),append!,buci,length.(bucii))),dims=2))
            x = Array{Float64}(undef,length(uci),xin)
            foreach(i->x[i,:]=gray.(imagestimuli[:stimuli][uci[i]][xi]).-bx,1:size(x,1)) # transform each stimuli to deviation from mean blank stimuli
            xcond=condtable[uci,:]

            for u in eachindex(unitspike)
                uid = unitid[u]
                ys = Array{Float64}(undef,length(delays),length(uci))
                ŷs = Array{Float64}(undef,length(delays),length(uci))
                stas = Array{Float64}(undef,length(delays),xin)
                for d in eachindex(delays)
                    depochs = ref2sync(epochs.+delays[d],dataset,unitsync[u])
                    y = epochspiketrainresponse_ono(unitspike[u],depochs,israte=true,isnan2zero=true)
                    by = mean(y[buciis]) # mean response of all blank trials
                    y = map(i->mean(y[i])-by,ucii) # transform each stimuli mean response to deviation from blank stimuli mean response
                    ys[d,:] = y
                    k = sta(x,y)
                    stas[d,:] = k
                    ŷs[d,:] = x*k
                end
                ugood[uid]=unitgood[u]
                uy[uid]=ys
                usta[uid]=stas
                uŷ[uid]=ŷs
                ugof[uid]=[goodnessoffit(ys[d,:],ŷs[d,:],k=xin) for d in eachindex(delays)]
            end
        else
            uci = unique(condidx)
            ucii = map(i->findall(condidx.==i),uci)
            x = Array{Float64}(undef,length(uci),xin)
            foreach(i->x[i,:]=gray.(imagestimuli[:stimuli][uci[i]][xi]),1:size(x,1))
            xcond=condtable[uci,:]

            for u in eachindex(unitspike)
                uid = unitid[u]
                ys = Array{Float64}(undef,length(delays),length(uci))
                ŷs = Array{Float64}(undef,length(delays),length(uci))
                stas = Array{Float64}(undef,length(delays),xin)
                for d in eachindex(delays)
                    depochs = ref2sync(epochs.+delays[d],dataset,unitsync[u])
                    y = epochspiketrainresponse_ono(unitspike[u],depochs,israte=true,isnan2zero=true)
                    y = map(i->mean(y[i]),ucii)
                    ys[d,:] = y
                    k = sta(x,y)
                    stas[d,:] = k
                    ŷs[d,:] = x*k
                end
                ugood[uid]=unitgood[u]
                uy[uid]=ys
                usta[uid]=stas
                uŷ[uid]=ŷs
                ugof[uid]=[goodnessoffit(ys[d,:],ŷs[d,:],k=xin) for d in eachindex(delays)]
            end
        end

        jldsave(joinpath(resultdir,"sta.jld2");sizepx,x,xi,xcond,ugood,uy,uŷ,usta,ugof,delays,siteid,sizedeg,exenv)
    end

    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : timetounit(40)
    batchunit = haskey(param,:batchunit) ? param[:batchunit] : nothing
    maxxsize = haskey(param,:maxxsize) ? param[:maxxsize] : 32
    if :ePPR in param[:model]
        sizepx = imageset[:sizepx]
        xi = imagestimuli[:unmaskindex]
        xscale = 255
        bg = xscale*gray(bgcolor)
        epochs = [condon condoff] .+ responsedelay
        mp = haskey(param,:epprparam) ? param[:epprparam] : (ndelay=1, nft=[3,3,3], lambda=320000)
        uroi=Dict();ugood=Dict();umodel=Dict();umodels=Dict();uhp=Dict()

        for u in eachindex(unitspike)
            # unitid[u]==630 || continue

            if isnothing(batchunit)
                roi=missing
            else
                id = "$(siteid)_$(unitgood[u] ? "S" : "M")U$(unitid[u])"
                i = findfirst(j->j==id,batchunit.id)
                isnothing(i) && continue
                roi = batchunit.roi[i]
            end
            idxrange,xsize,roi=takeroi(sizepx,ppd;roimaxresize=maxxsize,roi,issquare=true)
            x = Array{Float64}(undef,nrow(ctc),prod(xsize))
            foreach(i->x[i,:]=xscale*gray.(imresize_antialiasing(imagestimuli[:stimuli][condidx[i]][idxrange...],xsize)),1:size(x,1))

            sepochs = ref2sync(epochs,dataset,unitsync[u])
            y = epochspiketrainresponse_ono(unitspike[u],sepochs,israte=true,isnan2zero=true)
            unitresultdir = joinpath(resultdir,"IMEC$(unitsync[u])$(unitgood[u] ? "S" : "M")U$(unitid[u])_ePPR_$responsedelay")
            rm(unitresultdir,force=true,recursive=true)
            mlog = ePPRLog(debug=true,plot=true,dir=unitresultdir)
            mhp = ePPRHyperParams(xsize...;blankcolor=bg,mp.ndelay,mp.nft,mp.lambda)
            model,models = epprcv(x,y,mhp,mlog)

            if true && !isnothing(model)
                mlog(plotmodel(model,mhp),file="Model_Final (λ=$(mhp.lambda)).png")
            end

            uroi[unitid[u]]=roi
            ugood[unitid[u]]=unitgood[u]
            umodel[unitid[u]] = isnothing(model) ? model : clean!(model)
            umodels[unitid[u]] = clean!.(models)
            uhp[unitid[u]] = mhp
        end
        jldsave(joinpath(resultdir,"eppr.jld2");umodel,umodels,uhp,responsedelay,siteid,ugood,uroi,exenv)
    end
end

function process_condtest_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files),spikesorter=param[:spikesorter])

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    datadir = joinpath(param[:dataroot],subject,siteid)
    siteresultdir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(siteresultdir,testid)
    isdir(resultdir) || mkpath(resultdir)

    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"]
    unitgood=spike["unitgood"];unitposition=spike["unitposition"];unitsync=spike["unitsync"]
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    jldsave(joinpath(resultdir,"spike.jld2");spike,siteid)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff.-condon)
    exenv=Dict()
    exenv["ID"] = ex["ID"]
    exenv["conddur"] = conddur
    if haskey(ex,"Eye")
        exenv["eye"] = ex["Eye"]
    else
        exenv["eye"] = ex["Log"]
    end
    exenv["color"] = "$(exparam["ColorSpace"])_$(exparam["Color"])"

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    
    blank = :SpatialFreq=>0
    if ex["ID"] == "Color"
        # factors = [:HueAngle]
        blank = :Color=>36
        factors = [:Angle]
    elseif ex["ID"] == "OriSF"
        exenv["minmaxcolor"] = (mincolor = RGBA(envparam["MinColor"]...), maxcolor = RGBA(envparam["MaxColor"]...))
        exenv["TemporalFreq"] = getparam(envparam,"TemporalFreq")
        exenv["GratingType"] = getparam(envparam,"GratingType")
    end
    haskey(param,:blank) && (blank = param[:blank])

    ci = ctc[!,blank.first].!=blank.second
    cctc = ctc[ci,factors]
    ccond = condin(cctc)
    ccondon=condon[ci]
    ccondoff=condoff[ci]
    bi = .!ci
    bctc = ctc[bi,factors]
    bcondon=condon[bi]
    bcondoff=condoff[bi]
    isblank = !isempty(bcondon)


    # for ii in dataset["imecindex"]
    #     # Prepare AP
    #     d = "ap$ii"
    #     meta = dataset[d]["meta"]
    #     binfile = spike["datapath"]
    #     chmapraw = spike["chmapraw"]
    #     nch = spike["nch"]
    #     nsample = spike["nsample"]
    #     hx,hy,hz = meta["probespacing"]
    #     t0 = spike["t0"]
    #     winvraw = spike["winvraw"]
    #     exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
    #     pnrow,pncol = size(exchmask)
    #     depths = hy*(0:(pnrow-1))
    #     mmbinfile = mmap(binfile,Matrix{Int16},(nch,nsample),0)
    #     syncccondon = ref2sync(ccondon,dataset,ii)
    #     exenv["hy"] = hy
    #     exenv["t0"] = t0
    #     exenv["synccondon"] = syncccondon
    #     exenv["exchmask"] = exchmask
    #     exenv["chmapraw"] = chmapraw
    #     exenv["nch"] = nch
    #     exenv["nsample"] = nsample
    #     exenv["winvraw"] = winvraw
    #     exenv["fs"] = spike["fs"]
    #     exenv["binfile"] = binfile

    #     # Depth AP RMS
    #     epoch = [-40 150]
    #     epochs = syncccondon.+epoch
    #     # All AP epochs(mapped to concat file), unwhiten, gain corrected(voltage), bandpass filtered,
    #     # with all channels mapped in the shape of probe where excluded channels are replaced with rand local average
    #     ys = fill2mask(epochsamplenp(mmbinfile,spike["fs"],epochs.+t0,1:nch;meta,bandpass=[300,3000],whiten=winvraw),exchmask,chmap=chmapraw,randreplace=true)
    #     # 3 fold downsampling(10kHz)
    #     ys = resample(ys,1//3,dims=3)
    #     fs = spike["fs"]/3
    #     baseindex = epoch2sampleindex([0 50],fs)

    #     # RMS of combined columns
    #     @views ys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    #     @views prms = [rms(ys[i,baseindex[end]+1:end,j])^2 for i in 1:size(ys,1),j in 1:size(ys,3)]
    #     @views pc,freqs = coherencespectrum(ys[:,baseindex[end]+1:end,:],fs,freqrange=[300,3000],ismean=true) # freq average
    #     @views pc = hcat(map(i->bandmean(pc[:,:,i],r=5,s=1),1:size(pc,3))...) # local average with gaussian weights of -100:100μm, σ=20μm
    #     @views pbc,freqs = coherencespectrum(ys[:,baseindex,:],fs,freqrange=[300,3000],ismean=true)
    #     @views pbc = hcat(map(i->bandmean(pbc[:,:,i],r=5,s=1),1:size(pbc,3))...)
    #     pdc = abs.(pc.-pbc)
    #     @views crms = Dict(condstring(r)=>
    #         stfilter(dropdims(mean(mapwindow(x->rms(x)^2,ys[:,:,vcat((r.i .+ c*nrow(cctc) for c in 0:pncol-1)...)],(1,101,1),border="symmetric"),
    #                     dims=3),dims=3),temporaltype=:rcb,ti=baseindex)
    #         for r in eachrow(ccond)) # 10ms rms window
    #     times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(ys,2))
    #     if plot
    #         plotanalog(prms;hy,color=:heat,n=mean(prms,dims=2),xlabel="Trial",xunit="",cunit=:fr)
    #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_RMS$ext")),figfmt)
    #         plotanalog(pc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(pc,dims=2))
    #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_Coherence$ext")),figfmt)
    #         plotanalog(pdc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(pdc,dims=2))
    #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_dCoherence$ext")),figfmt)
    #         clim = mapreduce(x->maximum(abs.(x)),max,values(crms))
    #         for k in keys(crms)
    #             plotanalog(crms[k];x=times,hy,clims=(-clim,clim))
    #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dRMS$ext")),figfmt)
    #         end
    #     end
    #     jldsave(joinpath(resultdir,"ap$(ii).jld2");prms,pc,pbc,crms,times,depths,fs,siteid,exenv)

    #     # # All AP epochs(mapped to concat file), remain whiten, bandpass filtered,
    #     # # with all channels mapped in the shape of probe where excluded channels are replaced with rand local average
    #     # ys = fill2mask(epochsamplenp(mmbinfile,spike["fs"],epochs.+t0,1:nch;meta=[],bandpass=[300,3000],whiten=nothing),exchmask,chmap=chmapraw,randreplace=true)
    #     # # 3 fold downsampling(10kHz)
    #     # ys = resample(ys,1//3,dims=3)
    #     # fs = spike["fs"]/3

    #     # # RMS of combined columns
    #     # @views ys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    #     # @views wprms = [rms(ys[i,baseindex[end]+1:end,j])^2 for i in 1:size(ys,1),j in 1:size(ys,3)]
    #     # @views wpc,freqs = coherencespectrum(ys[:,baseindex[end]+1:end,:],fs,freqrange=[300,3000],ismean=true) # freq average
    #     # @views wpc = hcat(map(i->bandmean(wpc[:,:,i],r=5,s=1),1:size(wpc,3))...) # local average with gaussian weights of -100:100μm, σ=20μm
    #     # @views wpbc,freqs = coherencespectrum(ys[:,baseindex,:],fs,freqrange=[300,3000],ismean=true)
    #     # @views wpbc = hcat(map(i->bandmean(wpbc[:,:,i],r=5,s=1),1:size(wpbc,3))...)
    #     # wpdc = abs.(wpc.-wpbc)
    #     # @views wcrms = Dict(condstring(r)=>
    #     #     stfilter(dropdims(mean(mapwindow(x->rms(x)^2,ys[:,:,vcat((r.i .+ c*nrow(cctc) for c in 0:pncol-1)...)],(1,101,1),border="symmetric"),
    #     #                 dims=3),dims=3),temporaltype=:rcb,ti=baseindex)
    #     #     for r in eachrow(ccond)) # 10ms rms window
    #     # if plot
    #     #     plotanalog(wprms;hy,color=:heat,n=mean(wprms,dims=2),xlabel="Trial",xunit="",cunit=:fr)
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_wRMS$ext")),figfmt)
    #     #     plotanalog(wpc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(wpc,dims=2))
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_wCoherence$ext")),figfmt)
    #     #     plotanalog(wpdc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(wpdc,dims=2))
    #     #     foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_wdCoherence$ext")),figfmt)
    #     #     clim = mapreduce(x->maximum(abs.(x)),max,values(wcrms))
    #     #     for k in keys(wcrms)
    #     #         plotanalog(wcrms[k];x=times,hy,clims=(-clim,clim))
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_wdRMS$ext")),figfmt)
    #     #     end
    #     # end
    #     # jldsave(joinpath(resultdir,"ap$(ii).jld2");prms,pc,pbc,crms,wprms,wpc,wpbc,wcrms,times,depths,fs,siteid,exenv)

    #     # # Depth Unit PSTH
    #     # epoch = [-40 150]
    #     # epochs = syncccondon.+epoch
    #     # bw = 2
    #     # psthbins = epoch[1]:bw:epoch[2]
    #     # baseindex = epoch2sampleindex([0 50],1/(bw*SecondPerUnit))
    #     # # All Unit
    #     # ui = unitsync.==ii
    #     # @views unitepochpsth = map(st->psthspiketrains(epochspiketrain(st,epochs,isminzero=true,shift=-epoch[1]).y,psthbins,israte=true,ismean=false),unitspike[ui])
    #     # @views cpsth = Dict(condstring(r)=>begin
    #     #                     p = spacepsth(map(p->(;vmeanse(p.mat[r.i,:]).m,p.x),unitepochpsth),unitposition[ui,:],
    #     #                         w=replace(unitgood[ui],0=>1.2),spacerange=depths,bw=2hy,step=hy) # multi-unit count = 1.2 for unit density
    #     #                     (;psth=stfilter(mapwindow(mean,p.psth,(1,5),border="symmetric"),temporaltype=:sub,ti=baseindex),p.x,p.y,p.n) # 10ms mean window
    #     #                 end  for r in eachrow(ccond))
    #     # if plot
    #     #     clim = mapreduce(i->maximum(abs.(i.psth)),max,values(cpsth))
    #     #     for k in keys(cpsth)
    #     #         plotanalog(cpsth[k].psth;cpsth[k].x,cpsth[k].y,clims=(-clim,clim))
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_All-Unit_$(k)_dPSTH$ext")),figfmt)
    #     #     end
    #     # end
    #     # jldsave(joinpath(resultdir,"depthpsth$(ii).jld2");cpsth,siteid,exenv)


    #     # Prepare LFP
    #     d = "lf$ii"
    #     meta = dataset[d]["meta"]
    #     binfile = meta["fileName"]
    #     # lfbin = matchfile(Regex("^$(testid)\\w*.imec.lf.bin"),dir = datadir,join=true)[1]
    #     nsavedch = meta["nSavedChans"]
    #     nsample = meta["nFileSamp"]
    #     nch = nsavedch-1 # exclude sync channel
    #     hx,hy,hz = meta["probespacing"]
    #     t0 = 0 # not concat drift correct binary yet
    #     exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
    #     pnrow,pncol = size(exchmask)
    #     depths = hy*(0:(pnrow-1))
    #     mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)
    #     synccondon = ref2sync(condon,dataset,ii)
    #     exenv["hy"] = hy
    #     exenv["t0"] = t0
    #     exenv["synccondon"] = synccondon
    #     exenv["exchmask"] = exchmask
    #     exenv["nch"] = nch
    #     exenv["nsample"] = nsample
    #     exenv["fs"] = meta["fs"]
    #     exenv["binfile"] = binfile

    #     # Depth LFP and CSD
    #     epoch = [-40 150]
    #     epochs = synccondon.+epoch
    #     # All LFP epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered,
    #     # with all channels mapped in the shape of probe where excluded channels are replaced with local average
    #     ys = fill2mask(epochsamplenp(mmbinfile,meta["fs"],epochs.+t0,1:nch;meta,bandpass=[1,100]),exchmask)
    #     # 2.5 fold downsampling(1kHz)
    #     ys = resample(ys,1/2.5,dims=3)
    #     fs = meta["fs"]/2.5
    #     baseindex = epoch2sampleindex([0 50],fs)

    #     # LFP and CSD of combined columns
    #     @views ys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    #     @views clfp = Dict(condstring(r)=>
    #         dropdims(mean(ys[:,:,vcat((r.i .+ c*nrow(ctc) for c in 0:pncol-1)...)],dims=3),dims=3)
    #         for r in eachrow(cond))
    #     ccsd = Dict(k=>stfilter(csd(v,h=hy),temporaltype=:sub,ti=baseindex) for (k,v) in clfp)
    #     times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(ys,2))
    #     if plot
    #         lfplim = mapreduce(x->maximum(abs.(x)),max,values(clfp))
    #         csdlim = mapreduce(x->maximum(abs.(x)),max,values(ccsd))
    #         for k in keys(clfp)
    #             plotanalog(clfp[k];x=times,hy,clims=(-lfplim,lfplim))
    #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_LFP$ext")),figfmt)
    #             plotanalog(ccsd[k];x=times,color=:RdBu,hy,clims=(-csdlim,csdlim))
    #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dCSD$ext")),figfmt)
    #         end
    #     end
    #     jldsave(joinpath(resultdir,"lf$(ii).jld2");clfp,ccsd,times,depths,fs,baseindex,siteid,exenv)

    #     # # Depth Power Spectrum
    #     # epochdur = timetounit(preicidur)
    #     # epoch = [0 epochdur]
    #     # epochs = ref2sync(ccondon.+epoch,dataset,ii)
    #     # nw = 2
    #     # ys = fill2mask(epochsamplenp(mmbinfile,meta["fs"],epochs.+t0,1:nch;meta,bandpass=[1,100]),exchmask)
    #     # ys = resample(ys,1/2.5,dims=3)
    #     # fs = meta["fs"]/2.5
    #     # @views ys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    #     # ps,freqs = powerspectrum(ys,fs;freqrange=[1,100],nw)

    #     # epoch = [-epochdur 0]
    #     # epochs = ref2sync(ccondon.+epoch,dataset,ii)
    #     # ys = fill2mask(epochsamplenp(mmbinfile,meta["fs"],epochs.+t0,1:nch;meta,bandpass=[1,100]),exchmask)
    #     # ys = resample(ys,1/2.5,dims=3)
    #     # @views ys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    #     # bps,freqs = powerspectrum(ys,fs;freqrange=[1,100],nw)

    #     # pc = log2.(ps./bps)
    #     # @views cpc = Dict(condstring(r)=>
    #     #     dropdims(mean(pc[:,:,vcat((r.i .+ c*nrow(ctc) for c in 0:pncol-1)...)],dims=3),dims=3)
    #     #     for r in eachrow(cond))
    #     # if plot
    #     #     pclim = mapreduce(x->maximum(abs.(x)),max,values(cpc))
    #     #     for k in keys(cpc)
    #     #         plotanalog(cpc[k];x=freqs,hy,xlabel="Frequency",xunit="Hz",clims=(-pclim,pclim),color=:vik)
    #     #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_PowerContrast$ext")),figfmt)
    #     #     end
    #     # end
    #     # jldsave(joinpath(resultdir,"lf$(ii).jld2");clfp,ccsd,cpc,times,depths,freqs,fs,baseindex,siteid,exenv)
    # end

    # # Unit Spike Trian of Epochs
    # if plot
    #     epochext = max(preicidur,suficidur)
    #     epochs = [condon.-epochext condoff.+epochext]
    #     for i in eachindex(unitspike)
    #         ys,ns,ws,is = epochspiketrain(unitspike[i],ref2sync(epochs,dataset,unitsync[i]),isminzero=true,shift=epochext)
    #         title = "IMEC$(unitsync[i])_$(unitgood[i] ? "S" : "M")U$(unitid[i])_SpikeTrian"
    #         plotspiketrain(ys,timeline=[0,conddur],title=title)
    #         foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
    #     end
    # end

    # @goto CORR

    # Condition Response
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : timetounit(15)
    minresponse = haskey(param,:minresponse) ? param[:minresponse] : 3
    rs=[];prs=[];srs=[];brs=[];responsive=[];modulative=[]
    cepochs = [ccondon ccondoff] .+ responsedelay
    pcepochs = [ccondon.-preicidur ccondon] .+ responsedelay
    scepochs = [ccondoff ccondoff.+suficidur] .+ responsedelay
    bepochs = [bcondon bcondoff] .+ responsedelay
    for i in eachindex(unitspike)
        r = epochspiketrainresponse_ono(unitspike[i],ref2sync(cepochs,dataset,unitsync[i]),israte=true)
        pr = epochspiketrainresponse_ono(unitspike[i],ref2sync(pcepochs,dataset,unitsync[i]),israte=true)
        sr = epochspiketrainresponse_ono(unitspike[i],ref2sync(scepochs,dataset,unitsync[i]),israte=true)
        push!(rs,r);push!(prs,pr);push!(srs,sr)
        push!(responsive,isresponsive(pr,r,ccond.i))
        push!(modulative,ismodulative(r,ccond.i))
        # if isblank
        #     br = subrvr(unitspike[i],ref2sync(bepochs,dataset,unitsync[i]))
        #     mseuc = condresponse(br,condin(bctc))
        #     push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
        # end
        if plot && length(factors)==2
            plotcondresponse(condresponse(r,ccond))
            foreach(ext->savefig(joinpath(resultdir,"IMEC$(unitsync[i])_$(unitgood[i] ? "S" : "M")U$(unitid[i])_CondResponse$ext")),figfmt)
        end
    end

    # Condition Response In Factor Space
    fi,fa = factorspace(ccond)
    frs = map(r->factorresponse(r,fi),rs)
    pfrs = map(r->factorresponse(r,fi),prs)
    sfrs = map(r->factorresponse(r,fi),srs)

    maxi=[];enoughresponse=[];frf = NamedTuple(f=>[] for f in keys(fa))
    fai = collect(zip(1:ndims(fi),keys(fa)))
    for i in eachindex(frs)
        oi = Any[Tuple(argmax(replace(frs[i].m, missing => -Inf)))...]
        push!(maxi,oi)
        push!(enoughresponse,frs[i].m[oi...] >= minresponse)
        u = "IMEC$(unitsync[i])$(unitgood[i] ? "S" : "M")U$(unitid[i])"
        # max response line for each factor dimemsion
        for (d,f) in fai
            li = deepcopy(oi)
            li[d]=(:)
            rf = (responsive[i] && enoughresponse[i]) ? factorresponsefeature(fa[f],frs[i].r[li...],fm=frs[i].m[li...],factor=f) : missing
            push!(frf[f],rf)

            if plot
                rdf=DataFrame(:m=>frs[i].m[li...],:se=>frs[i].se[li...],:u=>u,f=>fa[f])
                prdf=DataFrame(:m=>pfrs[i].m[li...],:se=>pfrs[i].se[li...],:u=>"Pre_$u",f=>fa[f])
                srdf=DataFrame(:m=>sfrs[i].m[li...],:se=>sfrs[i].se[li...],:u=>"Suf_$u",f=>fa[f])
                projection = iscircfactor(f) ? :polar : :none
                plotcondresponse(dropmissing([prdf;srdf;rdf]);color=[:gray85,:gray50,:gray5],linewidth=[1,2,3],grid=true,projection,response=isblank ? ubr[i:i] : [])
                if !ismissing(rf) && !isempty(rf.fit)
                    x = projection == :polar ? (0:0.02:2π) : range(extrema(fa[f])...,step=0.01*mean(diff(sort(fa[f]))))
                    plot!(x->rf.fit.mfit.fun(x,rf.fit.mfit.param),x,lw=2,color=:deepskyblue,label="$(rf.fit.mfit.model)_fit")
                end
                foreach(ext->savefig(joinpath(resultdir,"$(u)_$(f)_Tuning$ext")),figfmt)
            end
        end
    end

    # Condition PSTH
    psthbw = haskey(param,:psthbw) ? param[:psthbw] : 10
    epoch = [-preicidur conddur+suficidur]
    binedges = epoch[1]:psthbw:epoch[2]
    epochs = ccondon .+ epoch
    psths = map((st,si)->psthspiketrains(epochspiketrain(st,ref2sync(epochs,dataset,si),isminzero=true,shift=preicidur).y,binedges,israte=true,ismean=false),unitspike,unitsync)
    fpsth = (;psth=map(p->factorresponse(p.mat,fi),psths),first(psths).x)
    jldsave(joinpath(resultdir,"factorresponse.jld2");fi,fa,frs,pfrs,sfrs,responsive,modulative,enoughresponse,maxi,frf,fpsth,siteid,unitid,unitgood,exenv)
return
    @label CORR
    # Single Unit Binary Spike Train of Condition Tests
    onsetshift = 300 # avoid transient onset response
    epoch = [onsetshift minconddur]
    epochs = ccondon.+epoch
    epochbins = epoch[1]:timetounit(1):epoch[2] # 1ms bin histogram convert spike times to binary sequence
    bst = map((st,si)->permutedims(psthspiketrains(epochspiketrain(st,ref2sync(epochs,dataset,si),isminzero=true,shift=-onsetshift).y,epochbins,israte=false,ismean=false).mat),
                unitspike[unitgood],unitsync[unitgood]) # nBin x nEpoch for each single unit
              
    # Single Unit Correlogram and Projection
    lag=100
    c = correlogramprojection(bst;lag,correction=(shufflejitter=true,l=25),maxprojlag=10,minbaselag=50,unitid=unitid[unitgood],condis=ccond.i)

    if !isempty(c.ccgs)
        if plot
            for i in eachindex(c.ccgs)
                title = "Correlogram of SU$(c.ccgis[i][1]) and SU$(c.ccgis[i][2])"
                hline(c.ccgths[i];lw=0.6,color=:gray)
                bar!(c.x,c.ccgs[i];bar_width=1,legend=false,color=:gray15,linecolor=:match,title,xlabel="Lag (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,-50,-10,0,10,50,lag])
                foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
            end
        end
    end
    jldsave(joinpath(resultdir,"projection.jld2");c...,siteid,exenv)
end

function process_cycle_spikeglx(files,param;uuid="",log=nothing,plot=true)
end
