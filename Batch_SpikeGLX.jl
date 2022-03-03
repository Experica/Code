using StatsBase,StatsPlots,Images,DSP,ePPR

function process_flash_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

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
    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    # jldsave(joinpath(resultdir,"spike.jld2");spike,siteid)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff-condon)
    if haskey(ex,"Eye")
        eye = ex["Eye"]
    else
        eye = ex["Log"]
    end

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)

    for ii in dataset["imecindex"]
        # Prepare AP
        d = "ap$ii"
        meta = dataset[d]["meta"]
        binfile = meta["fileName"]
        nsavedch = meta["nSavedChans"]
        nsample = meta["nFileSamp"]
        fs = meta["fs"]
        nch = nsavedch-1 # exclude sync channel
        hx,hy,hz = meta["probespacing"]
        exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
        pnrow,pncol = size(exchmask)
        depths = hy*(0:(pnrow-1))
        mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)

        # Depth AP RMS
        # epochdur = timetounit(150)
        # epoch = [0 epochdur]
        # epochs = ref2sync(condon.+epoch,dataset,ii)
        # # All AP epochs, gain corrected(voltage), bandpass filtered,
        # # in the same shape of probe where excluded channels are replaced with local average
        # ys = reshape2mask(epochsamplenp(mmbinfile,fs,epochs,1:nch,meta=meta,bandpass=[300,3000]),exchmask)
        # # 3 fold downsampling(10kHz)
        # ys = ys[:,:,1:3:end,:]
        # fs /= 3
        #
        # # ΔRMS of combined columns
        # @views pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        # basedur = timetounit(15)
        # baseindex = epoch2sampleindex([0 basedur],fs)
        # @views cdrms = Dict(condstring(r)=>
        #     stfilter(dropdims(mean(mapwindow(rms,pys[:,:,[r.i;r.i.+nrow(ctc)]],(1,101,1),border="symmetric"),dims=3),dims=3),temporaltype=:sub,ti=baseindex)
        #     for r in eachrow(cond)) # 10ms rms window
        # time = (1/fs/SecondPerUnit)*(0:(size(pys,2)-1))
        # if plot
        #     ylim = mapreduce(y->maximum(abs.(y)),max,values(cdrms))
        #     for k in keys(cdrms)
        #         plotanalog(cdrms[k];fs,hy,clims=(-ylim,ylim))
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dMeanRMS$ext")),figfmt)
        #     end
        # end
        # jldsave(joinpath(resultdir,"ap$(ii).jld2");cdrms,time,depths,fs,siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")



        # Depth Unit PSTH
        # epochdur = timetounit(150)
        # epoch = [0 epochdur]
        # epochs = ref2sync(condon.+epoch,dataset,ii)
        # bw = timetounit(2)
        # psthbins = epoch[1]:bw:epoch[2]
        # basedur = timetounit(15)
        # baseindex = epoch2sampleindex([0 basedur],1/(bw*SecondPerUnit))
        # # All Unit
        # ui = unitsync.==ii
        # @views unitepochpsth = map(st->psthspiketrains(epochspiketrain(st,epochs,isminzero=true).y,psthbins,israte=true,ismean=false),unitspike[ui])
        # @views cdpsth = Dict(condstring(r)=>begin
        #                     p = spacepsth(map(p->(;vmeanse(p.mat[r.i,:]).m,p.x),unitepochpsth),unitposition[ui,:],
        #                         w=replace(unitgood[ui],0=>1.2),lim=(0,hy*(pnrow-1)),bw=2hy,step=hy) # multi-unit count = 1.2 for unit density
        #                     (;psth=stfilter(mapwindow(mean,p.psth,(1,5),border="symmetric"),temporaltype=:sub,ti=baseindex),p.x,p.y,p.n) # 10ms mean window
        #                 end  for r in eachrow(cond))
        # if plot
        #     lim = mapreduce(i->maximum(abs.(i.psth)),max,values(cdpsth))
        #     for k in keys(cdpsth)
        #         plotanalog(cdpsth[k].psth;cdpsth[k].x,cdpsth[k].y,clims=(-lim,lim))
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_All-Unit_$(k)_dPSTH$ext")),figfmt)
        #     end
        # end
        # jldsave(joinpath(resultdir,"depthpsth$(ii).jld2");cdpsth,siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")


        # Prepare LFP
        d = "lf$ii"
        meta = dataset[d]["meta"]
        binfile = meta["fileName"]
        # lfbin = matchfile(Regex("^$(testid)\\w*.imec.lf.bin"),dir = datadir,join=true)[1]
        nsavedch = meta["nSavedChans"]
        nsample = meta["nFileSamp"]
        fs = meta["fs"]
        nch = nsavedch-1 # exclude sync channel
        hx,hy,hz = meta["probespacing"]
        exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
        pnrow,pncol = size(exchmask)
        depths = hy*(0:(pnrow-1))
        mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)

        # Depth LFP and CSD
        epochdur = timetounit(150)
        epoch = [0 epochdur]
        epochs = ref2sync(condon.+epoch,dataset,ii)
        # All LFP epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered,
        # in the same shape of probe where excluded channels are replaced with local average
        ys = reshape2mask(epochsamplenp(mmbinfile,fs,epochs,1:nch,meta=meta,bandpass=[1,100]),exchmask)

        # if plot
        #     for c in 1:pncol
        #         cmcys = Dict(condstring(r)=>dropdims(mean(ys[:,c,:,r.i],dims=3),dims=3) for r in eachrow(cond))
        #         cmccsd = Dict(condstring(r)=>dropdims(mean(csd(ys[:,c,:,r.i],h=hy),dims=3),dims=3) for r in eachrow(cond))
        #         for k in keys(cmcys)
        #             plotanalog(cmcys[k];fs,cunit=:uv,hy)
        #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_Column_$(c)_$(k)_MeanLFP$ext")),figfmt)
        #
        #             plotanalog(imfilter(cmccsd[k],Kernel.gaussian((1,1)));fs,hy)
        #             foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_Column_$(c)_$(k)_MeanCSD$ext")),figfmt)
        #         end
        #     end
        # end

        # LFP and ΔCSD of combined columns
        @views pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        cmlfp = Dict(condstring(r)=>dropdims(mean(pys[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
        basedur = timetounit(15)
        baseindex = epoch2sampleindex([0 basedur],fs)
        cdcsd = Dict(k=>stfilter(csd(v,h=hy),temporaltype=:sub,ti=baseindex) for (k,v) in cmlfp)
        time = (1/fs/SecondPerUnit)*(0:(size(pys,2)-1))
        if plot
            llim = mapreduce(y->maximum(abs.(y)),max,values(cmlfp))
            clim = mapreduce(y->maximum(abs.(y)),max,values(cdcsd))
            for k in keys(cmlfp)
                plotanalog(cmlfp[k];fs,hy,clims=(-llim,llim))
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_LFP$ext")),figfmt)

                plotanalog(cdcsd[k];fs,color=:RdBu,hy,clims=(-clim,clim))
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dCSD$ext")),figfmt)
            end
        end
        jldsave(joinpath(resultdir,"lfp$(ii).jld2");cmlfp,cdcsd,time,depths,fs,siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")
        return

        # Depth Power Spectrum
        epochdur = timetounit(300)
        epoch = [0 epochdur]
        epochs = ref2sync(condon.+epoch,dataset,ii)
        nw = 2
        ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
        pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)

        epoch = [-epochdur 0]
        epochs = ref2sync(condon.+epoch,dataset,ii)
        ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
        pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        bps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)

        pcs = ps./bps.-1
        cmpcs = Dict(condstring(r)=>dropdims(mean(pcs[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
        if plot
            fcmpcs = Dict(k=>imfilter(cmpcs[k],Kernel.gaussian((1,1))) for k in keys(cmpcs))
            lim = mapreduce(pc->maximum(abs.(pc)),max,values(fcmpcs))
            for k in keys(fcmpcs)
                plotanalog(fcmpcs[k],x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],clims=(-lim,lim),color=:vik)
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_PowerContrast$ext")),figfmt)
            end
        end
        save(joinpath(resultdir,"lfp$(ii).jld2"),"cmlfp",cmys,"cmcsd",cmdcsd,"cmpc",cmpcs,"freq",freq,"time",x,"depth",depths,"fs",fs,
        "siteid",siteid,"eye",eye,"color","$(exparam["ColorSpace"])_$(exparam["Color"])")

    end

    # Unit Spike Trian of Condition Epochs
    if plot
        epochext = max(preicidur,suficidur)
        epochs = [condon.-epochext condoff.+epochext]
        for u in eachindex(unitspike)
            ys,ns,ws,is = epochspiketrain(unitspike[u],ref2sync(epochs,dataset,unitsync[u]),isminzero=true,shift=epochext)
            title = "$(ugs[u])Unit_$(unitid[u])_SpikeTrian"
            plotspiketrain(ys,timeline=[0,conddur],title=title)
            foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
        end
    end

    # Single Unit Binary Spike Train of Condition Tests
    bepochext = timetounit(-300)
    bepoch = [-bepochext minconddur]
    bepochs = condon.+bepoch
    spikebins = bepoch[1]:timetounit(1):bepoch[2] # 1ms bin histogram will convert spike times to binary spike train
    bst = map((st,si)->permutedims(psthspiketrains(epochspiketrain(st,ref2sync(bepochs,dataset,si),isminzero=true,shift=bepochext).y,spikebins,israte=false,ismean=false).mat),unitspike[unitgood],unitsync[unitgood]) # nSpikeBin x nEpoch for each unit
    # Single Unit Correlogram and Circuit
    lag=50
    ccgs,x,ccgis,projs,eunits,iunits,projweights = circuitestimate(bst,lag=lag,unitid=unitid[unitgood],condis=cond.i)

    if !isempty(projs)
        if plot
            for i in eachindex(ccgs)
                title = "Correlogram between single unit $(ccgis[i][1]) and $(ccgis[i][2])"
                bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
                foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
            end
            plotcircuit(unitposition,unitid,projs,unitgood=unitgood,eunits=eunits,iunits=iunits,layer=layer)
            foreach(ext->savefig(joinpath(resultdir,"UnitPosition_Circuit$ext")),figfmt)
        end
        jldsave(joinpath(resultdir,"circuit.jld2");projs,eunits,iunits,projweights,siteid)
    end
end

function process_hartley_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

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
    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    # jldsave(joinpath(resultdir,"spike.jld2");spike,siteid)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    if haskey(ex,"Eye")
        eye = ex["Eye"]
    else
        eye = ex["Log"]
    end

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
    diameter = envparam["Diameter"]
    ppd = haskey(param,:ppd) ? param[:ppd] : 45 # pixel/deg, 4.5 pixel for a SF=0.1 RF
    # diameter = 5
    # ii = round(Int,4.5ppd):round(Int,8.5ppd)
    # jj = range(round(Int,5.5ppd),length=length(ii))
    # d=round(length(ii)/ppd,digits=1)

    sizedeg = (diameter,diameter)
    imagesetname = splitext(splitdir(ex["CondPath"])[2])[1] * "_size$(sizedeg)_ppd$ppd" # hartley subspace, degree size and ppd define a unique image set
    if !haskey(param,imagesetname)
        imageset = Dict{Any,Any}(:image => map(i->Gray.(grating(θ=deg2rad(i.Ori),sf=i.SpatialFreq,phase=i.SpatialPhase,size=sizedeg,ppd=ppd)),eachrow(condtable)))
        # imageset = Dict{Any,Any}(:image => map(i->Gray.(grating(θ=deg2rad(i.Ori),sf=i.SpatialFreq,phase=i.SpatialPhase,size=sizedeg,ppd=ppd)[ii,jj]),eachrow(condtable)))
        imageset[:sizepx] = size(imageset[:image][1])
        param[imagesetname] = imageset
    end
    # sizedeg=(d,d)

    # Prepare Image Stimuli
    imageset = param[imagesetname]
    bgcolor = oftype(imageset[:image][1][1],bgcolor)
    imagestimuliname = "bgcolor$(bgcolor)_masktype$(masktype)_maskradius$(maskradius)_masksigma$(masksigma)" # bgcolor and mask define a unique masking on an image set
    if !haskey(imageset,imagestimuliname)
        imagestimuli = Dict{Any,Any}(:stimuli => map(i->alphablend.(alphamask(i,radius=maskradius,sigma=masksigma,masktype=masktype).y,bgcolor),imageset[:image]))
        imagestimuli[:unmaskindex] = alphamask(imageset[:image][1],radius=maskradius,sigma=masksigma,masktype=masktype).i
        imageset[imagestimuliname] = imagestimuli
    end
    imagestimuli = imageset[imagestimuliname]


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
        uy=Dict();usta = Dict();uŷ=Dict();ugof=Dict()
        epochs = [condon condoff]
        delays = -40:5:200

        if isblank
            uci = unique(ccondidx)
            ucii = map(i->findall(condidx.==i),uci)
            buci = unique(bcondidx)
            bucii = mapreduce(i->findall(condidx.==i),append!,buci)
            bx = vec(mean(mapreduce(i->gray.(imagestimuli[:stimuli][i][xi]),hcat,buci),dims=2))
            x = Array{Float64}(undef,length(uci),xin)
            foreach(i->x[i,:]=gray.(imagestimuli[:stimuli][uci[i]][xi]).-bx,1:size(x,1))
            xcond=condtable[uci,:]

            for u in eachindex(unitspike)
                unitgood[u] || continue # only for single unit
                ys = Array{Float64}(undef,length(delays),length(uci))
                ŷs = Array{Float64}(undef,length(delays),length(uci))
                stas = Array{Float64}(undef,length(delays),xin)
                for d in eachindex(delays)
                    depochs = ref2sync(epochs.+delays[d],dataset,unitsync[u])
                    y = epochspiketrainresponse_ono(unitspike[u],depochs,israte=true,isnan2zero=true)
                    by = mean(y[bucii])
                    y = map(i->mean(y[i])-by,ucii)
                    ys[d,:] = y
                    k = sta(x,y)
                    stas[d,:] = k
                    ŷs[d,:] = x*k
                end
                uy[unitid[u]]=ys
                usta[unitid[u]]=stas
                uŷ[unitid[u]]=ŷs
                ugof[unitid[u]]=[goodnessoffit(ys[d,:],ŷs[d,:],k=xin) for d in eachindex(delays)]
            end
        else
            uci = unique(condidx)
            ucii = map(i->findall(condidx.==i),uci)
            x = Array{Float64}(undef,length(uci),xin)
            foreach(i->x[i,:]=gray.(imagestimuli[:stimuli][uci[i]][xi]),1:size(x,1))
            xcond=condtable[uci,:]

            for u in eachindex(unitspike)
                unitgood[u] || continue # only for single unit
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
                uy[unitid[u]]=ys
                usta[unitid[u]]=stas
                uŷ[unitid[u]]=ŷs
                ugof[unitid[u]]=[goodnessoffit(ys[d,:],ŷs[d,:],k=xin) for d in eachindex(delays)]
            end
        end

        jldsave(joinpath(resultdir,"sta.jld2");sizepx,x,xi,xcond,uy,uŷ,usta,ugof,delays,
        siteid,sizedeg,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])",maxcolor,mincolor)
    end

    if :ePPR in param[:model]
        sizepx = imageset[:sizepx]
        xi = imagestimuli[:unmaskindex]
        cell = haskey(param,:cell) ? param[:cell] : nothing
        roir = haskey(param,:roir) ? param[:roir] : 16
        xscale=255
        bg = xscale*gray(bgcolor)
        d = 60
        umodel=Dict();umodels=Dict();uhp=Dict()

        for u in eachindex(unitspike)
            unitgood[u] || continue
            # unitid[u]!=127 && continue

            if isnothing(cell)
                c = round.(Int,sizepx./2)
                r = minimum(c)
                r = clamproi(c,r,sizepx)
            else
                id = "$(siteid)_SU$(unitid[u])"
                i = findfirst(cell.id.==id)
                if isnothing(i) || ismissing(cell.roicenter[i])
                    c = round.(Int,sizepx./2)
                    r = minimum(c)
                    r = clamproi(c,r,sizepx)
                else
                    c = round.(Int,cell.roicenter[i]*ppd)
                    r = round(Int,cell.roiradius[i]*ppd)
                    r = clamproi(c,r,sizepx)
                end
            end
            xr = min(r,roir)
            xsize = (2xr+1,2xr+1)
            x = Array{Float64}(undef,nrow(ctc),prod(xsize))
            foreach(i->x[i,:]=xscale*gray.(imresize_antialiasing(imagestimuli[:stimuli][condidx[i]][map(i->(-r:r).+i,c)...],xsize)),1:size(x,1))

            y = epochspiketrainresponse_ono(unitspike[u],condon.+d,condoff.+d,israte=true,isnan2zero=true)
            unitresultdir = joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_ePPR_$d")
            rm(unitresultdir,force=true,recursive=true)
            log = ePPRLog(debug=true,plot=plot,dir=unitresultdir)
            hp = ePPRHyperParams(xsize...,blankcolor=bg,ndelay=param[:eppr_ndelay],nft=param[:eppr_nft],lambda=param[:eppr_lambda])
            model,models = epprcv(x,y,hp,log)

            if plot && !isnothing(model)
                log(plotmodel(model,hp),file="Model_Final (λ=$(hp.lambda)).png")
            end

            umodel[unitid[u]] = isnothing(model) ? model : clean!(model)
            umodels[unitid[u]] = clean!.(models)
            uhp[unitid[u]] = hp
        end
        jldsave(joinpath(resultdir,"eppr.jld2");umodel,umodels,uhp,delay=d,
        siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])",maxcolor,mincolor)
    end
end

function process_condtest_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

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
    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    # jldsave(joinpath(resultdir,"spike.jld2");spike,siteid)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff-condon)
    if haskey(ex,"Eye")
        eye = ex["Eye"]
    else
        eye = ex["Log"]
    end

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    blank = :SpatialFreq=>0
    if ex["ID"] == "Color"
        # factors = [:HueAngle]
        blank = :Color=>36
        factors = [:Angle]
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

    @goto here

    for ii in dataset["imecindex"]
        # Prepare AP
        d = "ap$ii"
        meta = dataset[d]["meta"]
        binfile = meta["fileName"]
        nsavedch = meta["nSavedChans"]
        nsample = meta["nFileSamp"]
        fs = meta["fs"]
        nch = nsavedch-1 # exclude sync channel
        hx,hy,hz = meta["probespacing"]
        exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
        pnrow,pncol = size(exchmask)
        depths = hy*(0:(pnrow-1))
        mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)

        # Depth AP RMS
        epochdur = timetounit(150)
        epoch = [0 epochdur]
        epochs = ref2sync(ccondon.+epoch,dataset,ii)
        # All AP epochs, gain corrected(voltage), bandpass filtered,
        # in the same shape of probe where excluded channels are replaced with local average
        ys = reshape2mask(epochsamplenp(mmbinfile,fs,epochs,1:nch,meta=meta,bandpass=[300,3000]),exchmask)
        # 3 fold downsampling(10kHz)
        ys = ys[:,:,1:3:end,:]
        fs /= 3

        # ΔRMS of combined columns
        @views pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        basedur = timetounit(15)
        baseindex = epoch2sampleindex([0 basedur],fs)
        @views cdrms = Dict(condstring(r)=>
            stfilter(dropdims(mean(mapwindow(rms,pys[:,:,[r.i;r.i.+nrow(ctc)]],(1,101,1),border="symmetric"),dims=3),dims=3),temporaltype=:sub,ti=baseindex)
            for r in eachrow(ccond)) # 10ms rms window
        time = (1/fs/SecondPerUnit)*(0:(size(pys,2)-1))
        if plot
            ylim = mapreduce(y->maximum(abs.(y)),max,values(cdrms))
            for k in keys(cdrms)
                plotanalog(cdrms[k];fs,hy,clims=(-ylim,ylim))
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dMeanRMS$ext")),figfmt)
            end
        end
        jldsave(joinpath(resultdir,"ap$(ii).jld2");cdrms,time,depths,fs,siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")

        # Depth Unit PSTH
        epochdur = timetounit(150)
        epoch = [0 epochdur]
        epochs = ref2sync(ccondon.+epoch,dataset,ii)
        bw = timetounit(2)
        psthbins = epoch[1]:bw:epoch[2]
        basedur = timetounit(15)
        baseindex = epoch2sampleindex([0 basedur],1/(bw*SecondPerUnit))
        # All Unit
        ui = unitsync.==ii
        @views unitepochpsth = map(st->psthspiketrains(epochspiketrain(st,epochs,isminzero=true).y,psthbins,israte=true,ismean=false),unitspike[ui])
        @views cdpsth = Dict(condstring(r)=>begin
                            p = spacepsth(map(p->(;vmeanse(p.mat[r.i,:]).m,p.x),unitepochpsth),unitposition[ui,:],
                                w=replace(unitgood[ui],0=>1.2),lim=(0,hy*(pnrow-1)),bw=2hy,step=hy) # multi-unit count = 1.2 for unit density
                            (;psth=stfilter(mapwindow(mean,p.psth,(1,5),border="symmetric"),temporaltype=:sub,ti=baseindex),p.x,p.y,p.n) # 10ms mean window
                        end  for r in eachrow(ccond))
        if plot
            lim = mapreduce(i->maximum(abs.(i.psth)),max,values(cdpsth))
            for k in keys(cdpsth)
                plotanalog(cdpsth[k].psth;cdpsth[k].x,cdpsth[k].y,clims=(-lim,lim))
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_All-Unit_$(k)_dPSTH$ext")),figfmt)
            end
        end
        jldsave(joinpath(resultdir,"depthpsth$(ii).jld2");cdpsth,siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")

        # Prepare LFP
        d = "lf$ii"
        meta = dataset[d]["meta"]
        binfile = meta["fileName"]
        # lfbin = matchfile(Regex("^$(testid)\\w*.imec.lf.bin"),dir = datadir,join=true)[1]
        nsavedch = meta["nSavedChans"]
        nsample = meta["nFileSamp"]
        fs = meta["fs"]
        nch = nsavedch-1 # exclude sync channel
        hx,hy,hz = meta["probespacing"]
        exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
        pnrow,pncol = size(exchmask)
        depths = hy*(0:(pnrow-1))
        mmbinfile = mmap(binfile,Matrix{Int16},(nsavedch,nsample),0)

        # Depth LFP and CSD
        epochdur = timetounit(150)
        epoch = [0 epochdur]
        epochs = ref2sync(ccondon.+epoch,dataset,ii)
        # All LFP epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered,
        # in the same shape of probe where excluded channels are replaced with local average
        ys = reshape2mask(epochsamplenp(mmbinfile,fs,epochs,1:nch,meta=meta,bandpass=[1,100]),exchmask)

        # LFP and ΔCSD of combined columns
        @views pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        cmlfp = Dict(condstring(r)=>dropdims(mean(pys[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(ccond))
        basedur = timetounit(15)
        baseindex = epoch2sampleindex([0 basedur],fs)
        cdcsd = Dict(k=>stfilter(csd(v,h=hy),temporaltype=:sub,ti=baseindex) for (k,v) in cmlfp)
        time = (1/fs/SecondPerUnit)*(0:(size(pys,2)-1))
        if plot
            llim = mapreduce(y->maximum(abs.(y)),max,values(cmlfp))
            clim = mapreduce(y->maximum(abs.(y)),max,values(cdcsd))
            for k in keys(cmlfp)
                plotanalog(cmlfp[k];fs,hy,clims=(-llim,llim))
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_LFP$ext")),figfmt)

                plotanalog(cdcsd[k];fs,color=:RdBu,hy,clims=(-clim,clim))
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_dCSD$ext")),figfmt)
            end
        end
        jldsave(joinpath(resultdir,"lfp$(ii).jld2");cmlfp,cdcsd,time,depths,fs,siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")

        # Depth Power Spectrum
        # epochdur = timetounit(preicidur)
        # epoch = [0 epochdur]
        # epochs = ref2sync(ccondon.+epoch,dataset,ii)
        # nw = 2
        # ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
        # pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        # ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)
        #
        # epoch = [-epochdur 0]
        # epochs = ref2sync(ccondon.+epoch,dataset,ii)
        # ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
        # pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        # bps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)
        #
        # pcs = ps./bps.-1
        # cmpcs = Dict(condstring(r)=>dropdims(mean(pcs[:,:,[r.i;r.i.+nrow(cctc)]],dims=3),dims=3) for r in eachrow(ccond))
        # if plot
        #     fcmpcs = Dict(k=>imfilter(cmpcs[k],Kernel.gaussian((1,1))) for k in keys(cmpcs))
        #     lim = mapreduce(pc->maximum(abs.(pc)),max,values(fcmpcs))
        #     for k in keys(fcmpcs)
        #         plotanalog(fcmpcs[k],x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],clims=(-lim,lim),color=:vik)
        #         foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_PowerContrast$ext")),figfmt)
        #     end
        # end
        # jldsave(joinpath(resultdir,"lfp$(ii).jld2");cmpcs,depth=depths,freq,siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")
    end

    @label here
    # Unit Spike Train of Condition Epochs
    # if plot
    #     epochext = max(preicidur,suficidur)
    #     epochs = [condon.-epochext condoff.+epochext]
    #     for u in eachindex(unitspike)
    #         ys,ns,ws,is = epochspiketrain(unitspike[u],ref2sync(epochs,dataset,unitsync[u]),isminzero=true,shift=epochext)
    #         title = "$(ugs[u])Unit_$(unitid[u])_SpikeTrain"
    #         plotspiketrain(ys;timeline=[0,conddur],title)
    #         foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
    #     end
    # end

    # Condition Response
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : timetounit(15)
    minresponse = haskey(param,:minresponse) ? param[:minresponse] : 3
    urs=[];uprs=[];usrs=[];ubrs=[];uresponsive=[];umodulative=[]
    cepochs = [ccondon.+responsedelay ccondoff.+responsedelay]
    pcepochs = [ccondon.+(responsedelay-preicidur) ccondon.+responsedelay]
    scepochs = [ccondoff.+responsedelay ccondoff.+(responsedelay+suficidur)]
    bepochs = [bcondon.+responsedelay bcondoff.+responsedelay]
    for u in eachindex(unitspike)
        rs = epochspiketrainresponse_ono(unitspike[u],ref2sync(cepochs,dataset,unitsync[u]),israte=true)
        prs = epochspiketrainresponse_ono(unitspike[u],ref2sync(pcepochs,dataset,unitsync[u]),israte=true)
        srs = epochspiketrainresponse_ono(unitspike[u],ref2sync(scepochs,dataset,unitsync[u]),israte=true)
        push!(urs,rs);push!(uprs,prs);push!(usrs,srs)
        push!(uresponsive,isresponsive(prs,rs,ccond.i))
        push!(umodulative,ismodulative([DataFrame(Y=rs) cctc])) # n-way ANOVA
        if isblank
            br = subrvr(unitspike[u],ref2sync(bepochs,dataset,unitsync[u]))
            mseuc = condresponse(br,condin(bctc))
            push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
        end
        if plot && length(factors)==2
            plotcondresponse(rs,cctc)
            foreach(ext->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_CondResponse$ext")),figfmt)
        end
    end

    # Condition Response in Factor Space
    fis,frs,fms,fses,fa=factorresponse(urs,ccond)
    _,pfrs,pfms,pfses,_=factorresponse(uprs,ccond)
    _,sfrs,sfms,sfses,_=factorresponse(usrs,ccond)

    umaxi=[];uenoughresponse=[];ufrf = Dict(k=>[] for k in keys(fa));umaxf=deepcopy(ufrf);umaxfri=deepcopy(ufrf)
    for u in eachindex(fms)
        oi = Any[Tuple(argmax(replace(fms[u],missing=>-Inf)))...]
        push!(umaxi,oi)
        push!(uenoughresponse,fms[u][oi...]>=minresponse)
        ut = "$(ugs[u][1])U"
        # max response slice for each factor
        for f in keys(fa)
            fdi = findfirst(f.==keys(fa))
            push!(umaxf[f],fa[f][oi[fdi]])
            fdn = length(fa[f])
            fi = deepcopy(oi)
            fi[fdi]=1:fdn
            push!(umaxfri[f],fi)
            rdf=DataFrame(m=fms[u][fi...],se=fses[u][fi...],u=fill(unitid[u],fdn),ug=fill(ut,fdn))
            rdf[:,f]=fa[f]
            prdf=DataFrame(m=pfms[u][fi...],se=pfses[u][fi...],u=fill(unitid[u],fdn),ug=fill("Pre_$ut",fdn))
            prdf[:,f]=fa[f]
            srdf=DataFrame(m=sfms[u][fi...],se=sfses[u][fi...],u=fill(unitid[u],fdn),ug=fill("Suf_$ut",fdn))
            srdf[:,f]=fa[f]
            frf = uenoughresponse[u] ? factorresponsefeature(fa[f],frs[u][fi...],fm=fms[u][fi...],factor=f) : missing
            # frf = factorresponsefeature(fa[f],frs[u][fi...],fm=fms[u][fi...],factor=f)
            push!(ufrf[f],frf)
            if plot
                projection = iscircfactor(f) ? :polar : :none
                plotcondresponse(dropmissing([prdf;srdf;rdf]);color=[:gray85,:gray50,:gray5],linewidth=[1,2,3],grid=true,projection,response=isblank ? ubr[u:u] : [])
                if !ismissing(frf) && !isempty(frf.fit)
                    x = projection == :polar ? (0:0.02:2π) : range(extrema(fa[f])...,step=0.01*mean(diff(sort(fa[f]))))
                    plot!(x->frf.fit.mfit.fun(x,frf.fit.mfit.param),x,lw=2,color=:deepskyblue,label="$(frf.fit.mfit.model)_fit")
                end
                foreach(ext->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_$(f)_Tuning$ext")),figfmt)
            end
        end
    end

    # Simple and Complex Cell
    f1f0=[];tf = getparam(envparam,"TemporalFreq")
    if getparam(envparam,"GratingType") == "Sinusoidal" && tf > 0
        epoch = [0 minconddur]
        bw = 10;fs = 1/(bw/1000)
        binedges = epoch[1]:timetounit(bw):epoch[2]
        for u in eachindex(fms)
            # epochs of max mean response condition
            epochs = ccondon[fis[u][umaxi[u]...]] .+ epoch .+ responsedelay
            m = psthspiketrains(epochspiketrain(unitspike[u],ref2sync(epochs,dataset,unitsync[u]),isminzero=true).y,binedges,israte=true,ismean=true).m
            pm = pfms[u][umaxi[u]...]
            Fs = dft(m.-pm,fs,tf,0)
            push!(f1f0,mapreduce(abs,/,Fs))
        end
    end
    jldsave(joinpath(resultdir,"factorresponse.jld2");fis,frs,fms,fses,pfrs,pfms,pfses,sfrs,sfms,sfses,fa,
    responsive=uresponsive,modulative=umodulative,enoughresponse=uenoughresponse,maxi=umaxi,maxf=umaxf,maxfri=umaxfri,exid=ex["ID"],
    frf=ufrf,f1f0,siteid,unitid,unitgood,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")

    # Single Unit Binary Spike Trian of Condition Tests
    # bepochext = timetounit(-300)
    # bepoch = [-bepochext minconddur]
    # bepochs = ccondon.+bepoch
    # spikebins = bepoch[1]:timetounit(1):bepoch[2]
    # bst = map((st,si)->permutedims(psthspiketrains(epochspiketrain(st,ref2sync(bepochs,dataset,si),isminzero=true,shift=bepochext).y,spikebins,israte=false,ismean=false).mat),unitspike[unitgood],unitsync[unitgood]) # nSpikeBin x nEpoch for each unit
    # # Single Unit Correlogram and Circuit
    # lag=50
    # ccgs,x,ccgis,projs,eunits,iunits,projweights = circuitestimate(bst,lag=lag,unitid=unitid[unitgood],condis=ccond.i)
    #
    # if !isempty(projs)
    #     if plot
    #         for i in eachindex(ccgs)
    #             title = "Correlogram between single unit $(ccgis[i][1]) and $(ccgis[i][2])"
    #             bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
    #             foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
    #         end
    #         plotcircuit(unitposition,unitid,projs,unitgood=unitgood,eunits=eunits,iunits=iunits,layer=layer)
    #         foreach(ext->savefig(joinpath(resultdir,"UnitPosition_Circuit$ext")),figfmt)
    #     end
    #     jldsave(joinpath(resultdir,"circuit.jld2");projs,eunits,iunits,projweights,siteid)
    # end
end

function process_cycle_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

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
    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    # jldsave(joinpath(resultdir,"spike.jld2");spike,siteid)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff-condon)
    if haskey(ex,"Eye")
        eye = ex["Eye"]
    else
        eye = ex["Log"]
    end

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    blank = :SpatialFreq=>0
    if ex["ID"] == "Color"
        # factors = [:HueAngle]
        factors = [:Angle]
        blank = :Color=>36
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

    for ii in dataset["imecindex"]
        # Prepare LFP
        d = "lf$ii"
        lfmeta = dataset[d]["meta"]
        lfbin = lfmeta["fileName"]
        nsavedch = lfmeta["nSavedChans"]
        nsample = lfmeta["nFileSamp"]
        fs = lfmeta["fs"]
        nch = lfmeta["snsApLfSy"][2]
        hx,hy,hz = lfmeta["probespacing"]
        exchmask = exchmasknp(dataset,imecindex=ii,datatype=d)
        pnrow,pncol = size(exchmask)
        depths = hy*(0:(pnrow-1))
        mmlf = Mmap.mmap(lfbin,Matrix{Int16},(nsavedch,nsample),0)

        # Depth Power Spectrum
        epochdur = timetounit(preicidur)
        epoch = [0 epochdur]
        epochs = ref2sync(ccondon.+epoch,dataset,ii)
        nw = 2
        ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
        pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)

        epoch = [-epochdur 0]
        epochs = ref2sync(ccondon.+epoch,dataset,ii)
        ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
        pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
        bps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)

        pcs = ps./bps.-1
        cmpcs = Dict(condstring(r)=>dropdims(mean(pcs[:,:,[r.i;r.i.+nrow(cctc)]],dims=3),dims=3) for r in eachrow(ccond))
        if plot
            fcmpcs = Dict(k=>imfilter(cmpcs[k],Kernel.gaussian((1,1))) for k in keys(cmpcs))
            lim = mapreduce(pc->maximum(abs.(pc)),max,values(fcmpcs))
            for k in keys(fcmpcs)
                plotanalog(fcmpcs[k],x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],clims=(-lim,lim),color=:vik)
                foreach(ext->savefig(joinpath(resultdir,"IMEC$(ii)_$(k)_PowerContrast$ext")),figfmt)
            end
        end
        jldsave(joinpath(resultdir,"lfp$(ii).jld2");cmpcs,depth=depths,freq,siteid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")
    end

    # Unit Position
    if plot
        plotunitposition(unitposition;unitgood,chposition=spike["chposition"],unitid,layer)
        foreach(ext->savefig(joinpath(resultdir,"UnitPosition$ext")),figfmt)
    end

    # Unit Spike Trian of Condition Epochs
    if plot
        epochext = max(preicidur,suficidur)
        epochs = [condon.-epochext condoff.+epochext]
        for u in eachindex(unitspike)
            ys,ns,ws,is = epochspiketrain(unitspike[u],ref2sync(epochs,dataset,unitsync[u]),isminzero=true,shift=epochext)
            title = "$(ugs[u])Unit_$(unitid[u])_SpikeTrian"
            plotspiketrain(ys,timeline=[0,conddur],title=title)
            foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
        end
    end

    # Condition Response
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : timetounit(15)
    minresponse = haskey(param,:minresponse) ? param[:minresponse] : 5
    ubr=[];uresponsive=[];umodulative=[]
    cepochs = [ccondon.+responsedelay ccondoff.+responsedelay]
    pcepochs = [ccondon.+(responsedelay-preicidur) ccondon.+responsedelay]
    scepochs = [ccondoff.+responsedelay ccondoff.+(responsedelay+suficidur)]
    bepochs = [bcondon.+responsedelay bcondoff.+responsedelay]
    for u in eachindex(unitspike)
        rs = epochspiketrainresponse_ono(unitspike[u],ref2sync(cepochs,dataset,unitsync[u]),israte=true)
        prs = epochspiketrainresponse_ono(unitspike[u],ref2sync(pcepochs,dataset,unitsync[u]),israte=true)
        push!(uresponsive,isresponsive(prs,rs,ccond.i))
        push!(umodulative,ismodulative([DataFrame(Y=rs) cctc]))
        if isblank
            br = subrvr(unitspike[u],ref2sync(bepochs,dataset,unitsync[u]))
            mseuc = condresponse(br,condin(bctc))
            push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
        end
        if plot && length(factors)==2
            plotcondresponse(rs,cctc)
            foreach(ext->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_CondResponse$ext")),figfmt)
        end
    end

    # Condition Response in Factor Space
    fms,fses,fa=factorresponse(unitspike,ccond,cepochs,dataset,unitsync)
    pfms,pfses,fa=factorresponse(unitspike,ccond,pcepochs,dataset,unitsync)
    sfms,sfses,fa=factorresponse(unitspike,ccond,scepochs,dataset,unitsync)

    uenoughresponse=[];ufrf = Dict(k=>[] for k in keys(fa));uoptf=deepcopy(ufrf);uoptfri=deepcopy(ufrf)
    for u in eachindex(fms)
        oi = Any[Tuple(argmax(replace(fms[u],missing=>-Inf)))...]
        push!(uenoughresponse,fms[u][oi...]>=minresponse)
        ut = "$(ugs[u][1])U"
        for f in keys(fa)
            fdd = findfirst(f.==keys(fa))
            push!(uoptf[f],fa[f][oi[fdd]])
            fdn = length(fa[f])
            fi = deepcopy(oi)
            fi[fdd]=1:fdn
            push!(uoptfri[f],fi)
            rdf=DataFrame(m=fms[u][fi...],se=fses[u][fi...],u=fill(unitid[u],fdn),ug=fill(ut,fdn))
            rdf[:,f]=fa[f]
            prdf=DataFrame(m=pfms[u][fi...],se=pfses[u][fi...],u=fill(unitid[u],fdn),ug=fill("Pre_$ut",fdn))
            prdf[:,f]=fa[f]
            srdf=DataFrame(m=sfms[u][fi...],se=sfses[u][fi...],u=fill(unitid[u],fdn),ug=fill("Suf_$ut",fdn))
            srdf[:,f]=fa[f]
            push!(ufrf[f],umodulative[u] ? factorresponsefeature(rdf[!,f],rdf[!,:m],factor=f) : missing)
            if plot
                proj = f in [:Ori,:Ori_Final,:HueAngle,:Angle] ? :polar : :cartesian
                df = [rdf;prdf;srdf]
                plotcondresponse(dropmissing(df),color=[:black,:gray70,:gray35],linewidth=[3,1,3],grid=true,projection=proj,response=isblank ? ubr[u:u] : [])
                foreach(ext->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_$(f)_Tuning$ext")),figfmt)
            end
        end
    end

    # Simple and Complex Cell
    f1f0=[];tf = getparam(envparam,"TemporalFreq")
    if getparam(envparam,"GratingType") == "Sinusoidal" && tf > 0
        epoch = [0 minconddur]
        bw = 15;fs = 1/(bw/1000)
        binedges = epoch[1]:timetounit(bw):epoch[2]
        for u in eachindex(fms)
            optcond = findcond(ccond;(;(k=>uoptf[k][u] for k in keys(uoptf))...)...)
            epochs = ccondon[optcond.i[1]] .+ epoch .+ responsedelay
            pepochs = [epochs[:,1].-preicidur epochs[:,1]]
            m = psthspiketrains(epochspiketrain(unitspike[u],ref2sync(epochs,dataset,unitsync[u]),isminzero=true).y,binedges,israte=true,ismean=true).m
            bm = mean(epochspiketrainresponse_ono(unitspike[u],ref2sync(pepochs,dataset,unitsync[u]),israte=true))
            Fs = dft(m.-bm,fs,tf,0)
            push!(f1f0,mapreduce(abs,/,Fs))

            # ps,freq = powerspectrum(m.-bm,fs,freqrange=[0,2tf],nw=2)
            # f0 = ps[1]
            # f1 = ps[argmin(abs.(freq[2:end].-tf))+1]
            # push!(f1f0,sqrt(f1/f0))
        end
    end
    jldsave(joinpath(resultdir,"factorresponse.jld2");fms,fses,pfms,pfses,sfms,sfses,fa,
    responsive=uresponsive,modulative=umodulative,enoughresponse=uenoughresponse,optf=uoptf,optfri=uoptfri,factorresponsefeature=ufrf,
    f1f0,siteid,unitid,eye,color="$(exparam["ColorSpace"])_$(exparam["Color"])")

    # Single Unit Binary Spike Trian of Condition Tests
    bepochext = timetounit(-300)
    bepoch = [-bepochext minconddur]
    bepochs = ccondon.+bepoch
    spikebins = bepoch[1]:timetounit(1):bepoch[2]
    bst = map((st,si)->permutedims(psthspiketrains(epochspiketrain(st,ref2sync(bepochs,dataset,si),isminzero=true,shift=bepochext).y,spikebins,israte=false,ismean=false).mat),unitspike[unitgood],unitsync[unitgood]) # nSpikeBin x nEpoch for each unit
    # Single Unit Correlogram and Circuit
    lag=50
    ccgs,x,ccgis,projs,eunits,iunits,projweights = circuitestimate(bst,lag=lag,unitid=unitid[unitgood],condis=ccond.i)

    if !isempty(projs)
        if plot
            for i in eachindex(ccgs)
                title = "Correlogram between single unit $(ccgis[i][1]) and $(ccgis[i][2])"
                bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
                foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
            end
            plotcircuit(unitposition,unitid,projs,unitgood=unitgood,eunits=eunits,iunits=iunits,layer=layer)
            foreach(ext->savefig(joinpath(resultdir,"UnitPosition_Circuit$ext")),figfmt)
        end
        jldsave(joinpath(resultdir,"circuit.jld2");projs,eunits,iunits,projweights,siteid)
    end
end
