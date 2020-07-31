using Statistics,StatsPlots,Mmap,Images,StatsBase

function process_flash_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    datadir = joinpath(param[:dataroot],subject,siteid)
    resultsitedir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(resultsitedir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood=spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff-condon)

    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)

    # Prepare LFP
    lfbin = matchfile(Regex("^$(testid)[A-Za-z0-9_]*.imec.lf.bin"),dir = datadir,join=true)[1]
    lfmeta = dataset["lf"]["meta"]
    nsavedch=lfmeta["nSavedChans"]
    nsample=lfmeta["nFileSamp"]
    fs=lfmeta["fs"]
    nch = lfmeta["snsApLfSy"][2]
    hx,hy,hz = lfmeta["probespacing"]
    exchmask = exchmasknp(dataset,type="lf")
    pnrow,pncol = size(exchmask)
    depths = hy*(0:(pnrow-1))
    mmlf = Mmap.mmap(lfbin,Matrix{Int16},(nsavedch,nsample),0)

    # Depth LFP and CSD
    epochdur = timetounit(150)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    # All LFP epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered, in the same shape of probe, where excluded channels are replaced with local average
    ys=reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)

    if plot
        for c in 1:pncol
            cmcys = Dict(condstring(r)=>dropdims(mean(ys[:,c,:,r.i],dims=3),dims=3) for r in eachrow(cond))
            cmccsd = Dict(condstring(r)=>dropdims(mean(csd(ys[:,c,:,r.i],h=hy),dims=3),dims=3) for r in eachrow(cond))
            for k in keys(cmcys)
                plotanalog(cmcys[k],fs=fs,cunit=:uv)
                foreach(ext->savefig(joinpath(resultdir,"Column_$(c)_$(k)_MeanLFP$ext")),figfmt)

                plotanalog(imfilter(cmccsd[k],Kernel.gaussian((1,1))),fs=fs)
                foreach(ext->savefig(joinpath(resultdir,"Column_$(c)_$(k)_MeanCSD$ext")),figfmt)
            end
        end
    end

    # Mean LFP and ΔCSD of combined columns
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    mys=dropdims(mean(pys,dims=3),dims=3)
    cmys = Dict(condstring(r)=>dropdims(mean(pys[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
    pcsd = csd(pys,h=hy)
    basedur = timetounit(15)
    baseindex = epoch2sampleindex([0 basedur],fs)
    dcsd=stfilter(pcsd,temporaltype=:sub,ti=baseindex)
    mdcsd = dropdims(mean(dcsd,dims=3),dims=3)
    cmdcsd = Dict(condstring(r)=>dropdims(mean(dcsd[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
    x = (1/fs/SecondPerUnit)*(0:(size(mdcsd,2)-1))
    if plot
        plotanalog(mys,fs=fs,cunit=:uv)
        foreach(ext->savefig(joinpath(resultdir,"MeanLFP$ext")),figfmt)

        plotanalog(imfilter(mdcsd,Kernel.gaussian((1,1))),fs=fs,y=depths,color=:RdBu)
        foreach(ext->savefig(joinpath(resultdir,"MeandCSD$ext")),figfmt)

        for k in keys(cmys)
            plotanalog(cmys[k],fs=fs,cunit=:uv)
            foreach(ext->savefig(joinpath(resultdir,"$(k)_MeanLFP$ext")),figfmt)

            plotanalog(imfilter(cmdcsd[k],Kernel.gaussian((1,1))),fs=fs,y=depths,color=:RdBu)
            foreach(ext->savefig(joinpath(resultdir,"$(k)_MeandCSD$ext")),figfmt)
        end
    end
    save(joinpath(resultdir,"lfp.jld2"),"lfp",mys,"clfp",cmys,"time",x,"depth",depths,"fs",fs,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])")
    save(joinpath(resultdir,"csd.jld2"),"csd",mdcsd,"ccsd",cmdcsd,"time",x,"depth",depths,"fs",fs,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])")

    # Depth Power Spectrum
    epochdur = timetounit(500)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    bw = 4
    ys=reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=bw*epochdur*SecondPerUnit)

    epoch = [-epochdur 0]
    epochs = condon.+epoch
    ys=reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    bps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=bw*epochdur*SecondPerUnit)

    pcs = ps./bps.-1
    mpcs = dropdims(mean(pcs,dims=3),dims=3)
    cmpcs = Dict(condstring(r)=>dropdims(mean(pcs[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
    if plot
        mps = dropdims(mean(ps,dims=3),dims=3)
        plotanalog(imfilter(mps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
        foreach(ext->savefig(joinpath(resultdir,"PowerSpectrum$ext")),figfmt)

        mbps =dropdims(mean(bps,dims=3),dims=3)
        plotanalog(imfilter(mbps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
        foreach(ext->savefig(joinpath(resultdir,"PowerSpectrum_Baseline$ext")),figfmt)

        plotanalog(imfilter(mpcs,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
        foreach(ext->savefig(joinpath(resultdir,"PowerContrast$ext")),figfmt)

        cmps = Dict(condstring(r)=>dropdims(mean(ps[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
        cmbps = Dict(condstring(r)=>dropdims(mean(bps[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
        for k in keys(cmpcs)
            plotanalog(imfilter(cmps[k],Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
            foreach(ext->savefig(joinpath(resultdir,"$(k)_PowerSpectrum$ext")),figfmt)

            plotanalog(imfilter(cmbps[k],Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
            foreach(ext->savefig(joinpath(resultdir,"$(k)_PowerSpectrum_Baseline$ext")),figfmt)

            plotanalog(imfilter(cmpcs[k],Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
            foreach(ext->savefig(joinpath(resultdir,"$(k)_PowerContrast$ext")),figfmt)
        end
    end
    save(joinpath(resultdir,"powercontrast.jld2"),"pcs",mpcs,"cpcs",cmpcs,"depth",depths,"freq",freq)

    # Unit Position
    if plot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(ext->savefig(joinpath(resultdir,"UnitPosition$ext")),figfmt)
    end
    spike["siteid"] = siteid
    save(joinpath(resultdir,"spike.jld2"),"spike",spike)

    # Unit Spike Train
    if plot
        epochext = preicidur
        for u in eachindex(unitspike)
            ys,ns,ws,is = epochspiketrain(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            title = "$(ugs[u])Unit_$(unitid[u])_SpikeTrian"
            plotspiketrain(ys,timeline=[0,conddur],title=title)
            foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
        end
    end

    # Unit Depth PSTH
    epochdur = timetounit(150)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    bw = timetounit(2)
    psthbins = epoch[1]:bw:epoch[2]

    baseindex = epoch2sampleindex([0 basedur],1/(bw*SecondPerUnit))
    normfun = x->x.-mean(x[baseindex])
    unitepochpsth = map(st->psthspiketrains(epochspiketrain(st,epochs,isminzero=true).y,psthbins,israte=true,ismean=false),unitspike[unitgood])
    unitpsth = map(pm->(vmeanse(pm.mat,normfun=normfun)...,pm.x),unitepochpsth)
    depthpsth,x,depths,depthnunit = spacepsth(unitpsth,unitposition[unitgood,:],hy*(0:pnrow).-(hy/2))
    cdepthpsth = Dict(condstring(r)=>spacepsth(map(pm->(vmeanse(pm.mat[[r.i;r.i.+nrow(ctc)],:],normfun=normfun)...,pm.x),unitepochpsth),unitposition[unitgood,:],hy*(0:pnrow).-(hy/2))[1] for r in eachrow(cond))
    if plot
        plotpsth(imfilter(depthpsth,Kernel.gaussian((1,1))),x,depths,n=depthnunit)
        foreach(ext->savefig(joinpath(resultdir,"DepthPSTH$ext")),figfmt)

        for k in keys(cdepthpsth)
            plotpsth(imfilter(cdepthpsth[k],Kernel.gaussian((1,1))),x,depths,n=depthnunit)
            foreach(ext->savefig(joinpath(resultdir,"$(k)_DepthPSTH$ext")),figfmt)
        end
    end
    save(joinpath(resultdir,"depthpsth.jld2"),"depthpsth",depthpsth,"cdepthpsth",cdepthpsth,"depth",depths,"time",x,"n",depthnunit)

    # Single Unit Binary Spike Train of Condition Tests
    bepochext = timetounit(-300)
    bepoch = [-bepochext minconddur]
    bepochs = condon.+bepoch
    spikebins = bepoch[1]:timetounit(1):bepoch[2] # 1ms bin histogram will convert spike times to binary spike train
    bst = map(st->permutedims(psthspiketrains(epochspiketrain(st,bepochs,isminzero=true,shift=bepochext).y,spikebins,israte=false,ismean=false).mat),unitspike[unitgood]) # nSpikeBin x nEpoch for each unit
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
        save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits,"projweights",projweights)
    end
end

function process_hartley_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    datadir = joinpath(param[:dataroot],subject,siteid)
    resultsitedir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(resultsitedir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood=spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]

    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]

    # Prepare Conditions
    condtable = DataFrame(ex["Cond"])
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    blank = haskey(param,:blank) ? param[:blank] : (:SpatialFreq,0)
    ci = ctc[!,blank[1]].!==blank[2]
    cctc = ctc[ci,factors]
    ccondon = condon[ci]
    ccondoff = condoff[ci]
    ccondidx = condidx[ci]
    bi = .!ci
    bctc = ctc[bi,factors]
    bcondon = condon[bi]
    bcondoff = condoff[bi]
    bcondidx = condidx[bi]
    isblank = !isempty(bcondon)

    # Unit Position
    if plot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(ext->savefig(joinpath(resultdir,"UnitPosition$ext")),figfmt)
    end
    spike["siteid"] = siteid
    save(joinpath(resultdir,"spike.jld2"),"spike",spike)


    # Prepare Imageset
    nscale = haskey(param,:nscale) ? param[:nscale] : 2
    downsample = haskey(param,:downsample) ? param[:downsample] : 2
    sigma = haskey(param,:sigma) ? param[:sigma] : 1.5
    ppd = haskey(param,:ppd) ? param[:ppd] : 45
    bgcolor = RGBA(getparam(envparam,"BGColor")...)
    maxcolor = RGBA(getparam(envparam,"MaxColor")...)
    mincolor = RGBA(getparam(envparam,"MinColor")...)
    masktype = getparam(envparam,"MaskType")
    maskradius = getparam(envparam,"MaskRadius")
    masksigma = getparam(envparam,"Sigma")
    diameter = 5#10#getparam(envparam,"Diameter")
    # ii = round(Int,4.5ppd):round(Int,8.5ppd)
    # jj = range(round(Int,5.5ppd),length=length(ii))
    # d=round(length(ii)/ppd,digits=1)
    sizedeg = (diameter,diameter)
    imagesetname = splitext(splitdir(ex["CondPath"])[2])[1] * "_size$sizedeg"
    if !haskey(param,imagesetname)
        imageset = map(i->GrayA.(grating(θ=deg2rad(i.Ori),sf=i.SpatialFreq,phase=i.SpatialPhase,size=sizedeg,ppd=ppd)),eachrow(condtable))
        imageset = Dict{Symbol,Any}(:pyramid => map(i->gaussian_pyramid(i, nscale-1, downsample, sigma),imageset))
        imageset[:sizepx] = map(i->size(i),imageset[:pyramid][1])
        param[imagesetname] = imageset
    end
    # sizedeg=(d,d)
    # Prepare Image Stimuli
    imageset = param[imagesetname]
    bgcolor = oftype(imageset[:pyramid][1][1][1],bgcolor)
    unmaskindex = map(i->alphamask(i,radius=maskradius,sigma=masksigma,masktype=masktype)[2],imageset[:pyramid][1])
    imagestimuli = map(s->map(i->alphablend.(alphamask(i[s],radius=maskradius,sigma=masksigma,masktype=masktype)[1],[bgcolor]),imageset[:pyramid]),1:nscale)




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
        scaleindex=1
        sizepx = imageset[:sizepx][scaleindex]
        xi = unmaskindex[scaleindex]
        uy=Dict();usta = Dict()
        delays = -30:5:210

        if isblank
            uci = unique(ccondidx)
            ucii = map(i->findall(condidx.==i),uci)
            buci = unique(bcondidx)
            bucii = mapreduce(i->findall(condidx.==i),append!,buci)
            bx = mean(mapreduce(i->gray.(imagestimuli[scaleindex][i][xi]),hcat,buci),dims=2)[:]
            x = Array{Float64}(undef,length(uci),length(xi))
            foreach(i->x[i,:]=gray.(imagestimuli[scaleindex][uci[i]][xi]).-bx,1:size(x,1))

            for u in eachindex(unitspike)
                !unitgood[u] && continue
                ys = Array{Float64}(undef,length(delays),length(uci))
                stas = Array{Float64}(undef,length(delays),length(xi))
                for d in eachindex(delays)
                    y = epochspiketrainresponse_ono(unitspike[u],condon.+delays[d],condoff.+delays[d],israte=true,isnan2zero=true)
                    by = mean(y[bucii])
                    y = map(i->mean(y[i])-by,ucii)
                    ys[d,:] = y
                    stas[d,:]=sta(x,y)
                end
                uy[unitid[u]]=ys
                usta[unitid[u]]=stas

                if plot
                    r = [extrema(stas)...]
                    for d in eachindex(delays)
                        title = "$(ugs[u])Unit_$(unitid[u])_STA_$(delays[d])"
                        p = plotsta(stas[d,:],sizepx=sizepx,sizedeg=sizedeg,index=xi,title=title,r=r)
                        foreach(ext->save(joinpath(resultdir,"$title$ext"),p),figfmt)
                    end
                end
            end
        else
            uci = unique(condidx)
            ucii = map(i->findall(condidx.==i),uci)
            x = Array{Float64}(undef,length(uci),length(xi))
            foreach(i->x[i,:]=gray.(imagestimuli[scaleindex][uci[i]][xi]),1:size(x,1))

            for u in eachindex(unitspike)
                !unitgood[u] && continue
                ys = Array{Float64}(undef,length(delays),length(uci))
                stas = Array{Float64}(undef,length(delays),length(xi))
                for d in eachindex(delays)
                    y = epochspiketrainresponse_ono(unitspike[u],condon.+delays[d],condoff.+delays[d],israte=true,isnan2zero=true)
                    y = map(i->mean(y[i]),ucii)
                    ys[d,:] = y
                    stas[d,:]=sta(x,y)
                end
                uy[unitid[u]]=ys
                usta[unitid[u]]=stas

                if plot
                    r = [extrema(stas)...]
                    for d in eachindex(delays)
                        title = "$(ugs[u])Unit_$(unitid[u])_STA_$(delays[d])"
                        p = plotsta(stas[d,:],sizepx=sizepx,sizedeg=sizedeg,index=xi,title=title,r=r)
                        foreach(ext->save(joinpath(resultdir,"$title$ext"),p),figfmt)
                    end
                end
            end
        end
        save(joinpath(resultdir,"sta.jld2"),"sizepx",sizepx,"x",x,"xi",xi,"xcond",condtable[uci,:],"uy",uy,"usta",usta,"delays",delays,
        "sizedeg",sizedeg,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])","maxcolor",maxcolor,"mincolor",mincolor)
    end

    if :ePPR in param[:model]
        for u in 75:75, d in 0.08:0.01:0.08
            cellid = "$(siteid)_U$(unitid[u])"
            unitresultdir = joinpath(resultdir,"Unit_$(unitid[u])_ePPR_$d")
            isdir(unitresultdir) && rm(unitresultdir,recursive=true)
            y = subrvr(unitspike[u],ccondon.+d,ccondoff.+d)

            debug = plot ? ePPRDebugOptions(level=DebugVisual,logdir=unitresultdir) : ePPRDebugOptions()
            hp = ePPRHyperParams(imagesize...,ndelay=param[:epprndelay],blankcolor=128,nft=param[:epprnft],lambda=param[:epprlambda])
            model,models = epprhypercv(x,y,hp,debug)

            if plot && !isnothing(model)
                debug(plotmodel(model,hp),log="Model_Final (λ=$(hp.lambda))")
            end
        end
    end
end

function process_condtest_spikeglx(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    datadir = joinpath(param[:dataroot],subject,siteid)
    resultsitedir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(resultsitedir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # Condition Tests
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff-condon)

    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    blank = (:SpatialFreq,0)
    if ex["ID"] == "Color"
        factors = [:HueAngle]
        blank = (:Color,36)
    end
    haskey(param,:blank) && (blank = param[:blank])
    ci = ctc[!,blank[1]].!==blank[2]
    cctc = ctc[ci,factors]
    ccond = condin(cctc)
    ccondon=condon[ci]
    ccondoff=condoff[ci]
    bi = .!ci
    bctc = ctc[bi,factors]
    bcondon=condon[bi]
    bcondoff=condoff[bi]
    isblank = !isempty(bcondon)

    # Prepare LFP
    lfbin = matchfile(Regex("^$(testid)[A-Za-z0-9_]*.imec.lf.bin"),dir = datadir,join=true)[1]
    lfmeta = dataset["lf"]["meta"]
    nsavedch=lfmeta["nSavedChans"]
    nsample=lfmeta["nFileSamp"]
    fs=lfmeta["fs"]
    nch = lfmeta["snsApLfSy"][2]
    hx,hy,hz = lfmeta["probespacing"]
    exchmask = exchmasknp(dataset,type="lf")
    pnrow,pncol = size(exchmask)
    depths = hy*(0:(pnrow-1))
    mmlf = Mmap.mmap(lfbin,Matrix{Int16},(nsavedch,nsample),0)

    # Depth Power Spectrum
    epochdur = timetounit(preicidur)
    epoch = [0 epochdur]
    epochs = ccondon.+epoch
    nw = 2
    ys=reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)

    epoch = [-epochdur 0]
    epochs = ccondon.+epoch
    ys=reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    bps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)

    pcs = ps./bps.-1
    cmpcs = Dict(condstring(r)=>dropdims(mean(pcs[:,:,[r.i;r.i.+nrow(cctc)]],dims=3),dims=3) for r in eachrow(ccond))
    if plot
        fcmpcs = Dict(k=>imfilter(cmpcs[k],Kernel.gaussian((1,1))) for k in keys(cmpcs))
        lim = mapreduce(pc->maximum(abs.(pc)),max,values(fcmpcs))
        for k in keys(fcmpcs)
            plotanalog(fcmpcs[k],x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],clims=(-lim,lim),color=:vik)
            foreach(ext->savefig(joinpath(resultdir,"$(k)_PowerContrast$ext")),figfmt)
        end
    end
    save(joinpath(resultdir,"powercontrast.jld2"),"cpcs",cmpcs,"depth",depths,"freq",freq)

    # Unit Position
    if plot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(ext->savefig(joinpath(resultdir,"UnitPosition$ext")),figfmt)
    end
    spike["siteid"] = siteid
    save(joinpath(resultdir,"spike.jld2"),"spike",spike)

    # Unit Spike Trian
    if plot
        epochext = preicidur
        for u in eachindex(unitspike)
            ys,ns,ws,is = epochspiketrain(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            title = "$(ugs[u])Unit_$(unitid[u])_SpikeTrian"
            plotspiketrain(ys,timeline=[0,conddur],title=title)
            foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
        end
    end

    # Condition Response
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : timetounit(15)
    minresponse = haskey(param,:minresponse) ? param[:minresponse] : 5
    ubr=[];uresponsive=[];umodulative=[]
    for u in eachindex(unitspike)
        rs = epochspiketrainresponse_ono(unitspike[u],ccondon.+responsedelay,ccondoff.+responsedelay,israte=true)
        prs = epochspiketrainresponse_ono(unitspike[u],ccondon.+(responsedelay-preicidur),ccondon.+responsedelay,israte=true)
        push!(uresponsive,isresponsive(prs,rs,ccond.i))
        push!(umodulative,ismodulative([DataFrame(Y=rs) cctc]))
        if isblank
            br = subrvr(unitspike[u],bcondon.+responsedelay,bcondoff.+responsedelay)
            mseuc = condresponse(br,condin(bctc))
            push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
        end
        if plot && length(factors)==2
            plotcondresponse(rs,cctc)
            foreach(ext->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_CondResponse$ext")),figfmt)
        end
    end

    # Condition Response in Factor Space
    fms,fses,fa=factorresponse(unitspike,ccond,ccondon.+responsedelay,ccondoff.+responsedelay)
    pfms,pfses,fa=factorresponse(unitspike,ccond,ccondon.+(responsedelay-preicidur),ccondon.+responsedelay)
    sfms,sfses,fa=factorresponse(unitspike,ccond,ccondoff.+responsedelay,ccondoff.+(responsedelay+suficidur))

    uenoughresponse=[];ufrf = Dict(k=>[] for k in keys(fa));uoptcond=[]
    for u in eachindex(fms)
        oi = Any[Tuple(argmax(replace(fms[u],missing=>-Inf)))...]
        push!(uenoughresponse,fms[u][oi...]>=minresponse)
        push!(uoptcond, Dict(map((f,i)->f=>fa[f][i],keys(fa),oi)))
        for f in keys(fa)
            fd = findfirst(f.==keys(fa))
            fdn = length(fa[f])
            fi = deepcopy(oi)
            fi[fd]=1:fdn
            rdf=DataFrame(m=fms[u][fi...],se=fses[u][fi...],u=fill(unitid[u],fdn),ug=fill("$(ugs[u][1])U",fdn))
            rdf[:,f]=fa[f]
            prdf=DataFrame(m=pfms[u][fi...],se=pfses[u][fi...],u=fill(unitid[u],fdn),ug=fill("Pre_$(ugs[u][1])U",fdn))
            prdf[:,f]=fa[f]
            srdf=DataFrame(m=sfms[u][fi...],se=sfses[u][fi...],u=fill(unitid[u],fdn),ug=fill("Suf_$(ugs[u][1])U",fdn))
            srdf[:,f]=fa[f]
            push!(ufrf[f],factorresponsefeature(rdf[!,f],rdf[!,:m],factor=f))
            if plot
                proj = f in [:Ori,:Ori_Final,:HueAngle] ? :polar : :cartesian
                df = [rdf;prdf;srdf]
                plotcondresponse(dropmissing(df),color=[:black,:gray70,:gray35],linewidth=[3,1,3],grid=true,projection=proj,response=isblank ? ubr[u:u] : [])
                foreach(ext->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_$(f)_Tuning$ext")),figfmt)
            end
        end
    end
    save(joinpath(resultdir,"factorresponse.jld2"),"optcond",uoptcond,"factorresponsefeature",ufrf,"fms",fms,"fses",fses,"pfms",pfms,"pfses",pfses,"sfms",sfms,"sfses",sfses,"fa",fa,
    "responsive",uresponsive,"modulative",umodulative,"enoughresponse",uenoughresponse,"unitgood",unitgood,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])")

    # Single Unit Binary Spike Trian of Condition Tests
    bepochext = timetounit(-300)
    bepoch = [-bepochext minconddur]
    bepochs = ccondon.+bepoch
    spikebins = bepoch[1]:timetounit(1):bepoch[2]
    bst = map(st->permutedims(psthspiketrains(epochspiketrain(st,bepochs,isminzero=true,shift=bepochext).y,spikebins,israte=false,ismean=false).mat),unitspike[unitgood]) # nSpikeBin x nEpoch for each unit
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
        save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits,"projweights",projweights)
    end
end
