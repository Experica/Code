using Statistics,StatsBase,StatsPlots,Mmap,Images,ePPR

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
    lfbin = matchfile(Regex("^$(testid)\\w*.imec.lf.bin"),dir = datadir,join=true)[1]
    lfmeta = dataset["lf"]["meta"]
    nsavedch = lfmeta["nSavedChans"]
    nsample = lfmeta["nFileSamp"]
    fs = lfmeta["fs"]
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
    # All LFP epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered, in the same shape of probe where excluded channels are replaced with local average
    ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)

    if plot
        for c in 1:pncol
            cmcys = Dict(condstring(r)=>dropdims(mean(ys[:,c,:,r.i],dims=3),dims=3) for r in eachrow(cond))
            cmccsd = Dict(condstring(r)=>dropdims(mean(csd(ys[:,c,:,r.i],h=hy),dims=3),dims=3) for r in eachrow(cond))
            for k in keys(cmcys)
                plotanalog(cmcys[k];fs,cunit=:uv,hy)
                foreach(ext->savefig(joinpath(resultdir,"Column_$(c)_$(k)_MeanLFP$ext")),figfmt)

                plotanalog(imfilter(cmccsd[k],Kernel.gaussian((1,1)));fs,hy)
                foreach(ext->savefig(joinpath(resultdir,"Column_$(c)_$(k)_MeanCSD$ext")),figfmt)
            end
        end
    end

    # Mean LFP and ΔCSD of combined columns
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    cmys = Dict(condstring(r)=>dropdims(mean(pys[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
    pcsd = csd(pys,h=hy)
    basedur = timetounit(15)
    baseindex = epoch2sampleindex([0 basedur],fs)
    dcsd = stfilter(pcsd,temporaltype=:sub,ti=baseindex)
    cmdcsd = Dict(condstring(r)=>dropdims(mean(dcsd[:,:,[r.i;r.i.+nrow(ctc)]],dims=3),dims=3) for r in eachrow(cond))
    x = (1/fs/SecondPerUnit)*(0:(size(dcsd,2)-1))
    if plot
        ylim = mapreduce(y->maximum(abs.(y)),max,values(cmys))
        fcmdcsd = Dict(k=>imfilter(cmdcsd[k],Kernel.gaussian((1,1))) for k in keys(cmdcsd))
        clim = mapreduce(y->maximum(abs.(y)),max,values(fcmdcsd))
        for k in keys(cmys)
            plotanalog(cmys[k];fs,hy,clims=(-ylim,ylim))
            foreach(ext->savefig(joinpath(resultdir,"$(k)_MeanLFP$ext")),figfmt)

            plotanalog(fcmdcsd[k];fs,y=depths,color=:RdBu,hy,clims=(-clim,clim))
            foreach(ext->savefig(joinpath(resultdir,"$(k)_MeandCSD$ext")),figfmt)
        end
    end

    # Depth Power Spectrum
    epochdur = timetounit(300)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    nw = 2
    ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)

    epoch = [-epochdur 0]
    epochs = condon.+epoch
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
            foreach(ext->savefig(joinpath(resultdir,"$(k)_PowerContrast$ext")),figfmt)
        end
    end
    save(joinpath(resultdir,"lfp.jld2"),"cmlfp",cmys,"cmcsd",cmdcsd,"cmpc",cmpcs,"freq",freq,"time",x,"depth",depths,"fs",fs,
    "siteid",siteid,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])")

    # Unit Position
    if plot
        plotunitposition(unitposition;unitgood,chposition=spike["chposition"],unitid,layer)
        foreach(ext->savefig(joinpath(resultdir,"UnitPosition$ext")),figfmt)
    end
    save(joinpath(resultdir,"spike.jld2"),"spike",spike,"siteid",siteid)

    # Unit Spike Trian of Condition Epochs
    if plot
        epochext = max(preicidur,suficidur)
        for u in eachindex(unitspike)
            ys,ns,ws,is = epochspiketrain(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            title = "$(ugs[u])Unit_$(unitid[u])_SpikeTrian"
            plotspiketrain(ys,timeline=[0,conddur],title=title)
            foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
        end
    end

    # Single Unit Depth PSTH
    epochdur = timetounit(150)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    bw = timetounit(5)
    psthbins = epoch[1]:bw:epoch[2]

    baseindex = epoch2sampleindex([0 basedur],1/(bw*SecondPerUnit))
    normfun = x->x.-mean(x[baseindex])
    unitepochpsth = map(st->psthspiketrains(epochspiketrain(st,epochs,isminzero=true).y,psthbins,israte=true,ismean=false),unitspike[unitgood])
    cmdepthpsth = Dict(condstring(r)=>spacepsth(map(pm->(vmeanse(pm.mat[r.i,:],normfun=normfun)...,pm.x),unitepochpsth),unitposition[unitgood,:],hy*(0:pnrow).-(hy/2)) for r in eachrow(cond))
    if plot
        fcdp = Dict(k=>imfilter(cmdepthpsth[k].psth,Kernel.gaussian((1,1))) for k in keys(cmdepthpsth))
        lim = mapreduce(pc->maximum(abs.(pc)),max,values(fcdp))
        for k in keys(cmdepthpsth)
            plotanalog(fcdp[k];cmdepthpsth[k].x,cmdepthpsth[k].y,cmdepthpsth[k].n,timeline=[],clims=(-lim,lim))
            foreach(ext->savefig(joinpath(resultdir,"$(k)_DepthPSTH$ext")),figfmt)
        end
    end
    save(joinpath(resultdir,"depthpsth.jld2"),"cmdepthpsth",cmdepthpsth,"siteid",siteid,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])")

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
        save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits,"projweights",projweights,"siteid",siteid)
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

    # Unit Position
    if plot
        plotunitposition(unitposition;unitgood,chposition=spike["chposition"],unitid,layer)
        foreach(ext->savefig(joinpath(resultdir,"UnitPosition$ext")),figfmt)
    end
    save(joinpath(resultdir,"spike.jld2"),"spike",spike,"siteid",siteid)

    # Prepare Imageset
    bgcolor = RGBA(getparam(envparam,"BGColor")...)
    maxcolor = RGBA(getparam(envparam,"MaxColor")...)
    mincolor = RGBA(getparam(envparam,"MinColor")...)
    masktype = getparam(envparam,"MaskType")
    maskradius = getparam(envparam,"MaskRadius")
    masksigma = getparam(envparam,"Sigma")
    diameter = 5#getparam(envparam,"Diameter")
    ppd = haskey(param,:ppd) ? param[:ppd] : 45
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
        delays = -40:5:200

        if isblank
            uci = unique(ccondidx)
            ucii = map(i->findall(condidx.==i),uci)
            buci = unique(bcondidx)
            bucii = mapreduce(i->findall(condidx.==i),append!,buci)
            bx = vec(mean(mapreduce(i->gray.(imagestimuli[:stimuli][i][xi]),hcat,buci),dims=2))
            x = Array{Float64}(undef,length(uci),xin)
            foreach(i->x[i,:]=gray.(imagestimuli[:stimuli][uci[i]][xi]).-bx,1:size(x,1))

            for u in eachindex(unitspike)
                unitgood[u] || continue
                ys = Array{Float64}(undef,length(delays),length(uci))
                ŷs = Array{Float64}(undef,length(delays),length(uci))
                stas = Array{Float64}(undef,length(delays),xin)
                for d in eachindex(delays)
                    y = epochspiketrainresponse_ono(unitspike[u],condon.+delays[d],condoff.+delays[d],israte=true,isnan2zero=true)
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

            for u in eachindex(unitspike)
                unitgood[u] || continue
                ys = Array{Float64}(undef,length(delays),length(uci))
                ŷs = Array{Float64}(undef,length(delays),length(uci))
                stas = Array{Float64}(undef,length(delays),xin)
                for d in eachindex(delays)
                    y = epochspiketrainresponse_ono(unitspike[u],condon.+delays[d],condoff.+delays[d],israte=true,isnan2zero=true)
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

        save(joinpath(resultdir,"sta.jld2"),"sizepx",sizepx,"x",x,"xi",xi,"xcond",condtable[uci,:],"uy",uy,"uŷ",uŷ,"usta",usta,"ugof",ugof,"delays",delays,
        "siteid",siteid,"sizedeg",sizedeg,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])","maxcolor",maxcolor,"mincolor",mincolor)
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
        save(joinpath(resultdir,"eppr.jld2"),"umodel",umodel,"umodels",umodels,"uhp",uhp,"delay",d,
        "siteid",siteid,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])","maxcolor",maxcolor,"mincolor",mincolor)
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

    # Prepare LFP
    lfbin = matchfile(Regex("^$(testid)\\w*.imec.lf.bin"),dir = datadir,join=true)[1]
    lfmeta = dataset["lf"]["meta"]
    nsavedch = lfmeta["nSavedChans"]
    nsample = lfmeta["nFileSamp"]
    fs = lfmeta["fs"]
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
    ys = reshape2mask(epochsamplenp(mmlf,fs,epochs,1:nch,meta=lfmeta,bandpass=[1,100]),exchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=nw)

    epoch = [-epochdur 0]
    epochs = ccondon.+epoch
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
            foreach(ext->savefig(joinpath(resultdir,"$(k)_PowerContrast$ext")),figfmt)
        end
    end
    save(joinpath(resultdir,"lfp.jld2"),"cmpcs",cmpcs,"depth",depths,"freq",freq,"siteid",siteid,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])")

    # Unit Position
    if plot
        plotunitposition(unitposition;unitgood,chposition=spike["chposition"],unitid,layer)
        foreach(ext->savefig(joinpath(resultdir,"UnitPosition$ext")),figfmt)
    end
    save(joinpath(resultdir,"spike.jld2"),"spike",spike,"siteid",siteid)

    # Unit Spike Trian of Condition Epochs
    if plot
        epochext = max(preicidur,suficidur)
        for u in eachindex(unitspike)
            ys,ns,ws,is = epochspiketrain(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            title = "$(ugs[u])Unit_$(unitid[u])_SpikeTrian"
            plotspiketrain(ys,timeline=[0,conddur],title=title)
            foreach(ext->savefig(joinpath(resultdir,"$title$ext")),figfmt)
        end
    end

    # Single Unit Condition Response
    sui=findall(unitgood)
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : timetounit(15)
    minresponse = haskey(param,:minresponse) ? param[:minresponse] : 5
    ubr=[];uresponsive=[];umodulative=[]
    for i in sui
        rs = epochspiketrainresponse_ono(unitspike[i],ccondon.+responsedelay,ccondoff.+responsedelay,israte=true)
        prs = epochspiketrainresponse_ono(unitspike[i],ccondon.+(responsedelay-preicidur),ccondon.+responsedelay,israte=true)
        push!(uresponsive,isresponsive(prs,rs,ccond.i))
        push!(umodulative,ismodulative([DataFrame(Y=rs) cctc]))
        if isblank
            br = subrvr(unitspike[i],bcondon.+responsedelay,bcondoff.+responsedelay)
            mseuc = condresponse(br,condin(bctc))
            push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
        end
        if plot && length(factors)==2
            plotcondresponse(rs,cctc)
            foreach(ext->savefig(joinpath(resultdir,"Single-Unit_$(unitid[i])_CondResponse$ext")),figfmt)
        end
    end

    # Condition Response in Factor Space
    fms,fses,fa=factorresponse(unitspike[sui],ccond,ccondon.+responsedelay,ccondoff.+responsedelay)
    pfms,pfses,fa=factorresponse(unitspike[sui],ccond,ccondon.+(responsedelay-preicidur),ccondon.+responsedelay)
    sfms,sfses,fa=factorresponse(unitspike[sui],ccond,ccondoff.+responsedelay,ccondoff.+(responsedelay+suficidur))

    uenoughresponse=[];ufrf = Dict(k=>[] for k in keys(fa));uoptf=deepcopy(ufrf);uoptfri=deepcopy(ufrf)
    for i in eachindex(fms)
        oi = Any[Tuple(argmax(replace(fms[i],missing=>-Inf)))...]
        push!(uenoughresponse,fms[i][oi...]>=minresponse)
        for f in keys(fa)
            fd = findfirst(f.==keys(fa))
            push!(uoptf[f],fa[f][oi[fd]])
            fdn = length(fa[f])
            fi = deepcopy(oi)
            fi[fd]=1:fdn
            push!(uoptfri[f],fi)
            rdf=DataFrame(m=fms[i][fi...],se=fses[i][fi...],u=fill(unitid[sui[i]],fdn),ug=fill("SU",fdn))
            rdf[:,f]=fa[f]
            prdf=DataFrame(m=pfms[i][fi...],se=pfses[i][fi...],u=fill(unitid[sui[i]],fdn),ug=fill("Pre_SU",fdn))
            prdf[:,f]=fa[f]
            srdf=DataFrame(m=sfms[i][fi...],se=sfses[i][fi...],u=fill(unitid[sui[i]],fdn),ug=fill("Suf_SU",fdn))
            srdf[:,f]=fa[f]
            push!(ufrf[f],umodulative[i] ? factorresponsefeature(rdf[!,f],rdf[!,:m],factor=f) : missing)
            if plot
                proj = f in [:Ori,:Ori_Final,:HueAngle,:Angle] ? :polar : :cartesian
                df = [rdf;prdf;srdf]
                plotcondresponse(dropmissing(df),color=[:black,:gray70,:gray35],linewidth=[3,1,3],grid=true,projection=proj,response=isblank ? ubr[i:i] : [])
                foreach(ext->savefig(joinpath(resultdir,"Single-Unit_$(unitid[sui[i]])_$(f)_Tuning$ext")),figfmt)
            end
        end
    end

    # Simple and Complex Cell
    f1f0=[];tf = getparam(envparam,"TemporalFreq")
    if getparam(envparam,"GratingType") == "Sinusoidal" && tf > 0 && getparam(envparam,"Drifting")
        epoch = [0 minconddur]
        bw = 15;fs = 1/(bw/1000)
        binedges = epoch[1]:timetounit(bw):epoch[2]
        for i in eachindex(fms)
            optcond = findcond(ccond;(;(k=>uoptf[k][i] for k in keys(uoptf))...)...)
            epochs = ccondon[optcond.i[1]] .+ epoch .+ responsedelay
            m = psthspiketrains(epochspiketrain(unitspike[sui[i]],epochs,isminzero=true).y,binedges,israte=true,ismean=true).m
            bm = mean(epochspiketrainresponse_ono(unitspike[sui[i]],epochs[:,1].-preicidur,epochs[:,1],israte=true))
            ps,freq = powerspectrum(m.-bm,fs,freqrange=[0,2tf],nw=2)
            f0 = ps[1]
            f1 = ps[argmin(abs.(freq[2:end].-tf))+1]
            push!(f1f0,sqrt(f1/f0))
        end
    end
    save(joinpath(resultdir,"factorresponse.jld2"),"fms",fms,"fses",fses,"pfms",pfms,"pfses",pfses,"sfms",sfms,"sfses",sfses,"fa",fa,
    "responsive",uresponsive,"modulative",umodulative,"enoughresponse",uenoughresponse,"optf",uoptf,"optfri",uoptfri,"factorresponsefeature",ufrf,
    "f1f0",f1f0,"siteid",siteid,"unitid",unitid[sui],"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])")

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
        save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits,"projweights",projweights,"siteid",siteid)
    end
end
