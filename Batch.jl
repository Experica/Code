using NeuroAnalysis,ProgressMeter,Logging,Statistics,DataFrames,Plots,Mmap,Images,StatsBase,ePPR
import Base: close

function batchtests(tests::DataFrame,param::Dict{Any,Any}=Dict{Any,Any}();log::Dict{Any,AbstractLogger}=Dict{Any,AbstractLogger}(),isplot=true)
    rdf=[];cdf=[]
    p = ProgressMeter.Progress(nrow(tests),1,"Batch Tests ... ",50)
    for t in eachrow(tests)
        try
            id = t[:ID];r=nothing
            if id=="OriGrating"
                u,c=processori(dataset,resultroot,uuid=uuid,log=log,delay=delay,binwidth=binwidth,isplot=isplot)
            elseif id=="Laser"
                u,c=processlaser(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
            elseif id=="Image"
                u,c=processimage(dataset,resultroot,env,uuid=uuid,log=log,delay=delay,binwidth=binwidth,isplot=isplot)
            elseif id=="LaserImage"
                u,c=processlaserimage(dataset,condroot,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
            elseif id=="Flash"
                r,c=processflash(t[:files],param,uuid=t[:UUID],log=log,isplot=isplot)
            elseif id=="Hartley"
                u,c=processhartley(t[:files],param,uuid=t[:UUID],log=log,isplot=isplot)
            elseif id in ["OriSF","OriSFColor"]
                r,c=processcondtest(t[:files],param,uuid=t[:UUID],log=log,isplot=isplot)
            end
            if !isnothing(r)
                push!(rdf,r);push!(cdf,c)
            end
        catch exc
            display("============================================")
            display("Error In Processing: $(t[:files])")
            display("============================================")
            display.(stacktrace(catch_backtrace()))
        end
        next!(p)
    end
    close(log)
    return vcat(rdf...),vcat(cdf...)
end

function close(log::Dict{Any,AbstractLogger})
    for l in values(log)
        close(l.stream)
    end
end
function (log::Dict{Any,AbstractLogger})(key,f,logfile)
    if !haskey(log,key)
        log[key] = SimpleLogger(open(logfile,"w"))
    end
    logger = log[key]
    if f isa String
        println(logger.stream,f)
    else
        with_logger(f,logger)
    end
    flush(logger.stream)
end

function processlaserimage(dataset::Dict,condroot::Dict{String,Any},resultroot;uuid="",delay=20,binwidth=10,minpredur=10,mincondtest=12000,
    nscale=2,downsample=2,sigma=1.5,pixelscale=255,isplot=true,forceprocess=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    bgcolor=RGBA(getparam(envparam,"BGColor")...)
    imagesetname = replace(getparam(envparam,"ImageSet","ImageQuad"),"Ã—","_")
    imagemasktype = getparam(envparam,"MaskType","ImageQuad")
    imagemaskradius = getparam(envparam,"MaskRadius","ImageQuad")
    imagemasksigma = getparam(envparam,"Sigma","ImageQuad")
    ct,ctc=ctctc(ex)

    lfactor = filter(f->f!=:Image,names(ctc))
    lcond = condin(ctc[:,lfactor])
    all(lcond[:n] .< mincondtest) && return [],[]

    if !haskey(condroot,imagesetname) && haskey(condroot,"rootdir") && isdir(condroot["rootdir"])
        pyramid = Dict{Symbol,Any}(:pyramid => map(i->gaussian_pyramid(i, nscale-1, downsample, sigma),
        loadimageset(joinpath(condroot["rootdir"],imagesetname),alpha=true)))
        pyramid[:size] = map(i->size(i),pyramid[:pyramid][1])
        condroot[imagesetname] = pyramid
    end
    imageset = condroot[imagesetname]
    bgimagecolor = oftype(imageset[:pyramid][1][1][1],bgcolor)
    unmaskindex = map(i->alphamask(i,radius=imagemaskradius,sigma=imagemasksigma,masktype=imagemasktype)[2],imageset[:pyramid][1])
    imagestimuli = map(s->map(i->alphablend.(alphamask(i[s],radius=imagemaskradius,sigma=imagemasksigma,masktype=imagemasktype)[1],[bgimagecolor]),imageset[:pyramid]),1:nscale)

    s = 2
    ximagesize = imageset[:size][s]
    xi = unmaskindex[s]
    imagestimulimatrix = Array{Float64}(length(imagestimuli[s]),prod(ximagesize))
    for i in 1:size(imagestimulimatrix,1)
        imagestimulimatrix[i,:] = vec(gray.(imagestimuli[s][i]))
    end

    predur = max(preicidur,minpredur)
    resultroot=abspath(resultroot,ex["ID"])
    !isdir(resultroot) && mkpath(resultroot)

    udf = []
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                preurs = subrvr(est[esu.==u],ct[:CondOn]+delay-predur,ct[:CondOn]+delay)
                y = subrvr(est[esu.==u],ct[:CondOn]+delay,ct[:CondOff]+delay)
                !isresponsive(preurs,y,lcond[:i]) && continue

                udir = joinpath(resultroot,"$(uuid)_E$(e)_U$(u)")
                !forceprocess && isdir(udir) && continue
                plotdir = isplot ? udir : nothing
                mrs = []
                for l in eachrow(lcond)
                    debug = isplot ? ePPRDebugOptions(level=DebugVisual,logdir=joinpath(plotdir,condstring(l,lfactor))) : ePPRDebugOptions()
                    hp = ePPRHyperParams(ximagesize...,xindex=xi,ndelay=1,blankcolor=gray(bgimagecolor)*pixelscale)
                    hp.nft = [6]
                    hp.lambda = 30000
                    model,models = epprcv(imagestimulimatrix[Int.(ctc[:Image][l[:i]]),:]*pixelscale,y[l[:i]],hp,debug)

                    if isplot && model!=nothing
                        debug(plotmodel(model,hp),log="Model_Final")
                    end
                    push!(mrs,(clean(model),clean.(models),hp,l))
                end

                push!(udf,DataFrame(UUID=uuid,e=e,u=u,modelresults=mrs))
            end
        end
    end

    return vcat(udf...),[]
end

function processimage(dataset::Dict,resultroot,env::Dict{String,Any};uuid="",log::Dict{String,AbstractLogger}=Dict{String,AbstractLogger}(),delay=20,binwidth=10,minpredur=10,mincondtest=12000,isplot=true,
    nscale=2,downsample=2,sigma=1.5,pixelscale=255,epprimagescaleindex=2,epprndelay=1,epprnft=[4,4,4],epprlambda=10,fitmodels=[:ePPR])

    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    subject=ex["Subject_ID"];recordsession=ex["RecordSession"];recordsite=ex["RecordSite"]
    bgcolor=RGBA(getparam(envparam,"BGColor")...)
    imagesetname = replace(getparam(envparam,"ImageSet","ImageQuad"),"Ã—"=>"_")
    imagemasktype = getparam(envparam,"MaskType","ImageQuad")
    imagemaskradius = getparam(envparam,"MaskRadius","ImageQuad")
    imagemasksigma = getparam(envparam,"Sigma","ImageQuad")
    ct,ctc = ctctc(ex)
    vcti = .!(isnan.(ct[:CondOn]) .& isnan.(ct[:CondOff]));ct=ct[vcti,:];ctc=ctc[vcti,:]

    predur = max(preicidur,minpredur)
    resultroot=abspath(resultroot,ex["ID"])
    !isdir(resultroot) && mkpath(resultroot)
    logfile=joinpath(resultroot,"batchlog.txt")
    # ignore dataset with less than mincondtest condtest
    if nrow(ct) < mincondtest
        log(ex["ID"],()->@info("number of condition tests($(nrow(ct))) in dataset($uuid) is less than $mincondtest"),logfile)
        return nothing,nothing
    end

    # load and prepare imageset in env
    if !haskey(env,imagesetname) && haskey(env,"rootdir") && isdir(env["rootdir"])
        pyramid = Dict{Symbol,Any}(:pyramid => map(i->gaussian_pyramid(i, nscale-1, downsample, sigma),
        loadimageset(joinpath(env["rootdir"],imagesetname),alpha=true)))
        pyramid[:size] = map(i->size(i),pyramid[:pyramid][1])
        env[imagesetname] = pyramid
    end
    # recover image stimuli
    imageset = env[imagesetname]
    bgimagecolor = oftype(imageset[:pyramid][1][1][1],bgcolor)
    unmaskindex = map(i->alphamask(i,radius=imagemaskradius,sigma=imagemasksigma,masktype=imagemasktype)[2],imageset[:pyramid][1])
    imagestimuli = map(s->map(i->alphablend.(alphamask(i[s],radius=imagemaskradius,sigma=imagemasksigma,masktype=imagemasktype)[1],[bgimagecolor]),imageset[:pyramid]),1:nscale)

    udf = []
    # spike triggered average model
    if :STA in fitmodels
        scaleindex=1
        xsize=imageset[:size][scaleindex]
        x = Array{Float64}(undef,size(ctc,1),prod(xsize))
        for i in 1:size(ctc,1)
            x[i,:]=vec(gray.(imagestimuli[scaleindex][Int(ctc[:Image][i])]))
        end
        delays=20:30:180
        if haskey(dataset,"spike")
            spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
            for e in spikeeid
                ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
                for u in euuid
                    cellid = join([subject,recordsession,recordsite,"E$e","U$u"],"_")
                    udir = joinpath(resultroot,"$(uuid)_E$(e)_U$(u)_STA_S$(scaleindex)")
                    isdir(udir) && rm(udir,recursive=true)
                    ustas=[]
                    for d in delays
                        _,y,_,_ = subrv(est[esu.==u],ct[:CondOn].+d,ct[:CondOff].+d,israte=false)
                        # ignore irresponsive unit
                        if sum(y)==0
                            log(ex["ID"],"no $(d)ms delayed responses of cell($cellid)",logfile)
                            continue
                        end
                        r=sta(x,y,xsize,decor=false)
                        if isplot
                            plotsta(r,delay=d,decor=false,savedir=udir)
                        end
                        push!(ustas,(d,r))
                    end
                    if !isempty(ustas)
                        push!(udf,DataFrame(UUID=uuid,CellID=cellid,sta=ustas))
                    end
                end
            end
        end
    end

    # eppr model
    if :ePPR in fitmodels
        xsize=imageset[:size][epprimagescaleindex]
        xi = unmaskindex[epprimagescaleindex]
        x = Array{Float64}(undef,size(ctc,1),prod(xsize))
        for i in 1:size(ctc,1)
            x[i,:]=vec(gray.(imagestimuli[epprimagescaleindex][Int(ctc[:Image][i])]))*pixelscale
        end
        if haskey(dataset,"spike")
            spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
            for e in spikeeid
                ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
                for u in euuid
                    cellid = join([subject,recordsession,recordsite,"E$e","U$u"],"_")
                    y = subrvr(est[esu.==u],ct[:CondOn].+delay,ct[:CondOff].+delay)
                    # ignore irresponsive unit
                    if sum(y)==0
                        log(ex["ID"],"no $(delay)ms delayed responses of cell($cellid)",logfile)
                        continue
                    end
                    udir = joinpath(resultroot,"$(uuid)_E$(e)_U$(u)_ePPR_S$(epprimagescaleindex)_D$(epprndelay)")
                    isdir(udir) && rm(udir,recursive=true)

                    debug = isplot ? ePPRDebugOptions(level=DebugVisual,logdir=udir) : ePPRDebugOptions()
                    hp = ePPRHyperParams(xsize...,xindex=xi,ndelay=epprndelay,blankcolor=gray(bgimagecolor)*pixelscale,nft=epprnft,lambda=epprlambda)
                    model,models = epprhypercv(x,y,hp,debug)

                    if isplot && model!=nothing
                        debug(plotmodel(model,hp),log="Model_Final (λ=$(hp.lambda))")
                    end

                    push!(udf,DataFrame(UUID=uuid,CellID=cellid,model=clean!(model),models=clean!(models),hp=hp))
                end
            end
        end
    end

    return isempty(udf) ? nothing : vcat(udf...),nothing
end

function processori(dataset::Dict,resultroot;uuid="",delay=20,log::Dict{String,AbstractLogger}=Dict{String,AbstractLogger}(),binwidth=10,minpredur=100,minrepeat=3,minfiringrate=5,isplot=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    subject=ex["Subject_ID"];recordsession=ex["RecordSession"];recordsite=ex["RecordSite"]
    ct,ctc=ctctc(ex);vcti = .!(isnan.(ct[:CondOn]) .& isnan.(ct[:CondOff]));ct=ct[vcti,:];ctc=ctc[vcti,:]
    cond=condin(ctc)

    resultroot=abspath(resultroot,ex["ID"])
    !isdir(resultroot) && mkpath(resultroot)
    logfile=joinpath(resultroot,"batchlog.txt")
    # ignore dataset with less than minrepeat for any condition
    if any(cond[:n].<minrepeat)
        log(ex["ID"],()->@info("dataset($uuid) doesn't have enough condition repetition",condrepeat=cond[:n]),logfile)
        return nothing,nothing
    end
    ff = finalfactor(cond)[1]
    oris = ctc[ff]
    predur = max(preicidur,minpredur)

    udf = []
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                cellid = join([subject,recordsession,recordsite,"E$e","U$u"],"_")
                preurs = subrvr(est[esu.==u],ct[:CondOn].+(delay-predur),ct[:CondOn].+delay)
                urs = subrvr(est[esu.==u],ct[:CondOn].+delay,ct[:CondOff].+delay)
                # ignore irresponsive unit
                if !isresponsive(preurs,urs,cond[:i])
                    log(ex["ID"],"all condition responses of cell($cellid) are insignificantly different from baseline responses",logfile)
                    continue
                end
                # ignore when maximum mean response less than minfiringrate
                mse = condresponse(urs,cond[:i])
                if all(mse[:m].<minfiringrate)
                    log(ex["ID"],"all mean condition responses of cell($cellid) are less than $minfiringrate spike/s",logfile)
                    continue
                end

                stats = statsori(Float64.(oris),Float64.(urs))
                if isplot
                    plotname = "$(uuid)_E$(e)_U$(u)"
                    plotcondresponse(urs,ctc,ff,u=u,title=plotname,legend=:none)
                    png(joinpath(resultroot,plotname))
                end

                push!(udf,[DataFrame(UUID=uuid,CellID=cellid) DataFrame(stats)])
            end
        end
    end

    return isempty(udf) ? nothing : vcat(udf...),nothing
end

function processlaser(dataset::Dict,resultroot;uuid="",delay=20,binwidth=10,minpredur=100,isplot=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    ct,ctc=ctctc(ex)
    cond=condin(ctc)

    predur = max(preicidur,minpredur)
    resultroot=abspath(resultroot,ex["ID"])
    !isdir(resultroot) && mkpath(resultroot)

    udf = []
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                preurs = subrvr(est[esu.==u],ct[:CondOn]+delay-predur,ct[:CondOn]+delay)
                urs = subrvr(est[esu.==u],ct[:CondOn]+delay,ct[:CondOff]+delay)
                !isresponsive(preurs,urs,cond[:i]) && continue
                ksadp = pvalue(KSampleADTest(getindex.([urs],cond[:i])...))

                if isplot
                    for f in finalfactor(ctc)
                        plotname = "$(uuid)_E$(e)_U$(u)_$f"
                        plotcondresponse(urs,ctc,u,factor=f,title=plotname,legend=:none)
                        png(joinpath(resultroot,plotname))
                    end
                end

                push!(udf,DataFrame(UUID=uuid,e=e,u=u,ksadp=ksadp))
            end
        end
    end

    cond[:UUID]=uuid
    return vcat(udf...),cond
end

function processflash(files,param;uuid="",log::Dict{Any,AbstractLogger}=Dict{Any,AbstractLogger}(),isplot=true)
    dataset = prepare(abspath(param[:dataexportroot],files))

    # Condition Tests
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    minconddur=minimum(condoff-condon)

    subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    sitedir = joinpath(param[:resultroot],subject,siteid)
    testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    resultdir = joinpath(sitedir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # Prepare LFP
    lfregex = Regex("^$(uppercase(testid))[A-Za-z0-9_]*.imec.lf.bin")
    lffile = joinpath(param[:dataroot],filter(i->occursin(lfregex,i),readdir(param[:dataroot]))[1])
    nsavech=dataset["lf"]["meta"]["nSavedChans"]
    nsample=dataset["lf"]["meta"]["nFileSamp"]
    fs=dataset["lf"]["meta"]["fs"]
    nch = dataset["lf"]["meta"]["snsApLfSy"][2]
    refchmask = refchmaskim(dataset["ap"]["meta"])
    nrow,ncol = size(refchmask)
    mmlfp = Mmap.mmap(lffile,Matrix{Int16},(nsavech,nsample),0)

    # Laminar CSD
    epochdur = 0.15
    epoch = [0 epochdur]
    epochs = condon.+epoch
    # all epoch LFP segments for each channel of probe in shape
    ys=reshape2ref(subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),refchmask)

    if isplot
        for col in 1:ncol
            mcys=dropdims(mean(ys[:,col,:,:],dims=3),dims=3)
            plotanalog(mcys,fs=fs)
            foreach(i->savefig(joinpath(resultdir,"column_$(col)_lfp$i")),[".png",".svg"])

            mccsd = dropdims(mean(csd(ys[:,col,:,:]),dims=3),dims=3)
            plotanalog(imfilter(mccsd,Kernel.gaussian((1,1))),fs=fs)
            foreach(i->savefig(joinpath(resultdir,"column_$(col)_csd$i")),[".png",".svg"])
        end
    end

    # Column Mean Normalized CSD
    pys = cat((ys[:,i,:,:] for i in 1:ncol)...,dims=3)
    pcsd = csd(pys)
    baseline = pcsd[:,epoch2samplerange([0 0.015],fs),:] # 0-15ms csd as baseline
    ncsd=pcsd.-mean(baseline,dims=2) # Δcsd relative to baseline
    mncsd = dropdims(mean(ncsd,dims=3),dims=3)
    mncsd[1,:].=0;mncsd[end,:].=0
    mncsd = imfilter(mncsd,Kernel.gaussian((2,1)))
    h=20;depths = h*(1:size(mncsd,1))
    if isplot
        plotanalog(mncsd,fs=fs,y=depths,color=:RdBu)
        foreach(i->savefig(joinpath(resultdir,"column_mean_normalized_csd$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"csd.jld2"),"csd",mncsd,"depth",depths,"fs",fs)

    # Depth Power Spectrum
    epochdur = 0.2
    epoch = [0 epochdur]
    epochs = condon.+epoch
    ys=reshape2ref(subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),refchmask)
    pys = cat((ys[:,i,:,:] for i in 1:ncol)...,dims=3)
    ps,freq = powerspectrum(pys,fs)

    epoch = [-epochdur 0]
    epochs = condon.+epoch
    ys=reshape2ref(subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),refchmask)
    pys = cat((ys[:,i,:,:] for i in 1:ncol)...,dims=3)
    bps,freq = powerspectrum(pys,fs)

    rcps = ps./bps.-1
    mrcps = imfilter(dropdims(mean(rcps,dims=3),dims=3),Kernel.gaussian((1,1)))
    if isplot
        mps = dropdims(mean(ps,dims=3),dims=3)
        plotanalog(imfilter(mps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:db,color=:fire)
        foreach(i->savefig(joinpath(resultdir,"depth_powerspectrum$i")),[".png",".svg"])
        mbps =dropdims(mean(bps,dims=3),dims=3)
        plotanalog(imfilter(mbps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:db,color=:fire)
        foreach(i->savefig(joinpath(resultdir,"depth_powerspectrum_baseline$i")),[".png",".svg"])

        plotanalog(mrcps,x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:db,color=:fire)
        foreach(i->savefig(joinpath(resultdir,"depth_powerspectrum_rc$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"powerspectrum.jld2"),"rcps",mrcps,"depth",depths,"freq",freq)

    unitgood=trues(length(unitspike))
    # Unit Position
    layer = haskey(param,:layer) ? param[:layer] : nothing
    if isplot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(i->savefig(joinpath(resultdir,"Unit_Position$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"spike.jld2"),"spike",spike)

    # Unit Spike Trian
    if isplot
        epochext = preicidur
        for u in 1:length(unitspike)
            ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            plotspiketrain(ys,timeline=[0,minconddur],title="Unit_$(unitid[u])")
            foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_SpikeTrian$i")),[".png",".svg"])
        end
    end

    # Unit Depth PSTH
    epochdur = 0.15
    epoch = [0 epochdur]
    epochs = condon.+epoch
    bw = 0.002
    psthbins = epoch[1]:bw:epoch[2]

    baseindex = epoch2samplerange([0 0.015],1/bw) # 0-15ms of psth as baseline
    unitpsth = map(ust->psth(subrv(ust,epochs,isminzero=true)[1],psthbins,israte=true,normfun=x->x.-mean(x[baseindex])),unitspike)
    depthpsth,x,depths = spacepsth(unitpsth,unitposition,spacebinedges=h*(0:nrow+1))
    depthpsth = imfilter(depthpsth,Kernel.gaussian((1,1)))
    if isplot
        plotpsth(depthpsth,x,depths)
        foreach(i->savefig(joinpath(resultdir,"depthpsth$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"depthpsth.jld2"),"depthpsth",depthpsth,"depth",depths,"x",x)

    # Single Unit Binary Spike Trian of Conditions
    bepochext = -0.3
    bepoch = [-bepochext minconddur]
    bepochs = condon.+bepoch
    spikebins = bepoch[1]:0.001:bepoch[2]
    subst = map(ust->float.(histmatrix(subrv(ust,bepochs,isminzero=true,shift=bepochext)[1],spikebins)[1])',unitspike[unitgood])
    # Single Unit Correlogram and Circuit
    lag=50;suid = unitid[unitgood]
    ccgs,x,ccgis,projs,eunits,iunits = circuitestimate(subst,lag=lag)

    if isplot
        for i in 1:length(ccgs)
            title = "Correlogram between Unit $(suid[ccgis[i][1]]) and $(suid[ccgis[i][2]])"
            bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
            foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
        end
    end
    save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits)

    return nothing,nothing
end

function processhartley(files,param;uuid="",log::Dict{Any,AbstractLogger}=Dict{Any,AbstractLogger}(),isplot=true)
    dataset = prepare(abspath(param[:dataexportroot],files))

    # Condition Tests
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood=spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    minconddur=minimum(condoff-condon)

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    factors = finalfactor(ctc)
    blank = haskey(param,:blank) ? param[:blank] : (:Ori,NaN)
    ci = ctc[blank[1]].!=blank[2]
    cctc = ctc[ci,factors]
    ccondidx = condidx[ci]
    ccondon = condon[ci]
    ccondoff = condoff[ci]

    bi = .!ci
    bctc = ctc[bi,factors]
    bcondidx = condidx[bi]
    bcondon = condon[bi]
    bcondoff = condoff[bi]

    subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    resultdir = joinpath(param[:resultroot],subject,siteid,testid)
    isdir(resultdir) || mkpath(resultdir)

    unitgood=trues(length(unitspike))
    # Unit Position
    layer = haskey(param,:layer) ? param[:layer] : nothing
    if isplot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(i->savefig(joinpath(resultdir,"unitposition$i")),[".png",".svg"])
    end

    # Prepare Imageset
    imgfile = joinpath(param[:dataroot],"$(testid)_image.mat")
    imgs = readmat(imgfile)["imageArray"]
    mat2julia!(imgs)
    px=613;py=135;rp=32

    rawimageset = map(i->dropdims(mean(i[py-rp:py+rp,px-rp:px+rp,:],dims=3),dims=3),imgs)
    imageset = map(i->gaussian_pyramid(i, 1, 4, 1.5)[2],rawimageset)

    imagesize = size(imageset[1])
    x = Array{Float64}(undef,length(ccondidx),prod(imagesize))
    foreach(i->x[i,:]=vec(imageset[ccondidx[i]]),1:size(x,1))

    if :STA in param[:model]
        for u in 1:length(unitspike), d in 0.08:0.01:0.08
            _,y,_,_ = subrv(unitspike[u],ccondon.+d,ccondoff.+d,israte=false)
            r=sta(x,y)
            if isplot
                plotsta(r,imagesize,delay=d)
                foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_STA_$d$i")),[".png"])
            end
        end
    end

    if :ePPR in param[:model]
        for u in 75:75, d in 0.08:0.01:0.08
            cellid = "$(siteid)_U$(unitid[u])"
            unitresultdir = joinpath(resultdir,"Unit_$(unitid[u])_ePPR_$d")
            isdir(unitresultdir) && rm(unitresultdir,recursive=true)
            y = subrvr(unitspike[u],ccondon.+d,ccondoff.+d)

            debug = isplot ? ePPRDebugOptions(level=DebugVisual,logdir=unitresultdir) : ePPRDebugOptions()
            hp = ePPRHyperParams(imagesize...,ndelay=param[:epprndelay],blankcolor=128,nft=param[:epprnft],lambda=param[:epprlambda])
            model,models = epprhypercv(x,y,hp,debug)

            if isplot && !isnothing(model)
                debug(plotmodel(model,hp),log="Model_Final (λ=$(hp.lambda))")
            end
        end
    end

    return nothing,nothing
end

function processcondtest(files,param;uuid="",log::Dict{Any,AbstractLogger}=Dict{Any,AbstractLogger}(),isplot=true)
    dataset = prepare(abspath(param[:dataexportroot],files))

    # Condition Tests
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    minconddur=minimum(condoff-condon)

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    factors = finalfactor(ctc)
    blank = haskey(param,:blank) ? param[:blank] : (:Ori,"blank")
    ci = ctc[blank[1]].!=blank[2]
    cctc = ctc[ci,factors]
    ccondon=condon[ci]
    ccondoff=condoff[ci]

    bi = .!ci
    bctc = ctc[bi,factors]
    bcondon=condon[bi]
    bcondoff=condoff[bi]

    subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    sitedir = joinpath(param[:resultroot],subject,siteid)
    testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    resultdir = joinpath(sitedir,testid)
    isdir(resultdir) || mkpath(resultdir)

    unitgood=trues(length(unitspike))
    # Unit Position
    layer = haskey(param,:layer) ? param[:layer] : nothing
    if isplot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(i->savefig(joinpath(resultdir,"Unit_Position$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"spike.jld2"),"spike",spike)

    # Unit Spike Trian
    if isplot
        epochext = preicidur
        for u in 1:length(unitspike)
            ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            plotspiketrain(ys,timeline=[0,minconddur],title="Unit_$(unitid[u])")
            foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_SpikeTrian$i")),[".png",".svg"])
        end
    end

    # Condition Response
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : 0.015
    for u in 1:length(unitspike)
        rs = subrvr(unitspike[u],ccondon.+responsedelay,ccondoff.+responsedelay)
        if isplot
            plotcondresponse(rs,cctc)
            foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_CondResponse$i")),[".png",".svg"])
        end
    end

    # blank condition response
    ubr=[]
    if !isempty(bcondon)
        for u in 1:length(unitspike)
            br = subrvr(unitspike[u],bcondon.+responsedelay,bcondoff.+responsedelay)
            mseuc = condresponse(br,condin(bctc))
            push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
        end
    end

    # Condition Response in Factor Space
    fms,fses,fa=factorresponse(unitspike,cctc,ccondon,ccondoff,responsedelay=responsedelay)

    ufs = Dict(k=>[] for k in keys(fa));optconds=[]
    for u in 1:length(fms)
        p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
        push!(optconds, OrderedDict(map((f,i)->f=>fa[f][i],keys(fa),p)))
        for f in keys(fa)
            fd = findfirst(f.==keys(fa))
            fdn = length(fa[f])
            fi = copy(p)
            fi[fd]=1:fdn
            mseuc=DataFrame(m=fms[u][fi...],se=fses[u][fi...],u=fill(unitid[u],fdn))
            mseuc[f]=fa[f]
            push!(ufs[f],factorresponsestats(mseuc[f],mseuc[:m],factor=f))
            if isplot
                plotcondresponse(dropmissing(mseuc),colors=[:black],responseline=isempty(ubr) ? [] : ubr[u:u])
                foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_$(f)_Tuning$i")),[".png",".svg"])
            end
        end
    end
    save(joinpath(resultdir,"factorresponse.jld2"),"factorstats",ufs,"fms",fms,"fses",fses,"fa",fa)

    corrcond = haskey(param,:corrcond) ? param[:corrcond] : :allcond
    if typeof(corrcond) <: Union{Dict,NamedTuple}
        corrcondi = findcond(cctc,corrcond)
    elseif corrcond == :bestcond
        corrcondi = mapfoldl(uoc->findcond(cctc,uoc),union,optconds)
    else
        corrcondi = 1:length(ccondon)
    end
    # Single Unit Binary Spike Trian of Conditions
    bepochext = -0.3
    bepoch = [-bepochext minconddur]
    bepochs = ccondon[corrcondi].+bepoch
    spikebins = bepoch[1]:0.001:bepoch[2]
    subst = map(ust->float.(histmatrix(subrv(ust,bepochs,isminzero=true,shift=bepochext)[1],spikebins)[1])',unitspike[unitgood])
    # Single Unit Correlogram and Circuit
    lag=50;suid = unitid[unitgood]
    ccgs,x,ccgis,projs,eunits,iunits = circuitestimate(subst,lag=lag)

    if isplot
        for i in 1:length(ccgs)
            title = "Correlogram between Unit $(suid[ccgis[i][1]]) and $(suid[ccgis[i][2]])"
            bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
            foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
        end
    end
    save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits)

    return nothing,nothing
end
