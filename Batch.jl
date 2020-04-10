using NeuroAnalysis,Query,FileIO,ProgressMeter,Logging,Statistics,DataFrames,Plots,Mmap,Images,StatsBase#,ePPR
using DataFramesMeta,StatsPlots,Interact,CSV,MAT,DataStructures,HypothesisTests,StatsFuns,Random
import Base: close

function batchtests(tests::DataFrame,param::Dict{Any,Any}=Dict{Any,Any}();log::Dict{Any,AbstractLogger}=Dict{Any,AbstractLogger}(),plot::Bool=true)
    p = ProgressMeter.Progress(size(tests,1),desc="Batch Tests ... ")
    for t in eachrow(tests)
        try
            if t[:ID]=="OriGrating"
                u,c=processori(dataset,resultroot,uuid=uuid,log=log,delay=delay,binwidth=binwidth,plot=plot)
            elseif t[:ID]=="Laser"
                u,c=processlaser(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,plot=plot)
            elseif t[:ID]=="Image"
                u,c=processimage(dataset,resultroot,env,uuid=uuid,log=log,delay=delay,binwidth=binwidth,plot=plot)
            elseif t[:ID]=="LaserImage"
                u,c=processlaserimage(dataset,condroot,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,plot=plot)
            elseif t[:ID] in ["Flash","Flash2Color"] && t[:sourceformat]=="SpikeGLX"
                process_flash(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["Hartley","HartleySubspace"] && t[:sourceformat]=="SpikeGLX"
                process_hartley(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["OriSF","OriSFColor","Color"] && t[:sourceformat]=="SpikeGLX"
                process_condtest(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["DirSF"] && t[:sourceformat]=="Scanbox"
                process_2P_dirsf(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["DirSFColor"] && t[:sourceformat]=="Scanbox"
                process_2P_dirsfcolor(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
            elseif t[:ID] in ["Hartley"] && t[:sourceformat]=="Scanbox"
                # process_2P_hartleySTA(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
                process_2P_hartleyFourier(t[:files],param,uuid=t[:UUID],log=log,plot=plot)
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
    nscale=2,downsample=2,sigma=1.5,pixelscale=255,plot=true,forceprocess=true)
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
                plotdir = plot ? udir : nothing
                mrs = []
                for l in eachrow(lcond)
                    debug = plot ? ePPRDebugOptions(level=DebugVisual,logdir=joinpath(plotdir,condstring(l,lfactor))) : ePPRDebugOptions()
                    hp = ePPRHyperParams(ximagesize...,xindex=xi,ndelay=1,blankcolor=gray(bgimagecolor)*pixelscale)
                    hp.nft = [6]
                    hp.lambda = 30000
                    model,models = epprcv(imagestimulimatrix[Int.(ctc[:Image][l[:i]]),:]*pixelscale,y[l[:i]],hp,debug)

                    if plot && model!=nothing
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

function processimage(dataset::Dict,resultroot,env::Dict{String,Any};uuid="",log::Dict{String,AbstractLogger}=Dict{String,AbstractLogger}(),delay=20,binwidth=10,minpredur=10,mincondtest=12000,plot=true,
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
                        if plot
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

                    debug = plot ? ePPRDebugOptions(level=DebugVisual,logdir=udir) : ePPRDebugOptions()
                    hp = ePPRHyperParams(xsize...,xindex=xi,ndelay=epprndelay,blankcolor=gray(bgimagecolor)*pixelscale,nft=epprnft,lambda=epprlambda)
                    model,models = epprhypercv(x,y,hp,debug)

                    if plot && model!=nothing
                        debug(plotmodel(model,hp),log="Model_Final (λ=$(hp.lambda))")
                    end

                    push!(udf,DataFrame(UUID=uuid,CellID=cellid,model=clean!(model),models=clean!(models),hp=hp))
                end
            end
        end
    end

    return isempty(udf) ? nothing : vcat(udf...),nothing
end

function processori(dataset::Dict,resultroot;uuid="",delay=20,log::Dict{String,AbstractLogger}=Dict{String,AbstractLogger}(),binwidth=10,minpredur=100,minrepeat=3,minfiringrate=5,plot=true)
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
                if plot
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

function processlaser(dataset::Dict,resultroot;uuid="",delay=20,binwidth=10,minpredur=100,plot=true)
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

                if plot
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

function process_flash(files,param;uuid="",log=nothing,plot=true)
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
    envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood=spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    minconddur=minimum(condoff-condon)

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)

    # Prepare LFP
    lffile = matchfile(Regex("^$(testid)[A-Za-z0-9_]*.imec.lf.bin"),dir = datadir,adddir=true)[1]
    nsavech=dataset["lf"]["meta"]["nSavedChans"]
    nsample=dataset["lf"]["meta"]["nFileSamp"]
    fs=dataset["lf"]["meta"]["fs"]
    nch = dataset["lf"]["meta"]["snsApLfSy"][2]
    hx,hy,hz = dataset["lf"]["meta"]["probespacing"]
    badchmask = badchmaskim(dataset)
    pnrow,pncol = size(badchmask)
    mmlf = Mmap.mmap(lffile,Matrix{Int16},(nsavech,nsample),0)

    # Depth LFP and CSD
    epochdur = timetounit(150)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    # All LFP epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered, in the same shape of probe, where bad channels are replaced with local average
    ys=reshape2mask(subrm(mmlf,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"],bandpass=[1,100]),badchmask)

    if plot
        for c in 1:pncol
            mcys=dropdims(mean(ys[:,c,:,:],dims=3),dims=3)
            plotanalog(mcys,fs=fs,cunit=:uv)
            foreach(i->savefig(joinpath(resultdir,"Column_$(c)_MeanLFP$i")),[".png",".svg"])

            mccsd = dropdims(mean(csd(ys[:,c,:,:],h=hy),dims=3),dims=3)
            plotanalog(imfilter(mccsd,Kernel.gaussian((1,1))),fs=fs)
            foreach(i->savefig(joinpath(resultdir,"Column_$(c)_MeanCSD$i")),[".png",".svg"])
        end
    end

    # Column Mean ΔCSD
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    pcsd = csd(pys,h=hy)
    basedur = timetounit(15)
    baseindex = epoch2samplerange([0 basedur],fs)
    dcsd=stfilter(pcsd,temporaltype=:sub,ti=baseindex,hedgevalue=0)
    mdcsd = dropdims(mean(dcsd,dims=3),dims=3)
    depths = hy*(1:size(mdcsd,1))
    if plot
        plotanalog(imfilter(mdcsd,Kernel.gaussian((1,1))),fs=fs,y=depths,color=:RdBu)
        foreach(i->savefig(joinpath(resultdir,"Columns_MeandCSD$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"csd.jld2"),"csd",mdcsd,"depth",depths,"fs",fs,"log",ex["Log"],"color","$(ex["Param"]["ColorSpace"])_$(ex["Param"]["Color"])")

    # Depth Power Spectrum
    epochdur = timetounit(500)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    ys=reshape2mask(subrm(mmlf,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"],bandpass=[1,100]),badchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=4*epochdur*SecondPerUnit)

    epoch = [-epochdur 0]
    epochs = condon.+epoch
    ys=reshape2mask(subrm(mmlf,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"],bandpass=[1,100]),badchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    bps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=4*epochdur*SecondPerUnit)

    rcps = ps./bps.-1
    mrcps = dropdims(mean(rcps,dims=3),dims=3)
    if plot
        mps = dropdims(mean(ps,dims=3),dims=3)
        plotanalog(imfilter(mps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
        foreach(i->savefig(joinpath(resultdir,"DepthPowerSpectrum$i")),[".png",".svg"])
        mbps =dropdims(mean(bps,dims=3),dims=3)
        plotanalog(imfilter(mbps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
        foreach(i->savefig(joinpath(resultdir,"DepthPowerSpectrum_Baseline$i")),[".png",".svg"])

        plotanalog(imfilter(mrcps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:v2,color=:PuRd)
        foreach(i->savefig(joinpath(resultdir,"DepthPowerSpectrum_RelativeChange$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"powerspectrum.jld2"),"rcps",mrcps,"depth",depths,"freq",freq)

    # Unit Position
    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    if plot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(i->savefig(joinpath(resultdir,"UnitPosition$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"spike.jld2"),"spike",spike)

    # Unit Spike Trian
    if plot
        epochext = preicidur
        for u in eachindex(unitspike)
            ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            plotspiketrain(ys,timeline=[0,conddur],title="$(ugs[u])Unit_$(unitid[u])")
            foreach(i->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_SpikeTrian$i")),[".png",".svg"])
        end
    end

    # Unit Depth PSTH
    epochdur = timetounit(150)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    bw = timetounit(2)
    psthbins = epoch[1]:bw:epoch[2]

    baseindex = epoch2samplerange([0 basedur],1/(bw*SecondPerUnit))
    normfun = x->x.-mean(x[baseindex])
    unitpsth = map(ust->psth(subrv(ust,epochs,isminzero=true)[1],psthbins,israte=true,normfun=normfun),unitspike)
    depthpsth,x,depths,depthnunit = spacepsth(unitpsth,unitposition,spacebinedges=hy*(0:pnrow+1))
    if plot
        plotpsth(imfilter(depthpsth,Kernel.gaussian((1,1))),x,depths,color=:minmax,n=depthnunit)
        foreach(i->savefig(joinpath(resultdir,"DepthPSTH$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"depthpsth.jld2"),"depthpsth",depthpsth,"depth",depths,"x",x,"n",depthnunit)

    # Single Unit Binary Spike Trian of Condition Tests
    bepochext = timetounit(-300)
    bepoch = [-bepochext minconddur]
    bepochs = condon.+bepoch
    spikebins = bepoch[1]:timetounit(1):bepoch[2] # 1ms bin histogram will convert spike times to binary spike trains
    subst = map(ust->float.(histmatrix(subrv(ust,bepochs,isminzero=true,shift=bepochext)[1],spikebins)[1])',unitspike[unitgood]) # nSpikeBin x nEpoch
    # Single Unit Correlogram and Circuit
    lag=50
    ccgs,x,ccgis,projs,eunits,iunits,projweights = circuitestimate(subst,lag=lag,unitid=unitid[unitgood],condis=cond.i)

    if !isempty(projs)
        if plot
            for i in eachindex(ccgs)
                title = "Correlogram between single unit $(ccgis[i][1]) and $(ccgis[i][2])"
                bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
                foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
            end
            plotcircuit(unitposition,unitid,projs,unitgood=unitgood,eunits=eunits,iunits=iunits,layer=layer)
            foreach(i->savefig(joinpath(resultdir,"UnitPosition_Circuit$i")),[".png",".svg"])
        end
        save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits,"projweights",projweights)
    end
end

function process_hartley(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))
    process_hartley(dataset,param;uuid=uuid,log=log,plot=plot)
end
function process_hartley(dataset::Dict,param;uuid="",log=nothing,plot=true)
    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    resultsitedir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(resultsitedir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # Condition Tests
    envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood=spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]

    # Prepare Conditions
    condtable = DataFrame(ex["Cond"])
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    blank = haskey(param,:blank) ? param[:blank] : (:Ori_Final,NaN)
    ci = ctc[!,blank[1]].!==blank[2]
    cctc = ctc[ci,factors]
    ccondidx = condidx[ci]
    ccondon = condon[ci]
    ccondoff = condoff[ci]
    bi = .!ci
    bctc = ctc[bi,factors]
    bcondidx = condidx[bi]
    bcondon = condon[bi]
    bcondoff = condoff[bi]
    isblank = !isempty(bcondon)

    # Unit Position
    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    if plot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(i->savefig(joinpath(resultdir,"UnitPosition$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"spike.jld2"),"spike",spike)



    # Prepare Imageset
    nscale = haskey(param,:nscale) ? param[:nscale] : 2
    downsample = haskey(param,:downsample) ? param[:downsample] : 2
    sigma = haskey(param,:sigma) ? param[:sigma] : 1.5
    bgcolor = RGBA(getparam(envparam,"BGColor")...)
    maxcolor = RGBA(getparam(envparam,"MaxColor")...)
    mincolor = RGBA(getparam(envparam,"MinColor")...)
    masktype = getparam(envparam,"MaskType")
    maskradius = getparam(envparam,"MaskRadius")
    masksigma = getparam(envparam,"Sigma")
    diameter = 7#getparam(envparam,"Diameter")
    stisize = (diameter,diameter)
    imagesetname = splitext(splitdir(ex["CondPath"])[2])[1] * "_stisize$stisize"
    if !haskey(param,imagesetname)
        imageset = map(i->GrayA.(grating(θ=deg2rad(i.Ori),sf=i.SpatialFreq,phase=i.SpatialPhase,stisize=stisize,ppd=45)),eachrow(condtable))
        imageset = Dict{Symbol,Any}(:pyramid => map(i->gaussian_pyramid(i, nscale-1, downsample, sigma),imageset))
        imageset[:imagesize] = map(i->size(i),imageset[:pyramid][1])
        param[imagesetname] = imageset
    end

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
        imagesize = imageset[:imagesize][scaleindex]
        xi = unmaskindex[scaleindex]
        uci = unique(condidx)
        ucii = map(i->findall(condidx.==i),uci)
        x = Array{Float64}(undef,length(uci),length(xi))
        foreach(i->x[i,:]=gray.(imagestimuli[scaleindex][uci[i]][xi]),1:size(x,1))

        uy=Dict();usta = Dict()
        delays = -10:5:200
        for u in eachindex(unitspike)
            !unitgood[u] && continue
            ys = Array{Float64}(undef,length(delays),length(ucii))
            stas = Array{Float64}(undef,length(delays),length(xi))
            for d in eachindex(delays)
                y = subrvr_ono(unitspike[u],condon.+delays[d],condoff.+delays[d],israte=true,isnan2zero=true)
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
                    p = plotsta(stas[d,:],imagesize=imagesize,stisize=stisize,index=xi,title=title,r=r)
                    foreach(i->save(joinpath(resultdir,"$title$i"),p),[".png"])
                end
            end
        end
        save(joinpath(resultdir,"sta.jld2"),"imagesize",imagesize,"x",x,"xi",xi,"xcond",condtable[uci,:],"uy",uy,"usta",usta,"delays",delays,
        "stisize",stisize,"log",ex["Log"],"color","$(ex["Param"]["ColorSpace"])_$(ex["Param"]["Color"])","maxcolor",maxcolor,"mincolor",mincolor)
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

function process_condtest(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    resultsitedir = joinpath(param[:resultroot],subject,siteid)
    # testid = join(filter(!isempty,[siteid,ex["TestID"]]),"_")
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(resultsitedir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # Condition Tests
    envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    spike = dataset["spike"];unitspike = spike["unitspike"];unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition=spike["unitposition"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    minconddur=minimum(condoff-condon)

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    if ex["ID"] == "Color"
        factors = [:HueAngle]
        blank = (:Color,36)
    else
        blank = haskey(param,:blank) ? param[:blank] : (:SpatialFreq,0)
    end
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

    # Unit Position
    ugs = map(i->i ? "Single-" : "Multi-",unitgood)
    layer = haskey(param,:layer) ? param[:layer] : nothing
    if plot
        plotunitposition(unitposition,unitgood=unitgood,chposition=spike["chposition"],unitid=unitid,layer=layer)
        foreach(i->savefig(joinpath(resultdir,"UnitPosition$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"spike.jld2"),"spike",spike)

    # Unit Spike Trian
    if plot
        epochext = preicidur
        for u in eachindex(unitspike)
            ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            title = "$(ugs[u])Unit_$(unitid[u])_SpikeTrian"
            plotspiketrain(ys,timeline=[0,conddur],title=title)
            foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
        end
    end

    # Condition Response
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : timetounit(15)
    minresponse = haskey(param,:minresponse) ? param[:minresponse] : 5
    ubr=[];uresponsive=[];umodulative=[]
    for u in eachindex(unitspike)
        rs = subrvr(unitspike[u],ccondon.+responsedelay,ccondoff.+responsedelay,israte=true)
        prs = subrvr(unitspike[u],ccondon.+(responsedelay-preicidur),ccondon.+responsedelay,israte=true)
        push!(uresponsive,isresponsive(prs,rs,ccond.i))
        push!(umodulative,ismodulative([DataFrame(Y=rs) cctc]))
        if isblank
            br = subrvr(unitspike[u],bcondon.+responsedelay,bcondoff.+responsedelay)
            mseuc = condresponse(br,condin(bctc))
            push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
        end
        if plot && length(factors)==2
            plotcondresponse(rs,cctc)
            foreach(i->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_CondResponse$i")),[".png",".svg"])
        end
    end

    # Condition Response in Factor Space
    fms,fses,fa=factorresponse(unitspike,cctc,ccondon,ccondoff,responsedelay=responsedelay)
    pfms,pfses,fa=factorresponse(unitspike,cctc,ccondon.-preicidur,ccondon,responsedelay=responsedelay)
    sfms,sfses,fa=factorresponse(unitspike,cctc,ccondoff,ccondoff.+suficidur,responsedelay=responsedelay)

    ufs = Dict(k=>[] for k in keys(fa));optconds=[];uenoughresponse=[]
    for u in eachindex(fms)
        oi = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
        push!(uenoughresponse,fms[u][oi...]>=minresponse)
        push!(optconds, Dict(map((f,i)->f=>fa[f][i],keys(fa),oi)))
        for f in keys(fa)
            fd = findfirst(f.==keys(fa))
            fdn = length(fa[f])
            fi = copy(oi)
            fi[fd]=1:fdn
            rdf=DataFrame(m=fms[u][fi...],se=fses[u][fi...],u=fill(unitid[u],fdn),ug=fill("$(ugs[u][1])U",fdn))
            rdf[!,f]=fa[f]
            push!(ufs[f],factorresponsestats(rdf[!,f],rdf[!,:m],factor=f))
            if plot
                proj = f in [:Ori,:Ori_Final,:HueAngle] ? :polar : :cartesian
                prdf=DataFrame(m=pfms[u][fi...],se=pfses[u][fi...],u=fill(unitid[u],fdn),ug=fill("Pre_$(ugs[u][1])U",fdn))
                prdf[!,f]=fa[f]
                srdf=DataFrame(m=sfms[u][fi...],se=sfses[u][fi...],u=fill(unitid[u],fdn),ug=fill("Suf_$(ugs[u][1])U",fdn))
                srdf[!,f]=fa[f]
                mseugc = [rdf;prdf;srdf]
                plotcondresponse(dropmissing(mseugc),colors=[:gray40,:crimson,:navyblue],projection=proj,responseline=isblank ? ubr[u:u] : [])
                foreach(i->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_$(f)Tuning$i")),[".png",".svg"])
            end
        end
    end
    if plot
        vi = uresponsive.&umodulative.&unitgood.&uenoughresponse
        pfactor=intersect([:Ori,:Ori_Final],keys(ufs))
        if length(pfactor)>0
            pfactor=pfactor[1]
            colorspace = ex["Param"]["ColorSpace"]
            hues = ex["Param"]["Color"]
            title = "UnitPosition_$(colorspace)_$(hues)_OptOriDirSF"
            p=plotunitpositionproperty(unitposition[vi,:],title=title,ori=map(i->i.oo,ufs[pfactor][vi]),os=map(i->1-i.ocv,ufs[pfactor][vi]),
            dir = map(i->i.od,ufs[pfactor][vi]),ds=map(i->1-i.dcv,ufs[pfactor][vi]),sf=map(i->i.osf,ufs[:SpatialFreq][vi]),layer=layer)
            foreach(i->save(joinpath(resultdir,"$title$i"),p),[".png",".svg"])

            oos = map(i->i.oo,ufs[pfactor][vi])
            oos = [oos;oos.+180]
            ys,ns,ws,is = histrv(oos,0,360,nbins=20)
            oa = deg2rad.(mean.([ws;ws[1]]))
            oh = [ns;ns[1]]./2

            ys,ns,ws,is = histrv(map(i->i.od,ufs[pfactor][vi]),0,360,nbins=20)
            da = deg2rad.(mean.([ws;ws[1]]))
            dh = [ns;ns[1]]

            title = "$(colorspace)_$(hues)_OptOriDirHist"
            Plots.plot([oa da],[oh dh],title=title,projection=:polar,linewidth=2,label=["Ori" "Dir"])
            foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
        end
        if haskey(ufs,:ColorID)
            plotunitposition(unitposition,color=map((i,j)->j ? HSV(i.oh,1-i.hcv,1) : HSVA(0,0,0.5,0.2),ufs[:ColorID],umodulative))
            foreach(i->savefig(joinpath(resultdir,"UnitPosition_Hue_Tuning$i")),[".png",".svg"])
        end
        if haskey(ufs,:HueAngle)
            colorspace = ex["Param"]["ColorSpace"]
            hues = ex["Param"]["Color"]

            ha = ccond[!,:HueAngle]
            hac = map(i->clamp.(cond[findfirst(cond[!,:HueAngle].==i),:Color],0,1),ha)
            ha = deg2rad.(ha)
            title = "UnitPosition_$(colorspace)_$(hues)_OptHue"
            plotunitposition(unitposition[vi,:],color=map(i->RGBA(hac[findclosestangle(deg2rad(i.oh),ha)]...),ufs[:HueAngle][vi]),
            markersize=map(i->(1-i.hcv)*4+2,ufs[:HueAngle][vi]),layer=layer,title=title)
            foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])

            ys,ns,ws,is = histrv(map(i->i.oh,ufs[:HueAngle][vi]),0,360,nbins=20)
            title = "$(colorspace)_$(hues)_OptHueHist"
            Plots.plot(deg2rad.(mean.([ws;ws[1]])),[ns;ns[1]],title=title,projection=:polar,linewidth=2,label="Hue")
            foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
        end
    end
    save(joinpath(resultdir,"factorresponse.jld2"),"factorstats",ufs,"fms",fms,"fses",fses,"pfms",pfms,"pfses",pfses,"sfms",sfms,"sfses",sfses,"fa",fa,
    "responsive",uresponsive,"modulative",umodulative,"enoughresponse",uenoughresponse,"log",ex["Log"],"color","$(ex["Param"]["ColorSpace"])_$(ex["Param"]["Color"])")

    # Single Unit Binary Spike Trian of Condition Tests
    bepochext = timetounit(-300)
    bepoch = [-bepochext minconddur]
    bepochs = ccondon.+bepoch
    spikebins = bepoch[1]:timetounit(1):bepoch[2]
    subst = map(ust->float.(histmatrix(subrv(ust,bepochs,isminzero=true,shift=bepochext)[1],spikebins)[1])',unitspike[unitgood])
    # Single Unit Correlogram and Circuit
    lag=50
    ccgs,x,ccgis,projs,eunits,iunits,projweights = circuitestimate(subst,lag=lag,unitid=unitid[unitgood],condis=ccond.i)

    if !isempty(projs)
        if plot
            for i in eachindex(ccgs)
                title = "Correlogram between single unit $(ccgis[i][1]) and $(ccgis[i][2])"
                bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
                foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
            end
            plotcircuit(unitposition,unitid,projs,unitgood=unitgood,eunits=eunits,iunits=iunits,layer=layer)
            foreach(i->savefig(joinpath(resultdir,"UnitPosition_Circuit$i")),[".png",".svg"])
        end
        save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits,"projweights",projweights)
    end
end

function process_2P_dirsf(files,param;uuid="",log=nothing,plot=false)

    interpolatedData = haskey(param,:interpolatedData) ? param[:interpolatedData] : true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
    preOffset = haskey(param,:preOffset) ? param[:preOffset] : 0.1
    responseOffset = haskey(param,:responseOffset) ? param[:responseOffset] : 0.05  # in sec
    α = haskey(param,:α) ? param[:α] : 0.05   # p value
    sampnum = haskey(param,:sampnum) ? param[:sampnum] : 100   # random sampling 100 times
    fitThres = haskey(param,:fitThres) ? param[:fitThres] : 0.5

    # Expt info
    dataset = prepare(files)
    ex = dataset["ex"]
    disk=ex["source"][1:2];subject=ex["Subject_ID"];recordSession=uppercasefirst(ex["RecordSite"]);testId=ex["TestID"]

    ## Prepare data & result path
    exptId = join(filter(!isempty,[recordSession[2:end], testId]),"_")
    dataFolder = joinpath(disk,subject, "2P_data", recordSession, exptId)

    ## load expt, scanning parameters
    envparam = ex["EnvParam"]
    sbx = dataset["sbx"]["info"]

    ## Align Scan frames with stimulus
    # Calculate the scan parameters
    scanFreq = sbx["resfreq"]
    lineNum = sbx["sz"][1]
    if haskey(sbx, "recordsPerBuffer_bi")
       scanMode = 2  # bidirectional scanning   # if process Splitted data, =1
    else
       scanMode = 1  # unidirectional scanning
    end
    sbxfs = 1/(lineNum/scanFreq/scanMode)   # frame rate
    trialOnLine = sbx["line"][1:2:end]
    trialOnFrame = sbx["frame"][1:2:end] + round.(trialOnLine/lineNum)        # if process splitted data use frame_split
    trialOffLine = sbx["line"][2:2:end]
    trialOffFrame = sbx["frame"][2:2:end] + round.(trialOnLine/lineNum)    # if process splitted data use frame_split

    # On/off frame indces of trials
    trialEpoch = Int.(hcat(trialOnFrame, trialOffFrame))
    # minTrialDur = minimum(trialOffFrame-trialOnFrame)
    # histogram(trialOffFrame-trialOnFrame,nbins=20,title="Trial Duration(Set to $minTrialDur)")

    # Transform Trials ==> Condition
    ctc = DataFrame(ex["CondTestCond"])
    trialNum =  size(ctc,1)
    # Include blank as a condition, marked as Inf
    for factor in 1:size(ctc,2)
        ctc[factor] = replace(ctc[factor], "blank"=>Inf)
    end

    condition = condin(ctc)
    condNum = size(condition,1) # including blanks

    # On/off frame indces of condations/stimuli
    preStim = ex["PreICI"]; stim = ex["CondDur"]; postStim = ex["SufICI"]
    trialOnTime = fill(0, trialNum)
    condOfftime = preStim + stim
    preEpoch = [0 preStim-preOffset]
    condEpoch = [preStim+responseOffset condOfftime-responseOffset]
    preFrame=epoch2samplerange(preEpoch, sbxfs)
    condFrame=epoch2samplerange(condEpoch, sbxfs)
    # preOn = fill(preFrame.start, trialNum)
    # preOff = fill(preFrame.stop, trialNum)
    # condOn = fill(condFrame.start, trialNum)
    # condOff = fill(condFrame.stop, trialNum)

    ## Load data
    segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,adddir=true)[1]
    segment = prepare(segmentFile)
    signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,adddir=true)[1]
    signal = prepare(signalFile)
    sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
    # spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

    planeNum = size(segment["mask"],3)  # how many planes
    planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)

    ## Use for loop process each plane seperately
    for pn in 1:planeNum
        # pn=1  # for test
        # Initialize DataFrame for saving results
        recordPlane = string("00",pn-1)  # plane/depth, this notation only works for expt has less than 10 planes
        siteId = join(filter(!isempty,[recordSession[2:end], testId, recordPlane]),"_")
        dataExportFolder = joinpath(disk,subject, "2P_analysis", recordSession, siteId, "DataExport")
        resultFolder = joinpath(disk,subject, "2P_analysis", recordSession, siteId, "Plots")
        isdir(dataExportFolder) || mkpath(dataExportFolder)
        isdir(resultFolder) || mkpath(resultFolder)
        result = DataFrame()

        cellRoi = segment["seg_ot"]["vert"][pn]
        cellNum = length(cellRoi)
        display("plane: $pn")
        display("Cell Number: $cellNum")
        cellId = collect(range(1, step=1, stop=cellNum))  # Currently, the cellID is the row index of signal

        if interpolatedData
            rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
            # spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        else
            rawF = transpose(signal["sig_ot"]["sig"][pn])
            # spike = transpose(signal["sig_ot"]["spks"][pn])
        end
        result.py = 0:cellNum-1
        result.ani = fill(subject, cellNum)
        result.dataId = fill(siteId, cellNum)
        result.cellId = 1:cellNum

        ## Plot dF/F traces of all trials for all cells
        # Cut raw fluorescence traces according to trial on/off time and calculate dF/F
        cellTimeTrial = sbxsubrm(rawF,trialEpoch,cellId;fun=dFoF(preFrame))  # cellID x timeCourse x Trials
        # Mean response within stim time
        cellMeanTrial = dropdims(mean(cellTimeTrial[:,condFrame,:], dims=2), dims=2)  # cellID x Trials
        # Plot
        if plot
            @manipulate for cell in 1:cellNum
                plotanalog(transpose(cellTimeTrial[cell,:,:]), fs=sbxfs, timeline=condEpoch.-preStim, xunit=:s, ystep=1,cunit=:p, color=:fire,xext=preStim)
            end
        end

        ## Average over repeats, and put each cell's response in factor space (dir-sf...), and find the maximal level of each factor
        factors = finalfactor(ctc)
        fa = OrderedDict(f=>unique(condition[f]) for f in factors)  # factor levels, the last one of each factor maybe blank(Inf)
        fms=[];fses=[];  # mean ans sem of each condition of each cell
        ufm = Dict(k=>[] for k in keys(fa))  # maxi factor level of each cell
        for cell in 1:cellNum
            mseuc = condresponse(cellMeanTrial[cell,:],condition)  # condtion response, averaged over repeats
            fm,fse,_  = factorresponse(mseuc)  # put condition response into factor space
            p = Any[Tuple(argmax(coalesce.(fm,-Inf)))...]
            push!(fms,fm.*100);push!(fses,fse.*100)   # change to percentage (*100)
            for f in collect(keys(fa))
                fd = findfirst(f.==keys(fa))   # facotr dimention
                push!(ufm[f], fa[f][p[fd]])  # find the maximal level of each factor
            end
        end
        # Plot
        if plot
            @manipulate for cell in 1:cellNum
                heatmap(fms[cell])
            end
            @manipulate for cell in 1:cellNum
                blankResp = cellMeanTrial[cell,condition[end,:i]]  # Blank conditions
                histogram(abs.(blankResp), nbins=10)
            end
        end

        ## Get the responsive cells
        uresponsive=[];umodulative=[]
        cti = reduce(append!,condition[1:end-1, :i],init=Int[])   # Drop blank, only choose stim conditions
        for cell in 1:cellNum
            # cell = 1
            # display(cell)
            condResp = cellMeanTrial[cell,cti]
            push!(umodulative,ismodulative([DataFrame(Y=condResp) ctc[cti,:]], alpha=α, interact=true))  # Check molulativeness within stim conditions
            blankResp = cellMeanTrial[cell,condition[end,:i]]  # Choose blank conditions
            # isresp = []
            # for i in 1:condNum
            #     condResp = cellMeanTrial[cell,condition[i, :i]]
            #     push!(isresp, pvalue(UnequalVarianceTTest(blankResp,condResp))<α)   # Check responsiveness between stim condtions and blank conditions
            # end
            condResp = cellMeanTrial[cell,condition[(condition.Dir .==ufm[:Dir][cell]) .& (condition.SpatialFreq .==ufm[:SpatialFreq][cell]), :i][1]]
            push!(uresponsive, pvalue(UnequalVarianceTTest(blankResp,condResp))<α)   # Check responsiveness between stim condtions and blank conditions
            # push!(uresponsive,any(isresp))
                # plotcondresponse(condResp ctc)
                # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[cell])_CondResponse$i")),[".png",".svg"])
        end
        visResp = uresponsive .| umodulative   # Combine responsivenness and modulativeness as visual responsiveness
        display(["uresponsive:", count(uresponsive)])
        display(["umodulative:", count(umodulative)])
        display(["Responsive cells:", count(visResp)])
        result.visResp = visResp
        result.responsive = uresponsive
        result.modulative = umodulative

        ## Check which cell is significantly tuning by orientation or direction
        oriAUC=[]; dirAUC=[];
        for cell in 1:cellNum
            # cell=1  # for test
                # Get all trial Id of under maximal sf
                # mcti = @where(condition, :SpatialFreq .== ufm[:SpatialFreq][cell])
            mcti = condition[condition.SpatialFreq.==ufm[:SpatialFreq][cell], :]
            blankResp = cellMeanTrial[cell,condition[end,:i]]

            oridist=[];dirdist=[];blkoridist=[];blkdirdist=[];
            for k =1:2
                if k ==1
                    resp=[cellMeanTrial[cell,mcti.i[r][t]] for r in 1:nrow(mcti), t in 1:mcti.n[1]]
                elseif k ==2
                    resp = Array{Float64}(undef, nrow(mcti), mcti[1,:n])
                    sample!(blankResp, resp; replace=true, ordered=false)
                end

                for j = 1:sampnum    # Random sampling sampnum times
                    for i=1:size(resp,1)
                        shuffle!(@view resp[i,:])
                    end
                    resu= [factorresponsestats(mcti[:Dir],resp[:,t],factor=:Dir,isfit=false) for t in 1:mcti.n[1]]
                    orivec = reduce(vcat,[resu[t].om for t in 1:mcti.n[1]])
                    orivecmean = mean(orivec, dims=1)[1]  # final mean vec
                    oridistr = [real(orivec) imag(orivec)] * [real(orivecmean) imag(orivecmean)]'  # Project each vector to the final vector, so now it is 1D distribution

                    # check significance of direction selective
                   oriorth = angle(mean(-orivec, dims=1)[1])  # angel orthogonal to mean ori vector
                   orivecdir = exp(im*oriorth/2)   # dir axis vector (orthogonal to ori vector) in direction space
                   dirvec = reduce(vcat,[resu[t].dm for t in 1:mcti.n[1]])
                   dirdistr = [real(dirvec) imag(dirvec)] * [real(orivecdir) imag(orivecdir)]'

                   if k ==1
                       push!(oridist, oridistr);push!(dirdist, dirdistr);
                   elseif k==2
                       push!(blkoridist, oridistr); push!(blkdirdist, dirdistr);
                   end
                end
            end

            blkoridist = reduce(vcat, blkoridist)
            blkdirdist = reduce(vcat, blkdirdist)
            oridist = reduce(vcat, oridist)
            dirdist = reduce(vcat, dirdist)

            oriauc=roccurve(blkoridist, oridist)
            dirauc=roccurve(blkdirdist, dirdist)

            push!(oriAUC,oriauc);push!(dirAUC,dirauc);

        end
        result.oriauc = oriAUC
        result.dirauc = dirAUC

        ## Get the optimal factor level using Circular Variance for each cell
        ufs = Dict(k=>[] for k in keys(fa))
        for cell in 1:length(fms), f in collect(keys(fa))
            p = Any[Tuple(argmax(coalesce.(fms[cell],-Inf)))...] # Replace missing with -Inf, then find the x-y coordinates of max value.
            fd = findfirst(f.==keys(fa))   # facotr dimention
            fdn = length(fa[f])  # dimention length/number of factor level
            p[fd]=1:fdn   # pick up a slice for plotting tuning curve
            mseuc=DataFrame(m=fms[cell][p...],se=fses[cell][p...],u=fill(cellId[cell],fdn),ug=fill(parse(Int, recordPlane), fdn))  # make DataFrame for plotting
            mseuc[f]=fa[f]

            # The optimal dir, ori (based on circular variance) and sf (based on log2 fitting)
            push!(ufs[f],factorresponsestats(dropmissing(mseuc)[f],dropmissing(mseuc)[:m],factor=f, isfit=oriAUC[cell]>fitThres))
            # plotcondresponse(dropmissing(mseuc),colors=[:black],projection=[],responseline=[], responsetype=:ResponseF)
            # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_$(f)_Tuning$i")),[".png"]#,".svg"])
        end

        tempDF=DataFrame(ufs[:SpatialFreq])
        result.fitsf = tempDF.osf
        tempDF=DataFrame(ufs[:Dir])
        result.cvdir = tempDF.od   # preferred direction from cv
        result.dircv = tempDF.dcv
        result.fitdir =map(i->isempty(i) ? NaN : :pd in keys(i) ? i.pd : NaN,tempDF.fit)  # preferred direction from fitting
        result.dsi =map(i->isempty(i) ? NaN : :dsi1 in keys(i) ? i.dsi1 : NaN,tempDF.fit)
        result.cvori = tempDF.oo  # preferred orientation from cv
        result.oricv = tempDF.ocv
        result.fitori =map(i->isempty(i) ? NaN : :po in keys(i) ? i.po : NaN,tempDF.fit)  # preferred orientation from cv
        result.osi =map(i->isempty(i) ? NaN : :osi1 in keys(i) ? i.osi1 : NaN,tempDF.fit)

        # Plot tuning curve of each factor of each cell
        if plot
            @manipulate for u in 1:length(fms), f in collect(keys(fa))
                p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]  # Replace missing with -Inf, then find the x-y coordinates of max value.
                fd = findfirst(f.==keys(fa))   # facotr dimention
                fdn = length(fa[f])  # dimention length/number of factor level
                p[fd]=1:fdn  # pick up a slice for plotting tuning curve
                mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(cellId[u],fdn),ug=fill(parse(Int, recordPlane), fdn))  # make DataFrame for plotting
                mseuc[f]=fa[f]
                plotcondresponse(dropmissing(mseuc),colors=[:black],projection=:polar,responseline=[], responsetype=:ResponseF)
            end
        end

        #Save results
        CSV.write(joinpath(resultFolder,join([subject,"_",siteId,"_result.csv"])), result)
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_result.jld2"])), "result",result)
    end
end

function process_2P_dirsfcolor(files,param;uuid="",log=nothing,plot=false)

    interpolatedData = haskey(param,:interpolatedData) ? param[:interpolatedData] : true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
    preOffset = haskey(param,:preOffset) ? param[:preOffset] : 0.1
    responseOffset = haskey(param,:responseOffset) ? param[:responseOffset] : 0.05  # in sec
    α = haskey(param,:α) ? param[:α] : 0.05   # p value
    sampnum = haskey(param,:sampnum) ? param[:sampnum] : 100   # random sampling 100 times
    fitThres = haskey(param,:fitThres) ? param[:fitThres] : 0.5
    hueSpace = haskey(param,:hueSpace) ? param[:hueSpace] : "DKL"   # Color space used? DKL or HSL
    diraucThres = haskey(param,:diraucThres) ? param[:diraucThres] : 0.8   # if passed, calculate hue direction, otherwise calculate hue axis
    oriaucThres = haskey(param,:oriaucThres) ? param[:oriaucThres] : 0.5
    # Respthres = haskey(param,:Respthres) ? param[:Respthres] : 0.1  # Set a response threshold to filter out low response cells?
    blankId = haskey(param,:blankId) ? param[:blankId] : 36  # Blank Id
    excId = haskey(param,:excId) ? param[:excId] : [27,28]  # Exclude some condition?
    excId = vcat(excId,blankId)
    # Expt info
    dataset = prepare(files)
    ex = dataset["ex"]
    disk=ex["source"][1:2];subject=ex["Subject_ID"];recordSession=uppercasefirst(ex["RecordSite"]);testId=ex["TestID"]

    ## Prepare data & result path
    exptId = join(filter(!isempty,[recordSession[2:end], testId]),"_")
    dataFolder = joinpath(disk,subject, "2P_data", recordSession, exptId)

    ## load expt, scanning parameters
    envparam = ex["EnvParam"]
    sbx = dataset["sbx"]["info"]

    ## Align Scan frames with stimulus
    # Calculate the scan parameters
    scanFreq = sbx["resfreq"]
    lineNum = sbx["sz"][1]
    if haskey(sbx, "recordsPerBuffer_bi")
       scanMode = 2  # bidirectional scanning   # if process Splitted data, =1
    else
       scanMode = 1  # unidirectional scanning
    end
    sbxfs = 1/(lineNum/scanFreq/scanMode)   # frame rate
    trialOnLine = sbx["line"][1:2:end]
    trialOnFrame = sbx["frame"][1:2:end] + round.(trialOnLine/lineNum)        # if process splitted data use frame_split
    trialOffLine = sbx["line"][2:2:end]
    trialOffFrame = sbx["frame"][2:2:end] + round.(trialOnLine/lineNum)    # if process splitted data use frame_split

    # On/off frame indces of trials
    trialEpoch = Int.(hcat(trialOnFrame, trialOffFrame))
    # minTrialDur = minimum(trialOffFrame-trialOnFrame)
    # histogram(trialOffFrame-trialOnFrame,nbins=20,title="Trial Duration(Set to $minTrialDur)")

    # Transform Trials ==> Condition
    ctc = DataFrame(ex["CondTestCond"])
    trialNum =  size(ctc,1)
    conditionAll = condin(ctc)
    # Remove extra conditions (only for AF4), and seperate blanks
    others = (:ColorID, excId)
    condi = .!in.(conditionAll[!,others[1]],[others[2]])
    # factors = finalfactor(conditionAll)
    conditionCond = conditionAll[condi,:]
    condNum = size(conditionCond,1) # not including blanks

    # Extract blank condition
    blank = (:ColorID, blankId)
    bi = in.(conditionAll[!,blank[1]],[blank[2]])
    conditionBlank = conditionAll[bi,:]
    # replace!(bctc.ColorID, 36 =>Inf)

    # Change ColorID ot HueAngle if needed
    if hueSpace == "DKL"
        ucid = sort(unique(conditionCond.ColorID))
        hstep = 360/length(ucid)
        conditionCond.ColorID = (conditionCond.ColorID.-minimum(ucid)).*hstep
        conditionCond=rename(conditionCond, :ColorID => :HueAngle)
    end


    # On/off frame indces of condations/stimuli
    preStim = ex["PreICI"]; stim = ex["CondDur"]; postStim = ex["SufICI"]
    trialOnTime = fill(0, trialNum)
    condOfftime = preStim + stim
    preEpoch = [0 preStim-preOffset]
    condEpoch = [preStim+responseOffset condOfftime-responseOffset]
    preFrame=epoch2samplerange(preEpoch, sbxfs)
    condFrame=epoch2samplerange(condEpoch, sbxfs)
    # preOn = fill(preFrame.start, trialNum)
    # preOff = fill(preFrame.stop, trialNum)
    # condOn = fill(condFrame.start, trialNum)
    # condOff = fill(condFrame.stop, trialNum)

    ## Load data
    segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,adddir=true)[1]
    segment = prepare(segmentFile)
    signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,adddir=true)[1]
    signal = prepare(signalFile)
    sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
    spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

    planeNum = size(segment["mask"],3)  # how many planes
    # planeNum = 1
    planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)

    ## Use for loop process each plane seperately
    for pn in 1:planeNum
        # pn=1  # for test
        # Initialize DataFrame for saving results
        recordPlane = string("00",pn-1)  # plane/depth, this notation only works for expt has less than 10 planes
        siteId = join(filter(!isempty,[recordSession[2:end], testId, recordPlane]),"_")
        dataExportFolder = joinpath(disk,subject, "2P_analysis", recordSession, siteId, "DataExport")
        resultFolder = joinpath(disk,subject, "2P_analysis", recordSession, siteId, "Plots")
        isdir(dataExportFolder) || mkpath(dataExportFolder)
        isdir(resultFolder) || mkpath(resultFolder)
        result = DataFrame()

        cellRoi = segment["seg_ot"]["vert"][pn]
        cellNum = length(cellRoi)
        display("plane: $pn")
        display("Cell Number: $cellNum")
        cellId = collect(range(1, step=1, stop=cellNum))  # Currently, the cellID is the row index of signal

        if interpolatedData
            rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
            # spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        else
            rawF = transpose(signal["sig_ot"]["sig"][pn])
            # spike = transpose(signal["sig_ot"]["spks"][pn])
        end
        result.py = 0:cellNum-1
        result.ani = fill(subject, cellNum)
        result.dataId = fill(siteId, cellNum)
        result.cellId = 1:cellNum

        ## Plot dF/F traces of all trials for all cells
        # Cut raw fluorescence traces according to trial on/off time and calculate dF/F
        cellTimeTrial = sbxsubrm(rawF,trialEpoch,cellId;fun=dFoF(preFrame))  # cellID x timeCourse x Trials
        # Mean response within stim time
        cellMeanTrial = dropdims(mean(cellTimeTrial[:,condFrame,:], dims=2), dims=2)  # cellID x Trials
        # Plot
        if plot
            @manipulate for cell in 1:cellNum
                plotanalog(transpose(cellTimeTrial[cell,:,:]), fs=sbxfs, timeline=condEpoch.-preStim, xunit=:s, ystep=1,cunit=:p, color=:fire,xext=preStim)
            end
        end

        ## Average over repeats, and put each cell's response in factor space (hue-dir-sf...), and find the maximal level of each factor
        factors = finalfactor(conditionCond)
        fa = OrderedDict(f=>unique(conditionCond[f]) for f in factors)  # factor levels, the last one of each factor maybe blank(Inf)
        fms=[];fses=[];  # factor mean spac, factor sem space, mean and sem of each condition of each cell
        ufm = Dict(k=>[] for k in keys(fa))  # maxi factor level of each cell
        for cell in 1:cellNum
            # cell=1
            mseuc = condresponse(cellMeanTrial[cell,:],conditionCond)  # condtion response, averaged over repeats
            fm,fse,_  = factorresponse(mseuc)  # put condition response into factor space
            p = Any[Tuple(argmax(coalesce.(fm,-Inf)))...]
            push!(fms,fm.*100);push!(fses,fse.*100)   # change to percentage (*100)
            for f in collect(keys(fa))
                fd = findfirst(f.==keys(fa))   # facotr dimention
                push!(ufm[f], fa[f][p[fd]])  # find the maximal level of each factor
            end
        end
        # Plot
        if plot
            # @manipulate for cell in 1:cellNum
            #     heatmap(fms[cell])
            # end
            @manipulate for cell in 1:cellNum
                # blankResp = cellMeanTrial[cell,vcat(conditionBlank[:,:i]...)]  # Blank conditions
                # histogram(abs.(blankResp), nbins=10)
                condResp = cellMeanTrial[cell,vcat(conditionCond[:,:i]...)]  # Stim conditions
                histogram(abs.(condResp), nbins=10)
            end
        end

        ## Get the responsive cells & blank response
        mseub=[];uresponsive=[];umodulative=[];  #uhighResp=[];
        cti = reduce(append!,conditionCond[:, :i],init=Int[])  # Choose hue condition, exclude blanks and others
        for cell in 1:cellNum
            # cell=1
            condResp = cellMeanTrial[cell,cti]  #
            push!(umodulative,ismodulative([DataFrame(Y=condResp) ctc[cti,:]], alpha=α, interact=true))  # Check molulativeness within stim conditions
            blankResp = cellMeanTrial[cell,vcat(conditionBlank[:,:i]...)]  # Choose blank conditions
            mseuc = condresponse(cellMeanTrial[cell,:],[vcat(conditionBlank[:,:i]...)]) # Get the mean & sem of blank response for a cell
            push!(mseub, mseuc)
            # isresp = []
            # for i in 1:condNum
            #     condResp = cellMeanTrial[cell,condition[i, :i]]
            #     push!(isresp, pvalue(UnequalVarianceTTest(blankResp,condResp))<α)
            # end
            # condResp = cellMeanTrial[cell,condition[(condition.Dir .==ufm[:Dir][cell]) .& (condition.SpatialFreq .==ufm[:SpatialFreq][cell]), :i][1]]
            condResp = cellMeanTrial[cell,conditionCond[(conditionCond.HueAngle .==ufm[:HueAngle][cell]).& (conditionCond.Dir .==ufm[:Dir][cell]) .& (conditionCond.SpatialFreq .==ufm[:SpatialFreq][cell]), :i][1]]
            push!(uresponsive, pvalue(UnequalVarianceTTest(blankResp,condResp))<α)   # Check responsiveness between stim condtions and blank conditions
            # push!(uhighResp, mean(condResp,dims=1)[1]>Respthres)
            # plotcondresponse(condResp cctc)
            # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[cell])_CondResponse$i")),[".png",".svg"])
        end

        visResp = uresponsive .| umodulative   # Combine responsivenness and modulativeness as visual responsiveness
        display(["uresponsive:", count(uresponsive)])
        # display(["uhighResp:", count(uhighResp)])
        display(["umodulative:", count(umodulative)])
        display(["Responsive cells:", count(visResp)])
        result.visResp = visResp
        result.responsive = uresponsive
        # result.uhighResp = uhighResp
        result.modulative = umodulative

        ## Check which cell is significantly tuning by orientation or direction
        # oripvalue=[];orivec=[];dirpvalue=[];dirvec=[];
        # hueaxpvalue=[];huedirpvalue=[];opratioc=[];dpratioc=[];
        oriAUC=[]; dirAUC=[]; hueaxAUC=[]; huedirAUC=[];
        for cell in 1:cellNum
            # cell=1  # for test
            # Get all trial Id of under maximal sf
            # mcti = @where(condition, :SpatialFreq .== ufm[:SpatialFreq][cell])
            mcti = conditionCond[(conditionCond.HueAngle.==ufm[:HueAngle][cell]).&(conditionCond.SpatialFreq.==ufm[:SpatialFreq][cell]), :]
            mbti = conditionBlank[conditionBlank.SpatialFreq.==ufm[:SpatialFreq][cell], :]
            blankResp = [cellMeanTrial[cell,conditionBlank.i[r][t]] for r in 1:nrow(conditionBlank), t in 1:conditionBlank.n[1]]
            # resp = [cellMeanTrial[cell,mcti.i[r][t]] for r in 1:nrow(mcti), t in 1:mcti.n[1]]
            # resu= [factorresponsestats(mcti[:dir],resp[:,t],factor=:dir,isfit=false) for t in 1:mcti.n[1]]
            # orivec = reduce(vcat,[resu[t].oov for t in 1:mcti.n[1]])
            # pori=[];pdir=[];pbori=[];pbdir=[];

            oridist=[];dirdist=[];blkoridist=[];blkdirdist=[];
            for k =1:2
                if k ==1
                    resp = [cellMeanTrial[cell,mcti.i[r][t]] for r in 1:nrow(mcti), t in 1:mcti[1,:n]]
                elseif k==2
                    resp = Array{Float64}(undef, nrow(mcti), ncol(mcti))
                    sample!(blankResp, resp; replace=true, ordered=false)
                end

                for j = 1:sampnum    # Random sampling sampnum times
                    for i=1:size(resp,1)
                        shuffle!(@view resp[i,:])
                    end
                    resu= [factorresponsestats(mcti[:Dir],resp[:,t],factor=:Dir, isfit=false) for t in 1:mcti[1,:n]]
                    orivec = reduce(vcat,[resu[t].om for t in 1:mcti[1,:n]])
                    orivecmean = mean(orivec, dims=1)[1]  # final mean vec
                    oridistr = [real(orivec) imag(orivec)] * [real(orivecmean) imag(orivecmean)]'  # Project each vector to the final vector, so now it is 1D distribution

                    # orip = hotellingt2test([real(orivec) imag(orivec)],[0 0],0.05)
                    # push!(orivec, orivectemp)
                    # check significance of direction selective
                    oriorth = angle(mean(-orivec, dims=1)[1])  # angel orthogonal to mean ori vector
                    orivecdir = exp(im*oriorth/2)   # dir axis vector (orthogonal to ori vector) in direction space
                    dirvec = reduce(vcat,[resu[t].dm for t in 1:mcti[1,:n]])
                    dirdistr = [real(dirvec) imag(dirvec)] * [real(orivecdir) imag(orivecdir)]'
                #    dirp = dirsigtest(orivecdir, dirvec)
                   if k ==1
                    #    push!(pori, orip);push!(pdir, dirp);
                    push!(oridist, oridistr);push!(dirdist, dirdistr);
                   elseif k==2
                    #    push!(pbori, orip);push!(pbdir, dirp);
                    push!(blkoridist, oridistr); push!(blkdirdist, dirdistr);
                   end
                end
            end
            # yso,nso,wso,iso=histrv(float.(pori),0,1,nbins=20)
            # ysd,nsd,wsd,isd=histrv(float.(pdir),0,1,nbins=20)
            # opfi=mean(yso[findmax(nso)[2]])
            # dpfi=mean(ysd[findmax(nsd)[2]])
            # ysob,nsob,wsob,isob=histrv(float.(pbori),0,1,nbins=20)
            # ysdb,nsdb,wsdb,isdb=histrv(float.(pbdir),0,1,nbins=20)
            # opratio=sum(nsob[map(i-> mean(wsob[i])<opfi, range(1,size(wsob,1)))])/sum(nsob)
            # dpratio=sum(nsdb[map(i-> mean(wsdb[i])<dpfi, range(1,size(wsdb,1)))])/sum(nsdb)
            # push!(oripvalue,opfi); push!(dirpvalue,dpfi);
            # push!(opratioc,opratio); push!(dpratioc,dpratio);
            blkoridist = reduce(vcat, blkoridist)
            blkdirdist = reduce(vcat, blkdirdist)
            oridist = reduce(vcat, oridist)
            dirdist = reduce(vcat, dirdist)

            oriauc=roccurve(blkoridist, oridist)
            dirauc=roccurve(blkdirdist, dirdist)

            push!(oriAUC,oriauc);push!(dirAUC,dirauc);


            if dirauc >= diraucThres #dirpvalue[cell] < α  # direction selective
                mcti = conditionCond[(conditionCond.Dir.==ufm[:Dir][cell]).&(conditionCond.SpatialFreq.==ufm[:SpatialFreq][cell]), :]
            elseif (dirauc < diraucThres) .& (oriauc > oriaucThres) #(dirpvalue[cell] > α) .& (oripvalue[cell] < α)   # no direction selective, but has orientation selective
                mcti = conditionCond[((conditionCond.Dir.==ufm[:Dir][cell]) .| (conditionCond.Dir.==mod.(ufm[:Dir][cell]+180,360))).&(conditionCond.SpatialFreq.==ufm[:SpatialFreq][cell]), :]
                mcti = by(mcti, :HueAngle, n=:n=>sum, i=:i=>d->[reduce(hcat,d')])
            else  # neither direction nor orientation selective
                mcti = conditionCond[conditionCond.SpatialFreq.==ufm[:SpatialFreq][cell], :]
                mcti = by(mcti, :HueAngle, n=:n=>sum, i=:i=>d->[reduce(hcat,d')])
            end

            # phueax=[];phuedir=[];pbdir=[];
            hueaxdist=[];huedirdist=[];blkhueaxdist=[];blkhuedirdist=[];
            for k =1:2
                if k ==1  # hue condition
                    resp = [cellMeanTrial[cell,mcti.i[r][t]] for r in 1:nrow(mcti), t in 1:mcti[1,:n]]
                elseif k==2  # blank condition
                    resp = Array{Float64}(undef, nrow(mcti), mcti[1,:n])
                    sample!(blankResp, resp; replace=true, ordered=false)
                end

                for j = 1:sampnum    # Random sampling sampnum times
                    for i=1:size(resp,1)
                        shuffle!(@view resp[i,:])
                    end

                    resu= [factorresponsestats(mcti[:HueAngle],resp[:,t],factor=:HueAngle,isfit=false) for t in 1:mcti[1,:n]]
                    huevec = reduce(vcat,[resu[t].ham for t in 1:mcti.n[1]])  # hue axis
                    # hueaxp = hotellingt2test([real(huevec) imag(huevec)],[0 0],0.05)
                    huevecmean = mean(huevec, dims=1)[1]  # final mean vec
                    hueaxdistr = [real(huevec) imag(huevec)] * [real(huevecmean) imag(huevecmean)]'  # Project each vector to the final vector, so now it is 1D distribution

                    huevec = reduce(vcat,[resu[t].hm for t in 1:mcti.n[1]])  # hue direction
                    # huedirp = hotellingt2test([real(huevec) imag(huevec)],[0 0],0.05)
                    huevecmean = mean(huevec, dims=1)[1]  # final mean vec
                    huedirdistr = [real(huevec) imag(huevec)] * [real(huevecmean) imag(huevecmean)]'  # Project each vector to the final vector, so now it is 1D distribution
                    # push!(phueax,hueaxp)
                    # push!(phuedir,huedirp)
                    if k ==1
                        push!(hueaxdist, hueaxdistr);push!(huedirdist, huedirdistr);
                    elseif k==2
                        push!(blkhueaxdist, hueaxdistr);push!(blkhuedirdist, huedirdistr);
                    end
                end
            end
            # ysh,nsh,wsh,ish=histrv(float.(phueax),0,1,nbins=20)
            # push!(hueaxpvalue,mean(ysh[findmax(nsh)[2]]))
            # ysh,nsh,wsh,ish=histrv(float.(phuedir),0,1,nbins=20)
            # push!(huedirpvalue,mean(ysh[findmax(nsh)[2]]))
            blkhueaxdist = reduce(vcat, blkhueaxdist)
            blkhuedirdist = reduce(vcat, blkhuedirdist)
            hueaxdist = reduce(vcat, hueaxdist)
            huedirdist = reduce(vcat, huedirdist)

            hueaxauc=roccurve(blkhueaxdist, hueaxdist)
            huedirauc=roccurve(blkhuedirdist, huedirdist)

            push!(hueaxAUC,hueaxauc);push!(huedirAUC,huedirauc);
        end
        # result.orip = oripvalue
        # result.opratio=opratioc
        # result.dirp = dirpvalue
        # result.dpratio=dpratioc
        # result.hueaxp = hueaxpvalue
        # result.huedirp = huedirpvalue
        result.hueoriauc = oriAUC
        result.huedirauc = dirAUC
        result.hueaxauc = hueaxAUC
        result.huediauc = huedirAUC

        ## Get the optimal factor level using Circular Variance for each cell
        ufs = Dict(k=>[] for k in keys(fa))
        for cell in 1:length(fms), f in collect(keys(fa))
            p = Any[Tuple(argmax(coalesce.(fms[cell],-Inf)))...] # Replace missing with -Inf, then find the x-y coordinates of max value.
            fd = findfirst(f.==keys(fa))   # facotr dimention
            fdn = length(fa[f])  # dimention length/number of factor level
            p[fd]=1:fdn   # pick up a slice for plotting tuning curve
            mseuc=DataFrame(m=fms[cell][p...],se=fses[cell][p...],u=fill(cellId[cell],fdn),ug=fill(parse(Int, recordPlane), fdn))  # make DataFrame for plotting
            mseuc[f]=fa[f]

            # The optimal dir, ori (based on circular variance) and sf (based on log2 fitting)
            push!(ufs[f],factorresponsestats(dropmissing(mseuc)[f],dropmissing(mseuc)[:m],factor=f, isfit=max(oriAUC[cell], dirAUC[cell], hueaxAUC[cell], huedirAUC[cell])>fitThres))
            # plotcondresponse(dropmissing(mseuc),colors=[:black],projection=[],responseline=[], responsetype=:ResponseF)
            # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_$(f)_Tuning$i")),[".png"]#,".svg"])
        end

        tempDF=DataFrame(ufs[:SpatialFreq])
        result.fithuesf = tempDF.osf
        tempDF=DataFrame(ufs[:Dir])
        result.cvhuedir = tempDF.od   # cv
        result.huedircv = tempDF.dcv
        result.fithuedir =map(i->isempty(i) ? NaN : :pd in keys(i) ? i.pd : NaN,tempDF.fit)  # fitting
        result.huedsi =map(i->isempty(i) ? NaN : :dsi1 in keys(i) ? i.dsi1 : NaN,tempDF.fit)
        result.cvhueori = tempDF.oo  # cv
        result.hueoricv = tempDF.ocv
        result.fithueori =map(i->isempty(i) ? NaN : :po in keys(i) ? i.po : NaN,tempDF.fit)  # fitting
        result.hueosi =map(i->isempty(i) ? NaN : :osi1 in keys(i) ? i.osi1 : NaN,tempDF.fit)
        tempDF=DataFrame(ufs[:HueAngle])
        result.cvhueax = tempDF.oha # cv
        result.hueaxcv = tempDF.hacv
        result.fithueax =map(i->isempty(i) ? NaN : :pha in keys(i) ? i.pha : NaN,tempDF.fit)  # fitting
        result.hueaxsi =map(i->isempty(i) ? NaN : :hasi1 in keys(i) ? i.hasi1 : NaN,tempDF.fit)
        result.cvhuedi = tempDF.oh # cv
        result.huedicv = tempDF.hcv
        result.fithuedi =map(i->isempty(i) ? NaN : :ph in keys(i) ? i.ph : NaN,tempDF.fit)  # fitting
        result.huedisi =map(i->isempty(i) ? NaN : :hsi1 in keys(i) ? i.hsi1 : NaN,tempDF.fit)
        result.maxhue = tempDF.maxh
        result.maxhueresp = tempDF.maxr

        # Plot tuning curve of each factor of each cell
        # plot = true
        if plot
            @manipulate for u in 1:length(fms), f in collect(keys(fa))
                p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]  # Replace missing with -Inf, then find the x-y coordinates of max value.
                fd = findfirst(f.==keys(fa))   # facotr dimention
                fdn = length(fa[f])  # dimention length/number of factor level
                p[fd]=1:fdn  # pick up a slice for plotting tuning curve
                mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(cellId[u],fdn),ug=fill(parse(Int, recordPlane), fdn))  # make DataFrame for plotting
                mseuc[f]=fa[f]
                plotcondresponse(dropmissing(mseuc),colors=[:black],projection=:polar,responseline=[], responsetype=:ResponseF)
            end
        end

        #Save results
        CSV.write(joinpath(resultFolder,join([subject,"_",siteId,"_result.csv"])), result)
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_result.jld2"])), "result",result)
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_tuning.jld2"])), "tuning",tempDF)
    end

end

function process_2P_hartleyFourier(files,param;uuid="",log=nothing,plot=false)

    interpolatedData = haskey(param,:interpolatedData) ? param[:interpolatedData] : true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
    delays = 0:0.05:0.5
    ntau = length(collect(delays))
    print(collect(delays))

    # Expt info
    dataset = prepare(files)
    ex = dataset["ex"]
    disk=ex["source"][1:2];subject=ex["Subject_ID"];recordSession=uppercasefirst(ex["RecordSite"]);testId=ex["TestID"]

    ## Prepare data & result path
    exptId = join(filter(!isempty,[recordSession[2:end], testId]),"_")
    dataFolder = joinpath(disk,subject, "2P_data", recordSession, exptId)

    ## load expt, scanning parameters
    envparam = ex["EnvParam"]
    coneType = getparam(envparam,"colorspace")
    szhtly_visangle = envparam["x_size"]  # deg
    sbx = dataset["sbx"]["info"]
    sbxft = ex["frameTimeSer"]   # time series of sbx frame in whole recording
    # Condition Tests
    envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    nstim = size(condidx,1)
    # condtable = DataFrame(ex["Cond"])
    condtable = DataFrame(ex["raw"]["log"]["randlog_T1"]["domains"]["Cond"])
    rename!(condtable, [:oridom, :kx, :ky,:bwdom,:colordom])
    condtable[:kx] = [Int(x) for x in condtable[:kx]]
    condtable[:ky] = [Int(x) for x in condtable[:ky]]
    max_k = max(abs.(condtable.kx)...)
    # find out blanks and unique conditions
    blkidx = condidx.>5641  # blanks start from 5641
    cidx = .!blkidx
    condidx2 = condidx.*cidx + blkidx.* 5641
    conduniq = unique(condidx2)
    ## Load data
    segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,adddir=true)[1]
    segment = prepare(segmentFile)
    signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,adddir=true)[1]
    signal = prepare(signalFile)
    # sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
    spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

    ## Load data
    planeNum = size(segment["mask"],3)  # how many planes
    planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)

    ## Use for loop process each plane seperately
    for pn in 1:planeNum
        # pn=2  # for test
        # Initialize DataFrame for saving results
        recordPlane = string("00",pn-1)  # plane/depth, this notation only works for expt has less than 10 planes
        siteId = join(filter(!isempty,[recordSession[2:end], testId, recordPlane]),"_")
        dataExportFolder = joinpath(disk,subject, "2P_analysis", recordSession, siteId, "DataExport")
        resultFolder = joinpath(disk,subject, "2P_analysis", recordSession, siteId, "Plots")
        isdir(dataExportFolder) || mkpath(dataExportFolder)
        isdir(resultFolder) || mkpath(resultFolder)
        result = DataFrame()

        cellRoi = segment["seg_ot"]["vert"][pn]
        cellNum = length(cellRoi)
        display("plane: $pn")
        display("Cell Number: $cellNum")

        if interpolatedData
            # rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
            spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        else
            # rawF = transpose(signal["sig_ot"]["sig"][pn])
            spike = transpose(signal["sig_ot"]["spks"][pn])
        end
        result.py = 0:cellNum-1
        result.ani = fill(subject, cellNum)
        result.dataId = fill(siteId, cellNum)
        result.cellId = 1:cellNum
        ## Chop spk trains according delays
        spk=zeros(nstim,ntau,cellNum)
        for d in eachindex(delays)
            y,num,wind,idx = subrv(sbxft,condon.+delays[d], condoff.+delays[d],isminzero=false,ismaxzero=false,shift=0,israte=false)
            for i =1:nstim
                spkepo = @view spike[:,idx[i][1]:idx[i][end]]
                spk[i,d,:]= mean(spkepo, dims=2)
            end
        end

        ## Sum cell response of different repeats
        r = zeros(2*max_k+1,2*max_k+1,ntau,cellNum)   # Blank condition [0,0] is now at [max_k+1, max_k+1]
        for i=1:nstim
            r[-condtable.kx[condidx[i]]+1+max_k,condtable.ky[condidx[i]]+1+max_k,:,:] = r[-condtable.kx[condidx[i]]+1+max_k,condtable.ky[condidx[i]]+1+max_k,:,:]+spk[i,:,:]
        end

        #  Normalize by stim repeats.
        reps = zeros(size(conduniq,1))
        for i in 1:size(conduniq,1)
             rep = length(findall(x->x==conduniq[i],condidx2))
             reps[i] = rep
            r[-condtable.kx[conduniq[i]]+1+max_k,condtable.ky[conduniq[i]]+1+max_k,:,:] ./= rep
        end

        ## Filter 2D tuning map
        for t = 1:ntau
            for n = 1:cellNum
                rf = r[:,:,t,n]
                rf = rf + rot180(rf) # average over phases PL
                r[:,:,t,n] = imfilter(rf,Kernel.gaussian((1,1),(3,3)))
            end
        end

        ## PL: Build a complax plane with the same size as Hartley space (-kxmax:kxmax, -kymax:kymax) for sf and ori estimation
        szhtly = 2*max_k+1
        vect = collect(-max_k:max_k)
        xx = repeat(vect',szhtly,1)
        yy = repeat(reverse(vect),1,szhtly)
        zz= xx + 1im .* yy
        # heatmap(angle.(zz),yflip=true)  # -pi to pi

        ## find best kernel and estimate preferred sf and ori

        taumax=[];kur=[];kurmax=[];kernraw=[];kernnor=[];
        signif=[];sfest=[];oriest=[];
        for i = 1:cellNum
            # i=15
            z = r[:,:,:,i]
            q = reshape(z,szhtly^2,:)   # in this case, there are 61^2 pixels in the stimulus.
            k = dropdims(mapslices(kurtosis,q;dims=1).-3, dims=1) # The kurtosis of any univariate normal distribution is 3. It is common to compare the kurtosis of a distribution to this value.
            tmax = findall(x->x==max(k...),k)[1]
            kmax = max(k...)
            sig = kmax>7
            kernRaw = z[:,:,tmax]  # raw kernel without
            kern = log10.(z[:,:,tmax] ./ z[max_k+1,max_k+1,tmax])

            # estimate ori/sf
            # bw = kern .> (max(kern...) .* 0.95)
            bw = kern .== max(kern...)
            idx = findall(x->x==1,bw)
            foreach(x->if x[1]>31 bw[x]=0 end,idx)
            zzm = sum(sum(zz.*bw,dims=1),dims=2)[1] / (length(idx)/2)

            sf = abs(zzm)/szhtly_visangle  # cyc/deg
            ori = rad2deg(angle(zzm)) # deg

            push!(taumax,tmax);push!(kur,k);push!(kurmax, kmax); push!(kernraw,kernRaw);push!(kernnor,kern);
            push!(signif,sig);push!(oriest,ori); push!(sfest,sf);

            # if sig == true
            #     heatmap(kern,yflip=true, aspect_ratio=:equal,color=:coolwarm)
            #     # plot([0 real(zzm)],[0 imag(zzm)],'wo-','linewidth',3,'markersize',14);
            # end
        end
        result.sig = signif
        result.oriest = oriest
        result.sfest = sfest
        result.taumax = taumax
        result.kernnor = kernnor
        result.kernraw = kernraw
        result.kurmax=kurmax
        result.kur = kur

        #Save results
        CSV.write(joinpath(resultFolder,join([subject,"_",siteId,"_",coneType,"_tuning_result.csv"])), result)
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_",coneType,"_tuning_result.jld2"])), "result",result)

    end


end
