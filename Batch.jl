using ProgressMeter,Logging,DataFrames,Plots,Images,ePPR
import Base: close

function batchtests(tests::DataFrame,dataroot,resultroot,datatype...;env::Dict{String,Any}=Dict{String,Any}(),log::Dict{String,AbstractLogger}=Dict{String,AbstractLogger}(),delay=20,binwidth=10,isplot=true)
    udf=[];cdf=[]
    p = ProgressMeter.Progress(nrow(tests),1,"Batch Tests ... ",50)
    for t in eachrow(tests)
        try
            u,c=processtest(prepare(abspath(dataroot,t[:files]),datatype...),resultroot,uuid=t[:UUID],env=env,log=log,delay=delay,binwidth=binwidth,isplot=isplot)
            u!=nothing && push!(udf,u)
            c!=nothing && push!(cdf,c)
        catch exc
            display("============================================")
            display("Error In Processing: $(t[:files])")
            display("============================================")
            display.(stacktrace(catch_backtrace()))
        end
        next!(p)
    end
    close(log)
    return vcat(udf...),vcat(cdf...)
end

function processtest(dataset::Dict,resultroot;uuid="",env::Dict{String,Any}=Dict{String,Any}(),log::Dict{String,AbstractLogger}=Dict{String,AbstractLogger}(),delay=20,binwidth=10,isplot=true)
    if haskey(dataset,"ex")
        testid = dataset["ex"]["ID"]
        if testid=="OriGrating"
            processori(dataset,resultroot,uuid=uuid,log=log,delay=delay,binwidth=binwidth,isplot=isplot)
        elseif testid=="Laser"
            processlaser(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        elseif testid=="Image"
            processimage(dataset,resultroot,env,uuid=uuid,log=log,delay=delay,binwidth=binwidth,isplot=isplot)
        elseif testid=="LaserImage"
            processlaserimage(dataset,condroot,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        end
    end
end

function close(log::Dict{String,AbstractLogger})
    for l in values(log)
        close(l.stream)
    end
end
function (log::Dict{String,AbstractLogger})(key,f,logfile)
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
