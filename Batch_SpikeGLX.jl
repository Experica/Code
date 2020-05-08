
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
    minconddur=minimum(condoff-condon)

    # Prepare Conditions
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)

    # Prepare LFP
    lffile = matchfile(Regex("^$(testid)[A-Za-z0-9_]*.imec.lf.bin"),dir = datadir,join=true)[1]
    nsavech=dataset["lf"]["meta"]["nSavedChans"]
    nsample=dataset["lf"]["meta"]["nFileSamp"]
    fs=dataset["lf"]["meta"]["fs"]
    nch = dataset["lf"]["meta"]["snsApLfSy"][2]
    hx,hy,hz = dataset["lf"]["meta"]["probespacing"]
    exchmask = exchmasknp(dataset,type="lf")
    pnrow,pncol = size(exchmask)
    mmlf = Mmap.mmap(lffile,Matrix{Int16},(nsavech,nsample),0)

    # Depth LFP and CSD
    epochdur = timetounit(150)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    # All LFP epochs, gain corrected(voltage), line noise(60,120,180Hz) removed, bandpass filtered, in the same shape of probe, where excluded channels are replaced with local average
    ys=reshape2mask(subrm(mmlf,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"],bandpass=[1,100]),exchmask)

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
    baseindex = epoch2sample([0 basedur],fs)
    dcsd=stfilter(pcsd,temporaltype=:sub,ti=baseindex)
    mdcsd = dropdims(mean(dcsd,dims=3),dims=3)
    depths = hy*(0:(pnrow-1))
    x = (1/fs/SecondPerUnit)*(0:(size(mdcsd,2)-1))
    if plot
        plotanalog(imfilter(mdcsd,Kernel.gaussian((1,1))),fs=fs,y=depths,color=:RdBu)
        foreach(i->savefig(joinpath(resultdir,"Columns_MeandCSD$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"csd.jld2"),"csd",mdcsd,"time",x,"depth",depths,"fs",fs,"log",ex["Log"],"color","$(exparam["ColorSpace"])_$(exparam["Color"])")

    # Depth Power Spectrum
    epochdur = timetounit(300)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    ys=reshape2mask(subrm(mmlf,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"],bandpass=[1,100]),exchmask)
    pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
    ps,freq = powerspectrum(pys,fs,freqrange=[1,100],nw=4*epochdur*SecondPerUnit)

    epoch = [-epochdur 0]
    epochs = condon.+epoch
    ys=reshape2mask(subrm(mmlf,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"],bandpass=[1,100]),exchmask)
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

    # Unit Spike Train
    if plot
        epochext = preicidur
        for u in eachindex(unitspike)
            ys,ns,ws,is = subst(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
            plotspiketrain(ys,timeline=[0,conddur],title="$(ugs[u])Unit_$(unitid[u])")
            foreach(i->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_SpikeTrain$i")),[".png",".svg"])
        end
    end

    # Unit Depth PSTH
    epochdur = timetounit(150)
    epoch = [0 epochdur]
    epochs = condon.+epoch
    bw = timetounit(2)
    psthbins = epoch[1]:bw:epoch[2]

    baseindex = epoch2sample([0 basedur],1/(bw*SecondPerUnit))
    normfun = x->x.-mean(x[baseindex])
    unitpsth = map(st->psthsts(subst(st,epochs,isminzero=true).y,psthbins,israte=true,normfun=normfun),unitspike[unitgood])
    depthpsth,x,depths,depthnunit = spacepsth(unitpsth,unitposition[unitgood,:],hy*(0:pnrow).-(hy/2))
    if plot
        plotpsth(imfilter(depthpsth,Kernel.gaussian((1,1))),x,depths,color=:minmax,n=depthnunit)
        foreach(i->savefig(joinpath(resultdir,"DepthPSTH$i")),[".png",".svg"])
    end
    save(joinpath(resultdir,"depthpsth.jld2"),"depthpsth",depthpsth,"depth",depths,"time",x,"n",depthnunit)

    # Single Unit Binary Spike Train of Condition Tests
    bepochext = timetounit(-300)
    bepoch = [-bepochext minconddur]
    bepochs = condon.+bepoch
    spikebins = bepoch[1]:timetounit(1):bepoch[2] # 1ms bin histogram will convert spike times to binary spike train
    bst = map(st->permutedims(psthsts(subst(st,bepochs,isminzero=true,shift=bepochext).y,spikebins,israte=false,ismean=false).mat),unitspike[unitgood]) # nSpikeBin x nEpoch for each unit
    # Single Unit Correlogram and Circuit
    lag=50
    ccgs,x,ccgis,projs,eunits,iunits,projweights = circuitestimate(bst,lag=lag,unitid=unitid[unitgood],condis=cond.i)

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

function process_hartley_spikeglx(files,param;uuid="",log=nothing,plot=true)
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

function process_condtest_spikeglx(files,param;uuid="",log=nothing,plot=true)
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
