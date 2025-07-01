using Images,HypothesisTests

function process_cycle_imager(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    siteresultdir = joinpath(param[:resultroot],subject,siteid)
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(siteresultdir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # Prepare Conditions
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    modulatefreq = envparam["ModulateTemporalFreq"]
    modulatetype = envparam["ModulateType"]
    modulateduty = envparam["ModulateDuty"]
    ncycle = ex["CondRepeat"]
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    exenv=Dict{Any,Any}("eye"=>ex["Eye"])

    # Prepare Frame
    imagefile = dataset["imagefile"]
    # ds = "H:\\ImagerData"
    # imagefile = map(f->ds * splitdrive(f)[2],imagefile)
    nepoch = dataset["imagenepoch"]
    nframe = length(imagefile)
    framerate = dataset["meta"]["AcquisitionControl"]["AcquisitionFrameRate"]
    exposure = dataset["meta"]["AcquisitionControl"]["ExposureTime"]/1000 # ms
    fmt = dataset["meta"]["DataFormat"]
    w = dataset["meta"]["ImageFormat"]["Width"]
    h = dataset["meta"]["ImageFormat"]["Height"]
    pixfmt = dataset["meta"]["ImageFormat"]["PixelFormat"]
    pixmax = 2^12 - 1

    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : 500 # hemodynamic delay
    baseframeindex = epoch2sampleindex([0 preicidur],framerate,maxsampleindex=nframe)
    frameindex = epoch2sampleindex([0 1000ncycle/modulatefreq].+(preicidur+responsedelay),framerate,maxsampleindex=nframe)
    # For limited memory and analysis purpose, Discrete Fourier Transform are directly evaluated, at few frequencies of interest.
    baseresponse = dropdims(mean(readrawim_Mono12Packed(imagefile[baseframeindex],w,h),dims=3),dims=3) # Float64
    freqs = (f1=modulatefreq,f2=2modulatefreq)
    Fs = dft_imager(imagefile[frameindex],w,h,framerate,baseresponse,freqs...)

    F1phase = angle.(Fs[1]) # [-π, π]
    F1mag = abs.(Fs[1])
    F2phase = angle.(Fs[2]) # [-π, π]
    F2mag = abs.(Fs[2])
    F1phase01 = mod2pi.(F1phase) ./ (2π)
    F2phase01 = mod2pi.(F2phase) ./ (2π)
    F1mag01 = clampscale(F1mag,3)
    F2mag01 = clampscale(F2mag,3)

    if plot
        foreach(ext->save(joinpath(resultdir,"F1Phase$ext"),F1phase01),figfmt)
        foreach(ext->save(joinpath(resultdir,"F1Mag$ext"),F1mag01),figfmt)
        foreach(ext->save(joinpath(resultdir,"F2Phase$ext"),F2phase01),figfmt)
        foreach(ext->save(joinpath(resultdir,"F2Mag$ext"),F2mag01),figfmt)
    end

    if ex["ID"] == "ISICycleOri"
        exenv["cycdir"] = exparam["CycleDirection"]
        if sign(exenv["cycdir"]) < 0
            F1phase01 = 1 .- F1phase01
            F2phase01 = 1 .- F2phase01
        end

        if plot
            diranglemap = map(a->HSV(360a,1,1),F1phase01)
            dirpolarmap = map((a,m)->HSV(360a,1,m),F1phase01,F1mag01)
            orianglemap = map(a->HSV(360a,1,1),F2phase01)
            oripolarmap = map((a,m)->HSV(360a,1,m),F2phase01,F2mag01)

            foreach(ext->save(joinpath(resultdir,"dir_anglemap$ext"),diranglemap),figfmt)
            foreach(ext->save(joinpath(resultdir,"dir_polarmap$ext"),dirpolarmap),figfmt)
            foreach(ext->save(joinpath(resultdir,"ori_anglemap$ext"),orianglemap),figfmt)
            foreach(ext->save(joinpath(resultdir,"ori_polarmap$ext"),oripolarmap),figfmt)
        end
        jldsave(joinpath(resultdir,"isi.jld2");freqs,Fs,F1phase01,F2phase01,exenv,siteid)
    elseif ex["ID"] == "ISICycleColorPlane"
        exenv["cycdir"] = exparam["CycleDirection"]
        if sign(exenv["cycdir"]) < 0
            F1phase01 = 1 .- F1phase01
            F2phase01 = 1 .- F2phase01
        end
        exenv["plane"] = exparam["ModulateParam"]
        exenv["intercept"] = exparam["Intercept"]

        if plot
            if exenv["plane"] == "DKLIsoLum"
                cm = cgrad(ColorMaps["dkl_mcchue_l0"].colors)
                F1phase01_color = map(i->cm[i],F1phase01)
            end
            foreach(ext->save(joinpath(resultdir,"F1phase_color$ext"),F1phase01_color),figfmt)
        end
        jldsave(joinpath(resultdir,"isi.jld2");freqs,Fs,F1phase01,F2phase01,exenv,siteid)
    elseif ex["ID"] == "ISICycle2Color"
        # dft cos phase(max:0->min:π->max:2π) to normalized polarity(min:0->max:1)
        pha2pol(p,s=0) = abs.(mod2pi.(p .+ s) .- π) ./ π
        F1polarity = pha2pol(F1phase)
        F1polarity_dog = clampscale(dogfilter(F1polarity),2)
        F1polarity_ahe = ahe(F1polarity)

        F2polarity = pha2pol(F2phase)
        F2polarity_dog = clampscale(dogfilter(F2polarity),2)
        F2polarity_ahe = ahe(F2polarity)
        
        exenv["color"] = "$(exparam["ColorSpace"])_$(exparam["Color"])"
        mincolor = RGBA(envparam["MinColor"]...)
        maxcolor = RGBA(envparam["MaxColor"]...)
        exenv["minmaxcolor"] = (;mincolor,maxcolor)
        colordist = maxcolor-mincolor
        # if ISI response to cos min/max color modulation is a same phase cos signal, then polarity(0/1) maps to color(min/max)
        if modulatetype == "Sinusoidal"
            # Cos modulation between MinColor and MaxColor
            cmfun = (a;f=1,p=0) -> mincolor + colordist * ((cos(2π*(f*a+p))+1)/2)
        elseif modulatetype == "Square"
            # Cos Sign modulation between MinColor and MaxColor
            cmfun = (a;f=1,p=0) -> mincolor + colordist * ((sign(modulateduty-mod(f*a+p+0.25,1))+1)/2)
        end
        F1polarity_dog_color = map(a->HSV(cmfun(a/2,p=0.5)),F1polarity_dog)
        F1polarity_dog_colormag = map((a,m)->HSV(a.h,a.s,m),F1polarity_dog_color,F1mag01)

        if plot
            foreach(ext->save(joinpath(resultdir,"F1Polarity$ext"),F1polarity),figfmt)
            foreach(ext->save(joinpath(resultdir,"F1Polarity_dog$ext"),F1polarity_dog),figfmt)
            foreach(ext->save(joinpath(resultdir,"F1Polarity_ahe$ext"),F1polarity_ahe),figfmt)
            foreach(ext->save(joinpath(resultdir,"F2Polarity$ext"),F2polarity),figfmt)
            foreach(ext->save(joinpath(resultdir,"F2Polarity_dog$ext"),F2polarity_dog),figfmt)
            foreach(ext->save(joinpath(resultdir,"F2Polarity_ahe$ext"),F2polarity_ahe),figfmt)
            foreach(ext->save(joinpath(resultdir,"F1Polarity_dog_color$ext"),F1polarity_dog_color),figfmt)
            foreach(ext->save(joinpath(resultdir,"F1Polarity_dog_colormag$ext"),F1polarity_dog_colormag),figfmt)
        end
        jldsave(joinpath(resultdir,"isi.jld2");freqs,Fs,F1polarity,F2polarity,exenv,siteid)
    end

end

function process_epoch_imager(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    siteresultdir = joinpath(param[:resultroot],subject,siteid)
    testid = splitdir(splitext(ex["source"])[1])[2]
    resultdir = joinpath(siteresultdir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # Prepare Conditions
    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    isbalance = allequal(cond.n)
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    exenv=Dict{Any,Any}("eye"=>ex["Eye"])

    # Prepare Frame
    imagefile = dataset["imagefile"]
    nepoch = dataset["imagenepoch"]
    minframe = mapreduce(length,min,imagefile)
    framerate = dataset["meta"]["AcquisitionControl"]["AcquisitionFrameRate"]
    exposure = dataset["meta"]["AcquisitionControl"]["ExposureTime"]/1000 # ms
    fmt = dataset["meta"]["DataFormat"]
    w = dataset["meta"]["ImageFormat"]["Width"]
    h = dataset["meta"]["ImageFormat"]["Height"]
    pixfmt = dataset["meta"]["ImageFormat"]["PixelFormat"]
    pixmax = 2^12 - 1

    # Epoch Response
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : conddur/2 # hemodynamic delay(usually 500ms), here use the late near saturated responses
    baseframeindex = epoch2sampleindex([0 preicidur],framerate,maxsampleindex=minframe)
    frameindex = epoch2sampleindex([0 conddur].+(preicidur+responsedelay),framerate,maxsampleindex=minframe)
    epochresponse = Array{Float64}(undef,h,w,nepoch)
    p = ProgressMeter.Progress(nepoch,desc="Epoch Response ... ")
    # ds = "H:\\ImagerData"
    # ds = "D:\\"
    Threads.@threads for i in 1:nepoch
        # imgfile = map(f->ds * splitdrive(f)[2],imagefile[i])
        # epochresponse[:,:,i] = frameresponse_imager(imgfile,w,h,frameindex,baseframeindex)
        epochresponse[:,:,i] = frameresponse_imager(imagefile[i],w,h,frameindex,baseframeindex)
        next!(p)
    end

    if ex["ID"] in ["ISIEpochOri8","ISIEpochOri12"]
        ds = (cond.Ori.+90).%360
        os = cond.Ori.%180
        qs = cond.Ori.%90
        uo = unique(os)
        uq = unique(qs)

        dp = map(o->mod.([o,o+180].+90,360),uo)
        dpi = map(p->map(d->cond.i[d.==ds],p),dp)
        dpt = map(is->pairtest(epochresponse,vcat(is[1]...),vcat(is[2]...)).stat,dpi)
        dcmap,amap,mmap = complexmap([dpt;dpt],deg2rad.([uo;uo.+180].+90),rsign=repeat([-1,1],inner=length(dpt)))
        diranglemap = map(a->HSV(rad2deg(a),1,1),amap)
        dirpolarmap = map((a,m)->HSV(rad2deg(a),1,m),amap,mmap)

        op = map(q->[q,q+90],uq)
        opi = map(p->map(o->cond.i[o.==os],p),op)
        opt = map(is->pairtest(epochresponse,vcat(is[1]...),vcat(is[2]...)).stat,opi)
        ocmap,amap,mmap = complexmap([opt;opt],deg2rad.([uq;uq.+90]),n=2,rsign=repeat([-1,1],inner=length(opt)))
        orianglemap = map(a->HSV(rad2deg(2a),1,1),amap)
        oripolarmap = map((a,m)->HSV(rad2deg(2a),1,m),amap,mmap)

        # rs = condresponse(epochresponse,cond.i).m
        # dcmap,amap,mmap = complexmap(rs,deg2rad.(ds))
        # diranglemap = map(a->HSV(rad2deg(a),1,1),amap)
        # dirpolarmap = map((a,m)->HSV(rad2deg(a),1,m),amap,mmap)

        # or = @views map(o->dropdims(mean(rs[:,:,o.==os],dims=3),dims=3),uo)
        # ocmap,amap,mmap = complexmap(or,deg2rad.(uo),n=2)
        # orianglemap = map(a->HSV(rad2deg(2a),1,1),amap)
        # oripolarmap = map((a,m)->HSV(rad2deg(2a),1,m),amap,mmap)
        
        if plot
            for (p,t) in zip(dp,dpt)
                t01 = clampscale(dogfilter(t),2)
                foreach(ext->save(joinpath(resultdir,"dir_$(p[1])_$(p[2])$ext"),t01),figfmt)
            end
            for (p,t) in zip(op,opt)
                t01 = clampscale(dogfilter(t),2)
                foreach(ext->save(joinpath(resultdir,"ori_$(p[1])_$(p[2])$ext"),t01),figfmt)
            end

            foreach(ext->save(joinpath(resultdir,"dir_anglemap$ext"),diranglemap),figfmt)
            foreach(ext->save(joinpath(resultdir,"dir_polarmap$ext"),dirpolarmap),figfmt)
            foreach(ext->save(joinpath(resultdir,"ori_anglemap$ext"),orianglemap),figfmt)
            foreach(ext->save(joinpath(resultdir,"ori_polarmap$ext"),oripolarmap),figfmt)
        end

        jldsave(joinpath(resultdir,"isi.jld2");cond,epochresponse,dp,dpt,op,opt,dcmap,ocmap,exenv,siteid)
    elseif ex["ID"] in ["ISIEpoch2Color","ISIEpochFlash2Color"]
        t = pairtest(epochresponse,cond.i...).stat
        t_dog = clampscale(dogfilter(t),2)
        
        exenv["color"] = "$(exparam["ColorSpace"])_$(exparam["Color"])"
        mincolor = RGBA(cond.Color[1]...)
        maxcolor = RGBA(cond.Color[2]...)
        exenv["minmaxcolor"] = (;mincolor,maxcolor)
        cm = cgrad(mincolor,maxcolor)
        t_dog_color = map(i->cm[i],t_dog)     
        
        if plot
            foreach(ext->save(joinpath(resultdir,"t$ext"),t_dog),figfmt)
            foreach(ext->save(joinpath(resultdir,"t_color$ext"),t_dog_color),figfmt)
        end
        jldsave(joinpath(resultdir,"isi.jld2");cond,epochresponse,t,exenv,siteid)
    end
    
end
