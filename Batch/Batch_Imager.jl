using Images

function dft_imager(files,w,h,fs,base,f...)
    N = length(files)
    ks = round.(Int,f./fs.*N)
    Fs = [zeros(ComplexF64,h,w) for _ in f]
    Ω = [exp(-im*2π*n/N) for n in 0:(N-1)]
    p = ProgressMeter.Progress(N,desc="DFT at $f Hz ")
    ls = [ReentrantLock() for _ in f]
    @inbounds Threads.@threads for n in 0:(N-1)
        img = readrawim_Mono12Packed(files[n+1],w,h)
        # img = img ./ base .- 1 # reduce DC, increase SNR
        img = log2.(img./base) # reduce DC, increase SNR
        @inbounds for i in eachindex(f)
            lock(ls[i]) do
                Fs[i] .+= img .* Ω[((n*ks[i])%N)+1]
            end
        end
        next!(p)
    end
    return Fs
end

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
    condoff = ex["CondTest"]["CondOff"]
    modulatefreq = envparam["ModulateTemporalFreq"]
    modulatetype = envparam["ModulateType"]
    modulateduty = envparam["ModulateDuty"]
    ncycle = ex["CondRepeat"]
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    exenv=Dict()
    exenv["eye"] = ex["Eye"]

    # Prepare Frame
    imagefile = dataset["imagefile"]
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
    F1mag01 = F1mag./maximum(F1mag)
    F2mag01 = F2mag./maximum(F2mag)

    if plot
        foreach(ext->save(joinpath(resultdir,"F1Phase$ext"),F1phase01),figfmt)
        foreach(ext->save(joinpath(resultdir,"F1Mag$ext"),F1mag01),figfmt)
        foreach(ext->save(joinpath(resultdir,"F1PhaseMag$ext"),GrayA.(F1phase01,F1mag01)),figfmt)
        foreach(ext->save(joinpath(resultdir,"F2Phase$ext"),F2phase01),figfmt)
        foreach(ext->save(joinpath(resultdir,"F2Mag$ext"),F2mag01),figfmt)
        foreach(ext->save(joinpath(resultdir,"F2PhaseMag$ext"),GrayA.(F2phase01,F2mag01)),figfmt)
    end

    if ex["ID"] == "ISICycleOri"
        orianglemap = map(a->HSV(360a,1,1),F2phase01)
        oripolarmap = map((a,m)->HSV(360a,1,m),F2phase01,F2mag01)
        diranglemap = map(a->HSV(360a,1,1),F1phase01)
        dirpolarmap = map((a,m)->HSV(360a,1,m),F1phase01,F1mag01)
        
        F2mag01e = adjust_histogram(F2mag01, AdaptiveEqualization(nbins = 256, rblocks = 16, cblocks = 16, clip = 0.9))
        clamp01!(F2mag01e)
        F2phase01e = adjust_histogram(F2phase01, AdaptiveEqualization(nbins = 256, rblocks = 16, cblocks = 16, clip = 0.9))
        clamp01!(F2phase01e)

        F1mag01e = adjust_histogram(F1mag01, AdaptiveEqualization(nbins = 256, rblocks = 16, cblocks = 16, clip = 0.9))
        clamp01!(F1mag01e)
        F1phase01e = adjust_histogram(F1phase01, AdaptiveEqualization(nbins = 256, rblocks = 16, cblocks = 16, clip = 0.9))
        clamp01!(F1phase01e)

        orianglemape = map(a->HSV(360a,1,1),F2phase01e)
        oripolarmape = map((a,m)->HSV(360a,1,m),F2phase01e,F2mag01e)
        diranglemape = map(a->HSV(360a,1,1),F1phase01e)
        dirpolarmape = map((a,m)->HSV(360a,1,m),F1phase01e,F1mag01e)
        
        if plot
            foreach(ext->save(joinpath(resultdir,"ori_anglemap$ext"),orianglemap),figfmt)
            foreach(ext->save(joinpath(resultdir,"ori_polarmap$ext"),oripolarmap),figfmt)
            foreach(ext->save(joinpath(resultdir,"dir_anglemap$ext"),diranglemap),figfmt)
            foreach(ext->save(joinpath(resultdir,"dir_polarmap$ext"),dirpolarmap),figfmt)
            foreach(ext->save(joinpath(resultdir,"ori_anglemap_Enhanced$ext"),orianglemape),figfmt)
            foreach(ext->save(joinpath(resultdir,"ori_polarmap_Enhanced$ext"),oripolarmape),figfmt)
            foreach(ext->save(joinpath(resultdir,"dir_anglemap_Enhanced$ext"),diranglemape),figfmt)
            foreach(ext->save(joinpath(resultdir,"dir_polarmap_Enhanced$ext"),dirpolarmape),figfmt)
        end
        jldsave(joinpath(resultdir,"isi.jld2");freqs,Fs,F1phase01,F2phase01,F1mag01,F2mag01,exenv,siteid)
    elseif ex["ID"] == "ISICycleColorPlane"

    elseif ex["ID"] == "ISICycle2Color"
        # cos phase(maxcolor:0->mincolor:π->maxcolor:2π) to linear normalized polarity(mincolor:0->maxcolor:1)
        pha2pol(p,s=0) = abs.(mod2pi.(p .+ s) .- π) ./ π
        F1polarity = pha2pol(F1phase)
        F1polaritye = adjust_histogram(F1polarity, AdaptiveEqualization(nbins = 256, rblocks = 20, cblocks = 20, clip = 0.1))
        clamp01!(F1polaritye)

        F2polarity = pha2pol(F2phase)
        F2polaritye = adjust_histogram(F2polarity, AdaptiveEqualization(nbins = 256, rblocks = 20, cblocks = 20, clip = 0.1))
        clamp01!(F2polaritye)
        
        exenv["color"] = "$(exparam["ColorSpace"])_$(exparam["Color"])"
        mincolor = RGBA(envparam["MinColor"]...)
        maxcolor = RGBA(envparam["MaxColor"]...)
        exenv["minmaxcolor"] = (;mincolor,maxcolor)
        colordist = maxcolor-mincolor
        if modulatetype == "Sinusoidal"
            # Cos modulation between MinColor and MaxColor
            cmfun = (a;f=1,p=0) -> mincolor + colordist * ((cos(2π*(f*a+p))+1)/2)
        elseif modulatetype == "Square"
            # Cos Sign modulation between MinColor and MaxColor
            cmfun = (a;f=1,p=0) -> mincolor + colordist * ((sign(modulateduty-mod(f*a+p+0.25,1))+1)/2)
        end
        F1polaritycolor = map(a->HSV(cmfun(a/2,p=0.5)),F1polaritye)
        F1polaritycolormag = map((a,m)->HSV(a.h,a.s,m),F1polaritycolor,F1mag01)

        if plot
            foreach(ext->save(joinpath(resultdir,"F1Polarity$ext"),F1polarity),figfmt)
            foreach(ext->save(joinpath(resultdir,"F1Polarity_Enhanced$ext"),F1polaritye),figfmt)
            foreach(ext->save(joinpath(resultdir,"F2Polarity$ext"),F2polarity),figfmt)
            foreach(ext->save(joinpath(resultdir,"F2Polarity_Enhanced$ext"),F2polaritye),figfmt)
            foreach(ext->save(joinpath(resultdir,"F1PolarityColor$ext"),F1polaritycolor),figfmt)
            foreach(ext->save(joinpath(resultdir,"F1PolarityColorMag$ext"),F1polaritycolormag),figfmt)
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
    condoff = ex["CondTest"]["CondOff"]
    ctc = condtestcond(ex["CondTestCond"])
    cond = condin(ctc)
    factors = finalfactor(ctc)
    isbalance = length(unique(cond.n)) == 1
    figfmt = haskey(param,:figfmt) ? param[:figfmt] : [".png"]
    exenv=Dict()
    exenv["eye"] = ex["Eye"]

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
    responsedelay = haskey(param,:responsedelay) ? param[:responsedelay] : 500 # hemodynamic delay
    baseframeindex = epoch2sampleindex([0 preicidur],framerate,maxsampleindex=minframe)
    frameindex = epoch2sampleindex([0 conddur].+(preicidur+responsedelay),framerate,maxsampleindex=minframe)
    epochresponse = Array{Float64}(undef,h,w,nepoch)
    p = ProgressMeter.Progress(nepoch,desc="Epoch Response ... ")
    @inbounds Threads.@threads for i in 1:nepoch
        epochframes = readrawim_Mono12Packed(imagefile[i],w,h)
        epochresponse[:,:,i] = frameresponse(epochframes;frameindex,baseframeindex)
        next!(p)
    end

    # if ex["ID"] == "ISIEpochOri8"
    #     factors = factors[1]
    #     if isbalance
    #         responses = epochresponse
    #         angles = deg2rad.(ctc[!,factors])
    #     else
    #         responses = condimageresponse(epochresponse,cond.i)
    #         angles = deg2rad.(cond[!,factors])
    #     end
    #     maps = complexmap(responses,angles)
    # end
    #
    # if plot
    #     normmmap = maps.mmap/maximum(maps.mmap)
    #     oriamap = mod.(maps.amap,π)
    #     diranglemap = map(a->HSV(rad2deg(a),1,1),maps.amap)
    #     dirpolarmap = map((a,m)->HSV(rad2deg(a),1,m),maps.amap,normmmap)
    #     orianglemap = map(a->HSV(rad2deg(2a),1,1),oriamap)
    #     oripolarmap = map((a,m)->HSV(rad2deg(2a),1,m),oriamap,normmmap)
    #     foreach(ext->save(joinpath(resultdir,"dir_anglemap$ext"),diranglemap),figfmt)
    #     foreach(ext->save(joinpath(resultdir,"dir_polarmap$ext"),dirpolarmap),figfmt)
    #     foreach(ext->save(joinpath(resultdir,"ori_anglemap$ext"),orianglemap),figfmt)
    #     foreach(ext->save(joinpath(resultdir,"ori_polarmap$ext"),oripolarmap),figfmt)
    # end
    jldsave(joinpath(resultdir,"isi.jld2");epochresponse,exenv,siteid)
end
