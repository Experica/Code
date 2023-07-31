using NeuroAnalysis,StatsBase,JLD2,YAML,FileWatching,DataFrames,Images,HypothesisTests

function online_epoch_imager(testroot,resultroot;dt=5,tout=15,figfmt = [".png"],showprogress=true)
    test = splitdir(testroot)[end]
    exfile = joinpath(testroot,"$test.yaml")
    metafile = joinpath(testroot,"$test.meta")

    printstyled("Start Online Analysis: $test\n",color=:yellow,reverse=true)
    while true
        timedwait(()->isfile(metafile),tout;pollint=dt) == :ok && break
        printstyled("No experiment and meta data detected, continue checking? (y|n)> ",color=:blue)
        if readline() != "y"
            printstyled("Abort Online Analysis: $test\n",color=:red,reverse=true); return
        end
    end

    printstyled("Preparing experiment and meta data ...\n",color=:cyan)
    mdata = load(Val(Imager),metafile)
    edata = load(Val(Experica),exfile)

    ex = edata["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
    siteresultdir = joinpath(resultroot,subject,siteid)
    resultdir = joinpath(siteresultdir,test);mkpath(resultdir)

    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"];conddur = ex["CondDur"]
    exenv=Dict{Any,Any}("ID"=>ex["ID"],"eye"=>ex["Eye"])
    conddesign = DataFrame(ex["Cond"]) |> sort!
    if ex["ID"] == "ISIEpochOri8"
        ds = conddesign.Ori.-90
        os = unique(conddesign.Ori.%180)
        qs = unique(os.%90)
        diranglemap=dirpolarmap=orianglemap=oripolarmap=nothing
    end
    
    fmt = mdata["meta"]["DataFormat"]
    w = mdata["meta"]["ImageFormat"]["Width"]
    h = mdata["meta"]["ImageFormat"]["Height"]
    pixfmt = mdata["meta"]["ImageFormat"]["PixelFormat"]
    fps = mdata["meta"]["AcquisitionControl"]["AcquisitionFrameRate"]
    delay = 500 # hemodynamic delay
    nbaseframe = round(Int,preicidur/1000*fps)
    nframepcond = round(Int,conddur/1000*fps)
    condstartframe = round(Int,(preicidur+delay)/1000*fps)
    baseframeindex = 1:nbaseframe
    frameindex = (0:nframepcond-1).+condstartframe
    ei = 0
    repeat = 1
    imagefile = Vector{String}[]
    epochresponse = Matrix{Float64}[]
    cr = Matrix{Float64}[]
    ct = DataFrame()
    ctc= DataFrame()

    while true
        epochpath = joinpath(testroot,".Epoch$ei")
        # assume experiment stopped when requested epoch file not saved in timeout
        timedwait(()->isfile(epochpath),tout;pollint=dt) == :ok || break
        printstyled("Checking Epoch$ei ...\n",color=:blue)
        t = YAML.load_file(epochpath);ci=t["CondIndex"]+1
        push!(ct,(CondIndex=ci,CondRepeat=t["CondRepeat"]))
        push!(ctc,(;(f=>conddesign[ci,f] for f in propertynames(conddesign))...))

        epochdir = joinpath(testroot,"Epoch$ei")
        epochfiles = [joinpath(epochdir,"$test-Frame$i.$fmt") for i in eachindex(readdir(epochdir))]
        push!(imagefile,epochfiles)
        push!(epochresponse, frameresponse_imager(epochfiles,w,h,frameindex,baseframeindex))
        ei+=1

        count(i->i==repeat,ct.CondRepeat) == nrow(conddesign) || continue
        printstyled("All conditions have repeated $repeat times ...\n",color=:red)
        ri = findall(ct.CondRepeat.==repeat)
        repeat+=1
        
        if ex["ID"] == "ISIEpochOri8"
            corder = sortperm(ctc.Ori[ri])
            rs = @views epochresponse[ri[corder]]
            if isempty(cr)
                append!(cr,rs);continue
            end
            erdir = joinpath(resultdir,"Epoch0-$ei");mkpath(erdir)
            foreach(i->cr[i] = (cr[i].+rs[i])/2,eachindex(cr))
            
            dcmap,amap,mmap = complexmap(cr,deg2rad.(ds))
            diranglemap = map(a->HSV(rad2deg(a),1,1),amap)
            dirpolarmap = map((a,m)->HSV(rad2deg(a),1,m),amap,mmap)

            ocr = @views map(i->.+(cr[indexin([i,i+180],conddesign.Ori)]...)/2,os)
            ocmap,amap,mmap = complexmap(ocr,2deg2rad.(os))
            orianglemap = map(a->HSV(rad2deg(a),1,1),amap)
            showprogress && display(orianglemap)
            oripolarmap = map((a,m)->HSV(rad2deg(a),1,m),amap,mmap)

            foreach(ext->save(joinpath(erdir,"dir_anglemap$ext"),diranglemap),figfmt)
            foreach(ext->save(joinpath(erdir,"dir_polarmap$ext"),dirpolarmap),figfmt)
            foreach(ext->save(joinpath(erdir,"ori_anglemap$ext"),orianglemap),figfmt)
            foreach(ext->save(joinpath(erdir,"ori_polarmap$ext"),oripolarmap),figfmt)
            jldsave(joinpath(erdir,"isi.jld2");ctc,epochresponse,exenv,siteid)
        end
        
    end

    if ex["ID"] == "ISIEpochOri8"
        foreach(ext->save(joinpath(resultdir,"dir_anglemap$ext"),diranglemap),figfmt)
        foreach(ext->save(joinpath(resultdir,"dir_polarmap$ext"),dirpolarmap),figfmt)
        foreach(ext->save(joinpath(resultdir,"ori_anglemap$ext"),orianglemap),figfmt)
        foreach(ext->save(joinpath(resultdir,"ori_polarmap$ext"),oripolarmap),figfmt)
        jldsave(joinpath(resultdir,"isi.jld2");ctc,epochresponse,exenv,siteid)
    end
    
    printstyled("Finish Online Analysis: $test\n",color=:green,reverse=true)
end

function online_cycle_imager(testroot,resultroot;dt=5,tout=15,eachcycle=true,figfmt = [".png"],showprogress=true)
    test = splitdir(testroot)[end]
    exfile = joinpath(testroot,"$test.yaml")
    metafile = joinpath(testroot,"$test.meta")

    printstyled("Start Online Analysis: $test\n",color=:yellow,reverse=true)
    while true
        timedwait(()->isfile(metafile),tout;pollint=dt) == :ok && break
        printstyled("No experiment and meta data detected, continue checking? (y|n)> ",color=:blue)
        if readline() != "y"
            printstyled("Abort Online Analysis: $test\n",color=:red,reverse=true); return
        end
    end

    printstyled("Preparing experiment and meta data ...\n",color=:cyan)
    mdata = load(Val(Imager),metafile)
    edata = load(Val(Experica),exfile)

    ex = edata["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
    siteresultdir = joinpath(resultroot,subject,siteid)
    resultdir = joinpath(siteresultdir,test);mkpath(resultdir)

    envparam = ex["EnvParam"];exparam = ex["Param"];preicidur = ex["PreICI"]
    modulatefreq = envparam["ModulateTemporalFreq"]
    modulatetype = envparam["ModulateType"]
    modulateduty = envparam["ModulateDuty"]
    exenv=Dict{Any,Any}("ID"=>ex["ID"],"eye"=>ex["Eye"])
    if ex["ID"] == "ISICycle2Color"
        F1polarity = F1polarity_dog = F1polarity_ahe = nothing
        # cos phase(maxcolor:0->mincolor:π->maxcolor:2π) to linear normalized polarity(mincolor:0->maxcolor:1)
        pha2pol(p,s=0) = abs.(mod2pi.(p .+ s) .- π) ./ π
        
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
    end

    epochdir = joinpath(testroot,"Epoch0")
    fmt = mdata["meta"]["DataFormat"]
    w = mdata["meta"]["ImageFormat"]["Width"]
    h = mdata["meta"]["ImageFormat"]["Height"]
    pixfmt = mdata["meta"]["ImageFormat"]["PixelFormat"]
    fps = mdata["meta"]["AcquisitionControl"]["AcquisitionFrameRate"]
    delay = 500 # hemodynamic delay
    freq = modulatefreq
    nbaseframe = round(Int,preicidur/1000*fps)
    nframepcycle = round(Int,1/freq*fps)
    cyclestartframe = round(Int,(preicidur+delay)/1000*fps)
    df = round(Int,dt*fps) # number of frames to check in each loop
    baseresponse = nothing
    fi = 0
    ncycle = 0
    cyclestopframe = cyclestartframe+nframepcycle-1
    imagefile = String[]
    F1=zeros(ComplexF64,h,w)
    
    while true
        if timedwait(()->isfile(joinpath(epochdir,"$test-Frame$(fi+df).$fmt")),tout;pollint=dt) == :ok
            fr = fi:fi+df
        else # assume experiment stopped when requested frame not saved in timeout
            lf = length(readdir(epochdir))-1
            lf>=fi ? (fr = fi:lf) : break
        end
        printstyled("Checking Saved Frames $fr ...\n",color=:blue)
        append!(imagefile,[joinpath(epochdir,"$test-Frame$f.$fmt") for f in fr])
        fi=fr[end]+1

        if isnothing(baseresponse)
            if nbaseframe-1 <= fr[end]
                baseresponse = dropdims(mean(readrawim_Mono12Packed(imagefile[1:nbaseframe],w,h),dims=3),dims=3) # Float64
            end
            continue
        end
        cyclestopframe > fr[end] && continue
        ncycle+=1

        printstyled("Cycle $ncycle Frames Ready ...\n",color=:red)
        crdir = joinpath(resultdir,"Frame$(cyclestartframe)-$(cyclestopframe)");mkpath(crdir)
        if eachcycle # only current cycle
            F1 .+= dft_imager(imagefile[cyclestartframe:cyclestopframe],w,h,fps,baseresponse,freq)[1]
            cyclestartframe=cyclestopframe+1
        else # accumulate to current cycle
            F1 = dft_imager(imagefile[cyclestartframe:cyclestopframe],w,h,fps,baseresponse,freq)[1]
        end
        cyclestopframe+=nframepcycle

        F1phase = angle.(F1) # [-π, π]
        F1mag = abs.(F1)
        F1phase01 = mod2pi.(F1phase) ./ (2π)
        F1mag01 = clampscale(F1mag,3)
        foreach(ext->save(joinpath(crdir,"F1Phase$ext"),F1phase01),figfmt)
        foreach(ext->save(joinpath(crdir,"F1Mag$ext"),F1mag01),figfmt)

        if ex["ID"] == "ISICycle2Color"
            F1polarity = pha2pol(F1phase)
            F1polarity_dog = clampscale(dogfilter(F1polarity),2)
            showprogress && display(Gray.(F1polarity_dog))
            F1polarity_ahe = ahe(F1polarity)

            foreach(ext->save(joinpath(crdir,"F1Polarity_dog$ext"),F1polarity_dog),figfmt)
            foreach(ext->save(joinpath(crdir,"F1Polarity_ahe$ext"),F1polarity_ahe),figfmt)
            jldsave(joinpath(crdir,"isi.jld2");freq,F1,F1polarity,exenv,siteid)
        end
    end

    if ex["ID"] == "ISICycle2Color"
        F1polarity_dog_color = map(a->HSV(cmfun(a/2,p=0.5)),F1polarity_dog)
        foreach(ext->save(joinpath(resultdir,"F1Polarity_dog_Color$ext"),F1polarity_dog_color),figfmt)
        foreach(ext->save(joinpath(resultdir,"F1Polarity_dog$ext"),F1polarity_dog),figfmt)
        foreach(ext->save(joinpath(resultdir,"F1Polarity_ahe$ext"),F1polarity_ahe),figfmt)
        jldsave(joinpath(resultdir,"isi.jld2");imagefile,freq,F1,F1polarity,exenv,siteid)
    end
    
    printstyled("Finish Online Analysis: $test\n",color=:green,reverse=true)
end
