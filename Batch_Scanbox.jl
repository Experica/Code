using DataFramesMeta,Interact,CSV,MAT,DataStructures,HypothesisTests,StatsFuns,Random

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
    disk=string(param[:dataexportroot][1:2])
    subject=uppercase(ex["Subject_ID"]);recordSession=uppercasefirst(ex["RecordSite"]);testId=ex["TestID"]

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

    if (sbx["line"][1] == 0.00) | (sbx["frame"][1] == 0.00)  # Sometimes there is extra pulse at start, need to remove it
        stNum = 2
    else
        stNum = 1
    end
    trialOnLine = sbx["line"][stNum:2:end]
    trialOnFrame = sbx["frame"][stNum:2:end] + round.(trialOnLine/lineNum)        # if process splitted data use frame_split
    trialOffLine = sbx["line"][stNum+1:2:end]
    trialOffFrame = sbx["frame"][stNum+1:2:end] + round.(trialOffLine/lineNum)    # if process splitted data use frame_split

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
    if interpolatedData
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,join=true)[1]
        signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,join=true)[1]
    else
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
        signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9].signals"),dir=dataFolder,join=true)[1]
    end

    segment = prepare(segmentFile)
    signal = prepare(signalFile)
    sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
    # spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

    planeNum = size(segment["mask"],3)  # how many planes
    if interpolatedData
        planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)
    end
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

        if interpolatedData
            cellRoi = segment["seg_ot"]["vert"][pn]
        else
            cellRoi = segment["vert"]
        end
        cellNum = length(cellRoi)
        display("plane: $pn")
        display("Cell Number: $cellNum")
        cellId = collect(range(1, step=1, stop=cellNum))  # Currently, the cellID is the row index of signal

        if interpolatedData
            rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
            # spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        else
            rawF = sig
            # spike = spks
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
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_result.jld2"])), "result",result,"params", param)
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
    disk=string(param[:dataexportroot][1:2])
    subject=uppercase(ex["Subject_ID"]);recordSession=uppercasefirst(ex["RecordSite"]);testId=ex["TestID"]

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
    trialOffFrame = sbx["frame"][2:2:end] + round.(trialOffLine/lineNum)    # if process splitted data use frame_split

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
    # if hueSpace == "DKL"
    ucid = sort(unique(conditionCond.ColorID))
    hstep = 360/length(ucid)
    conditionCond.ColorID = (conditionCond.ColorID.-minimum(ucid)).*hstep
    conditionCond=rename(conditionCond, :ColorID => :HueAngle)
    # end


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
    if interpolatedData
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,join=true)[1]
        signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,join=true)[1]
    else
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
        signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9].signals"),dir=dataFolder,join=true)[1]
    end

    segment = prepare(segmentFile)
    signal = prepare(signalFile)
    sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
    # spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

    planeNum = size(segment["mask"],3)  # how many planes
    # planeNum = 1
    if interpolatedData
        planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)
    end
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

        if interpolatedData
            cellRoi = segment["seg_ot"]["vert"][pn]
        else
            cellRoi = segment["vert"]
        end
        cellNum = length(cellRoi)
        display("plane: $pn")
        display("Cell Number: $cellNum")
        cellId = collect(range(1, step=1, stop=cellNum))  # Currently, the cellID is the row index of signal

        if interpolatedData
            rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
            # spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        else
            rawF = sig
            # spike = spks
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
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_result.jld2"])), "result",result,"params", param)
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_tuning.jld2"])), "tuning",tempDF)
    end

end

function process_2P_hartleySTA(files,param;uuid="",log=nothing,plot=false)

    # files=tests.files[3]  # for testing

    interpolatedData = haskey(param,:interpolatedData) ? param[:interpolatedData] : true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
    hartelyBlkId = haskey(param,:hartelyBlkId) ? param[:hartelyBlkId] : 5641
    delays = param[:delays]
    stanorm = param[:stanorm]
    stawhiten = param[:stawhiten]
    hartleyscale = param[:hartleyscale]
    print(collect(delays))

    # Expt info
    dataset = prepare(files)
    ex = dataset["ex"];
    disk=string(param[:dataexportroot][1:2])
    subject=uppercase(ex["Subject_ID"]);recordSession=uppercasefirst(ex["RecordSite"]);testId=ex["TestID"]

    ## Prepare data & result path
    exptId = join(filter(!isempty,[recordSession[2:end], testId]),"_")
    dataFolder = joinpath(disk,subject, "2P_data", recordSession, exptId)

    ## load expt, scanning parameters
    envparam = ex["EnvParam"]
    coneType = getparam(envparam,"colorspace")
    sbx = dataset["sbx"]["info"]
    sbxft = ex["frameTimeSer"]   # time series of sbx frame in whole recording
    # Condition Tests
    envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    condon = ex["CondTest"]["CondOn"]
    condoff = ex["CondTest"]["CondOff"]
    condidx = ex["CondTest"]["CondIndex"]
    # condtable = DataFrame(ex["Cond"])
    condtable =  DataFrame(ex["raw"]["log"]["randlog_T1"]["domains"]["Cond"])
    rename!(condtable, [:oridom, :kx, :ky,:bwdom,:colordom])

    # find out blanks and unique conditions
    blkidx = condidx .>= hartelyBlkId  # blanks start from 5641
    cidx = .!blkidx
    condidx2 = condidx.*cidx + blkidx.* hartelyBlkId   # all blanks now have the same id 5641

    ## Load data
    if interpolatedData
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,join=true)[1]
        signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,join=true)[1]
    else
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
        signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9].signals"),dir=dataFolder,join=true)[1]
    end
    segment = prepare(segmentFile)
    signal = prepare(signalFile)
    # sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
    spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

    ##
    # Prepare Imageset
    nscale = haskey(param,:nscale) ? param[:nscale] : 1
    downsample = haskey(param,:downsample) ? param[:downsample] : 2
    sigma = haskey(param,:sigma) ? param[:sigma] : 1.5
    # bgRGB = [getparam(envparam,"backgroundR"),getparam(envparam,"backgroundG"),getparam(envparam,"backgroundB")]
    bgcolor = RGBA([0.5,0.5,0.5,1]...)
    coneType = string(getparam(envparam,"colorspace"))
    masktype = getparam(envparam,"mask_type")
    maskradius = getparam(envparam,"mask_radius")
    masksigma = 1#getparam(envparam,"Sigma")
    xsize = getparam(envparam,"x_size")
    ysize = getparam(envparam,"y_size")
    stisize = xsize
    ppd = haskey(param,:ppd) ? param[:ppd] : 46
    imagesetname = "Hartley_stisize$(stisize)_scale$(hartleyscale)"
    maskradius = haskey(param,:maskradius) ? param[:maskradius] : 0.16 #maskradius/stisize
    # if coneType == "L"
    #     maxcolor = RGBA()
    #     mincolor = RGBA()
    # elseif coneType == "M"
    #
    # elseif coneType == "S"
    #
    # end

    if !haskey(param,imagesetname)
        imageset = map(i->GrayA.(hartley(kx=i.kx,ky=i.ky,bw=i.bwdom,stisize=stisize, ppd=ppd,norm=false,scale=hartleyscale)),eachrow(condtable))
        # imageset = map(i->GrayA.(grating(θ=deg2rad(i.Ori),sf=i.SpatialFreq,phase=rem(i.SpatialPhase+1,1)+0.02,stisize=stisize,ppd=23)),eachrow(condtable))
        imageset = Dict{Symbol,Any}(:pyramid => map(i->gaussian_pyramid(i, nscale-1, downsample, sigma),imageset))
        imageset[:imagesize] = map(i->size(i),imageset[:pyramid][1])
        param[imagesetname] = imageset
    end

    # Prepare Image Stimuli
    imageset = param[imagesetname]
    bgcolor = oftype(imageset[:pyramid][1][1][1],bgcolor)
    unmaskindex = map(i->alphamask(i,radius=maskradius,sigma=masksigma,masktype=masktype)[2],imageset[:pyramid][1])
    imagestimuli = map(s->map(i->alphablend.(alphamask(i[s],radius=maskradius,sigma=masksigma,masktype=masktype)[1],[bgcolor]),imageset[:pyramid]),1:nscale)


    ## Load data
    planeNum = size(segment["mask"],3)  # how many planes
    if interpolatedData
        planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)
    end

    ## Use for loop process each plane seperately
    for pn in 1:planeNum
        # pn=1  # for test
        display("plane: $pn")
        # Initialize DataFrame for saving results
        recordPlane = string("00",pn-1)  # plane/depth, this notation only works for expt has less than 10 planes
        siteId = join(filter(!isempty,[recordSession[2:end], testId, recordPlane]),"_")
        dataExportFolder = joinpath(disk,subject, "2P_analysis", recordSession, siteId, "DataExport")
        resultFolder = joinpath(disk,subject, "2P_analysis", recordSession, siteId, "Plots")
        isdir(dataExportFolder) || mkpath(dataExportFolder)
        isdir(resultFolder) || mkpath(resultFolder)

        if interpolatedData
            cellRoi = segment["seg_ot"]["vert"][pn]
        else
            cellRoi = segment["vert"]
        end
        cellNum = length(cellRoi)
        display("Cell Number: $cellNum")

        if interpolatedData
            # rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
            spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        else
            # rawF = sig
            spike = spks
        end

        ## Calculate STA
        if :STA in param[:model]
            scaleindex=1
            imagesize = imageset[:imagesize][scaleindex]
            xi = unmaskindex[scaleindex]
            uci = unique(condidx2)
            ucii = map(i->findall(condidx2.==i),deleteat!(uci,findall(isequal(hartelyBlkId),uci)))   # find the repeats of each unique condition
            ubii = map(i->findall(condidx2.==i), [hartelyBlkId])


            cx = Array{Float64}(undef,length(ucii),length(xi))
            foreach(i->cx[i,:]=gray.(imagestimuli[scaleindex][uci[i]][xi]),1:size(cx,1))

            uy = Array{Float64}(undef,cellNum,length(delays),length(ucii))
            ucy = Array{Float64}(undef,cellNum,length(delays),length(ucii))
            uby = Array{Float64}(undef,cellNum,length(delays),length(ubii))
            usta = Array{Float64}(undef,cellNum,length(delays),length(xi))

            for d in eachindex(delays)
                display("Processing delay: $d")
                y,num,wind,idx = epochspiketrain(sbxft,condon.+delays[d], condoff.+delays[d],isminzero=false,ismaxzero=false,shift=0,israte=false)
                spk=zeros(size(spike,1),length(idx))
                for i =1:length(idx)
                    spkepo = @view spike[:,idx[i][1]:idx[i][end]]
                    spk[:,i]= mean(spkepo, dims=2)
                end
                for cell in 1:cellNum
                    # display(cell)
                    cy = map(i->mean(spk[cell,:][i]),ucii)  # response to grating
                    bly = map(i->mean(spk[cell,:][i]),ubii) # response to blank, baseline
                    ry = cy.-bly  # remove baseline
                    csta = sta(cx,ry,norm=stanorm,whiten=stawhiten)  # calculate sta
                    ucy[cell,d,:]=cy
                    uby[cell,d,:]=bly
                    uy[cell,d,:]=ry
                    usta[cell,d,:]=csta
                    if plot
                        csta = (csta.+1)./2
                        csta = fill(0.5,imagesize)
                        csta[xi] = csta
                        r = [extrema(csta)...]
                        title = "$(ugs[cell])Unit_$(unitid[cell])_STA_$(delays[d])"
                        p = plotsta(csta,imagesize=imagesize,stisize=stisize,index=xi,title=title,r=r)
                        foreach(i->save(joinpath(resultFolder,"$title$i"),p),[".png"])
                    end
                end
            end
            save(joinpath(dataExportFolder,join([subject,"_",siteId,"_",coneType,"_sta.jld2"])),"imagesize",imagesize,"cx",cx,"xi",xi,"xcond",condtable[uci,:],"uy",uy,"ucy",ucy,"usta",usta,"uby",uby,"delays",delays,"maskradius",maskradius,"stisize",stisize,"color",coneType)
        end
    end

end


function process_2P_hartleyFourier(files,param;uuid="",log=nothing,plot=false)

    interpolatedData = haskey(param,:interpolatedData) ? param[:interpolatedData] : true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
    hartelyBlkId = haskey(param,:hartelyBlkId) ? param[:hartelyBlkId] : 5641
    delays = param[:delays]
    ntau = length(collect(delays))
    print(collect(delays))

    # Expt info
    dataset = prepare(files)
    ex = dataset["ex"]
    disk=string(param[:dataexportroot][1:2])
    subject=uppercase(ex["Subject_ID"]);recordSession=uppercasefirst(ex["RecordSite"]);testId=ex["TestID"]

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
    blkidx = condidx.>hartelyBlkId  # blanks start from 5641
    cidx = .!blkidx
    condidx2 = condidx.*cidx + blkidx.* hartelyBlkId
    conduniq = unique(condidx2)
    ## Load data
    if interpolatedData
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,join=true)[1]
        signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,join=true)[1]
    else
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
        signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9].signals"),dir=dataFolder,join=true)[1]
    end
    segment = prepare(segmentFile)
    signal = prepare(signalFile)
    # sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
    spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

    ## Load data
    planeNum = size(segment["mask"],3)  # how many planes
    if interpolatedData
        planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)
    end

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

        if interpolatedData
            cellRoi = segment["seg_ot"]["vert"][pn]
        else
            cellRoi = segment["vert"]
        end
        cellNum = length(cellRoi)
        display("plane: $pn")
        display("Cell Number: $cellNum")

        if interpolatedData
            # rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
            spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        else
            # rawF = sig
            spike = spks
        end
        result.py = 0:cellNum-1
        result.ani = fill(subject, cellNum)
        result.dataId = fill(siteId, cellNum)
        result.cellId = 1:cellNum
        ## Chop spk trains according delays
        spk=zeros(nstim,ntau,cellNum)
        for d in eachindex(delays)
            y,num,wind,idx = epochspiketrain(sbxft,condon.+delays[d], condoff.+delays[d],isminzero=false,ismaxzero=false,shift=0,israte=false)
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

        taumax=[];kstd=[];kstdmax=[];kernraw=[];kernnor=[];kernest=[];
        kdelta=[];signif=[];slambda=[];sfmax=[];orimax=[];sfmean=[];orimean=[];
        sfidx=[];sfcurve=[];oriidx=[];oricurve=[];
        for i = 1:cellNum
            # i=15
            z = r[:,:,:,i]
            q = reshape(z,szhtly^2,:)   # in this case, there are 61^2 pixels in the stimulus.
            # k = dropdims(mapslices(kurtosis,q;dims=1).-3, dims=1) # The kurtosis of any univariate normal distribution is 3. It is common to compare the kurtosis of a distribution to this value.
            k = [std(q[:,j]) for j in 1:size(q,2)]
            tmax = findall(x->x==max(k...),k)[1]
            kmax = max(k...)
            # sig = kmax>7
            kernRaw = z[:,:,tmax]  # raw kernel without blank normalization
            kern = log10.(z[:,:,tmax] ./ z[max_k+1,max_k+1,tmax])   # kernal normalized by blank

            # separability measure and estimate kernel
            u,s,v = svd(kernRaw)
            s = Diagonal(s)
            lambda = s[1,1]/s[2,2]
            q = s[1,1]
            s = zeros(size(s))
            s[1,1] = q
            kest = u*s*v'   # estimated kernel

            # energy measure
            delta = kmax / k[1] - 1
            sig = delta > 0.25

            # find the maxi/best condition
            # bw = kern .> (max(kern...) .* 0.95)
            bwmax = kern .== max(kern...)
            idxmax = findall(x->x==1,bwmax)
            foreach(x->if x[1]>(max_k+1) bwmax[x]=0 end,idxmax)   # choose upper quadrants
            # estimate ori/sf by max
            zzm = sum(sum(zz.*bwmax,dims=1),dims=2)[1] / (length(idxmax)/2)
            sf_max = abs(zzm)/szhtly_visangle  # cyc/deg
            ori_max = rad2deg(angle(zzm)) # deg

            # find the maxi/best condition
            bw = kern .> (max(kern...) .* 0.95)
            idx = findall(x->x==1,bw)
            foreach(x->if x[1]>(max_k+1) bw[x]=0 end,idx)   # choose upper quadrants
            # estimate ori/sf by mean
            zzm = sum(sum(zz.*bw,dims=1),dims=2)[1] / (length(idx)/2)
            sf_mean = abs(zzm)/szhtly_visangle  # cyc/deg
            ori_mean = rad2deg(angle(zzm)) # deg

            # Ori tuning curve
            sf_best = max((abs.(zz).*bwmax)...)
            idxsf = findall(x->x==sf_best,abs.(zz))
            filter!(x->x[1]<=(max_k+1),idxsf)

            ori_idx=rad2deg.(angle.(reverse(zz[idxsf])))
            ori_curve=reverse(kern[idxsf])

            # SF tuning curve
            ori_best = max((angle.(zz).*bwmax)...)
            idxori = findall(x->x==ori_best,angle.(zz))

            sf_idx = (abs.(zz[idxori]))./szhtly_visangle
            sf_curve = kern[idxori]
            idxinf = findall(x->!isinf(x),sf_curve)
            sf_idx = sort(sf_idx[idxinf])
            if idxori[end][2] < max_k
                sf_curve = reverse(sf_curve)
            end
            filter!(x->x!=-Inf,sf_curve)

            push!(taumax,tmax);push!(kstd,k);push!(kstdmax, kmax); push!(kernraw,kernRaw);push!(kernnor,kern);
            push!(kernest,kest);push!(signif,sig);push!(kdelta,delta); push!(slambda,lambda);
            push!(orimax,ori_max); push!(sfmax,sf_max);push!(orimean,ori_mean); push!(sfmean,sf_mean);
            push!(oriidx, ori_idx);push!(oricurve,ori_curve); push!(sfidx,sf_idx); push!(sfcurve,sf_curve)

            # if sig == true
            #     heatmap(kern,yflip=true, aspect_ratio=:equal,color=:coolwarm)
            #     # plot([0 real(zzm)],[0 imag(zzm)],'wo-','linewidth',3,'markersize',14);
            # end
        end
        result.signif = signif
        result.taumax = taumax
        result.kstdmax = kstdmax
        result.kdelta = kdelta
        result.slambda = slambda
        result.orimax = orimax
        result.sfmax = sfmax
        result.orimean = orimean
        result.sfmean = sfmean

        result1=copy(result)

        result.kstd = kstd
        result.kernnor = kernnor
        result.kernraw = kernraw
        result.kernest = kernest
        result.oriidx = oriidx
        result.oricurve = oricurve
        result.sfidx = sfidx
        result.sfcurve = sfcurve

        #Save results
        CSV.write(joinpath(resultFolder,join([subject,"_",siteId,"_",coneType,"_tuning_result.csv"])), result1)
        save(joinpath(dataExportFolder,join([subject,"_",siteId,"_",coneType,"_tuning_result.jld2"])), "result",result,"params",param)
    end
end
