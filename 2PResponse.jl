
# Peichao's Notes:
# 1. Code was written for 2P data (Direction-Spatial Frequency test) from Scanbox. Will export results (dataframe and csv) for plotting.
# 2. If you have multiple planes, it works with splited & interpolated dat. Note results are slightly different.

using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact, CSV,MAT, DataStructures, HypothesisTests, StatsFuns, Random

# Expt info
disk = "O:"
subject = "AF4"  # Animal
recordSession = "003" # Unit
testId = "010"  # Stimulus test

interpolatedData = true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
preOffset = 0.1
responseOffset = 0.05  # in sec
α = 0.05   # p value
sampnum = 100   # random sampling 100 times
fitThres = 0.5
isplot = false
# ismodulative([DataFrame(Y=condResp) ctc[cti,:]], alpha=α, interact=false)
## Prepare data & result path
exptId = join(filter(!isempty,[recordSession, testId]),"_")
dataFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), exptId)
metaFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), "metaFiles")

## load expt, scanning parameters
# metaFile=matchfile(Regex("$subject*_$recordSession*_$testId*_ot_meta.mat"),dir=metaFolder,adddir=true)[1]
metaFile=matchfile(Regex("$subject*_$recordSession*_$testId*_ot_meta.mat"),dir=metaFolder,adddir=true)[1]
dataset = prepare(metaFile)
ex = dataset["ex"]
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

## Load data
if interpolatedData
    segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,adddir=true)[1]
    signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,adddir=true)[1]
else
    segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*.segment"),dir=dataFolder,adddir=true)[1]
    signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9].signals"),dir=dataFolder,adddir=true)[1]
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
    siteId = join(filter(!isempty,[recordSession, testId, recordPlane]),"_")
    dataExportFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "DataExport")
    resultFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "Plots")
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
    if isplot
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
        # display(cell)
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
    if isplot
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
        # display(cell)
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
        push!(ufs[f],factorresponsefeature(dropmissing(mseuc)[f],dropmissing(mseuc)[:m],factor=f, isfit=oriAUC[cell]>fitThres))
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
    if isplot
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
display("Processing Is Done!!!!!") #

# Plot Spike Train for all trials of all cells
# epochext = preicidur
# @manipulate for cell in 1:cellNum
# ys,ns,ws,is = subrv(spike[cell,:],condOn,condOff,isminzero=true,shift=0)
# plotspiketrain(ys,timeline=[0,minCondDur],title="Unit_$(unitid[u])")
# end

# for u in 1:length(unitspike)
# ys,ns,ws,is = subrv(unitspike[u],condOn.-epochext,condOff.+epochext,isminzero=true,shift=epochext)
# plotspiketrain(ys,timeline=[0,minCondDur],title="Unit_$(unitid[u])")
# foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_SpikeTrian$i")),[".png",".svg"])
# end

# foreach(i->ctcli[condition[i,:i]] .= i, 1:condNum)
# vcat(condition[1:end-1, :i]...)

# Tuning map
# plotunitposition(unitposition,color=map(i->HSV(2*i.oo,1,1-i.ocv),ufs[:Ori]),alpha=1)
# foreach(i->savefig(joinpath(resultdir,"UnitPosition_OriTuning$i")),[".png",".svg"])
# plotunitposition(unitposition,color=map(i->HSV(i.od,1,1-i.dcv),ufs[:Ori]),alpha=1)
# foreach(i->savefig(joinpath(resultdir,"UnitPosition_DirTuning$i")),[".png",".svg"])
# save(joinpath(resultdir,"factorresponse.jld2"),"factorstats",ufs,"fms",fms,"fses",fses,"fa",fa)

# @df DataFrame(ufs[:Ori]) corrplot([:oo :od :ocv :dcv],nbins=30)
