
# Notes: 1. Test with splited data, so load splitted signal and mat files

using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact, CSV,MAT, DataStructures, HypothesisTests, StatsFuns

# Expt info
disk = "O:"
subject = "AF4"  # Animal
recordSession = "005" # Unit
recordPlane = "001" # plane/depth
testId = "007"  # Stimulus test

preOffset = 0.1
# responseDelay = 0.1
α = 0.05   # p value
isplot = false

## Prepare data & result path
siteId = join(filter(!isempty,[recordSession, testId, recordPlane]),"_")
dataFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), recordPlane, siteId)
metaFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), "metaFiles")
dataExportFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), recordPlane, siteId, "DataExport")
resultFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), recordPlane, siteId, "Results")
isdir(dataExportFolder) || mkpath(dataExportFolder)
isdir(resultFolder) || mkpath(resultFolder)

## Initialize DataFrame for saving results
result = DataFrame()

## load data, expt, scanning parameters
metaFile=matchfile(Regex("[A-Za-z0-9]*$testId[A-Za-z0-9]*_[A-Za-z0-9]*_meta.mat"),dir=metaFolder,adddir=true)[1]
dataset = prepare(metaFile)
ex = dataset["ex"]
envparam = ex["EnvParam"]
sbx = dataset["sbx"]["info"]

segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.segment"),dir=dataFolder,adddir=true)[1]
segment = prepare(segmentFile)
cellRoi = segment["vert"]   # ???? Note: this vert structure was saved for Python, col and row are reversed.
cellId = collect(range(1, step=1, stop=length(cellRoi)))  # Currently, the cellID is the row index of signal
# roiMap = segment["mask"]

signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.signals"),dir=dataFolder,adddir=true)[1]
signal = prepare(signalFile)
rawF = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
spike = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train
cellNum = size(rawF,1)
result.ani = fill(subject, cellNum)
result.dataId = fill(siteId, cellNum)
result.cellId = 1:cellNum

## Align Scan frames with stimulus
# Load the scan parameters
scanFreq = sbx["resfreq"]
lineNum = sbx["sz"][1]
if haskey(sbx, "recordsPerBuffer_bi")
   scanMode = 1  # bidirectional scanning   ###################### Split
else
   scanMode = 1  # unidirectional scanning
end
sbxfs = 1/(lineNum/scanFreq/scanMode)   # frame rate
trialOnLine = sbx["line"][1:2:end]
trialOnFrame = sbx["frame_split"][1:2:end] + round.(trialOnLine/lineNum)        ############################  split
trialOffLine = sbx["line"][2:2:end]
trialOffFrame = sbx["frame_split"][2:2:end] + round.(trialOnLine/lineNum)    ############################  split

# On/off frame indces of trials
trialEpoch = Int.(hcat(trialOnFrame, trialOffFrame))
# minTrialDur = minimum(trialOffFrame-trialOnFrame)
# histogram(trialOffFrame-trialOnFrame,nbins=20,title="Trial Duration(Set to $minTrialDur)")

# Trials ==> Condition
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
condEpoch = [preStim condOfftime]
preFrame=epoch2samplerange(preEpoch, sbxfs)
condFrame=epoch2samplerange(condEpoch, sbxfs)
# preOn = fill(preFrame.start, trialNum)
# preOff = fill(preFrame.stop, trialNum)
# condOn = fill(condFrame.start, trialNum)
# condOff = fill(condFrame.stop, trialNum)


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

## Put each cell's response in factor space (dir-sf...)
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
if isplot
    @manipulate for cell in 1:cellNum
        heatmap(fms[cell])
    end
end

## Get the responsive cells
uresponsive=[];umodulative=[]
cti = reduce(append!,condition[1:end-1, :i],init=Int[])
for cell in 1:cellNum
    condResp = cellMeanTrial[cell,cti]
    push!(umodulative,ismodulative([DataFrame(Y=condResp) ctc[cti,:]], alpha=α, interact=true))
    blankResp = cellMeanTrial[cell,condition[end,:i]]
    isresp = []
    for i in 1:condNum
        condResp = cellMeanTrial[cell,condition[i, :i]]
        push!(isresp, pvalue(UnequalVarianceTTest(blankResp,condResp))<α)
    end
    push!(uresponsive,any(isresp))
    # plotcondresponse(condResp ctc)
    # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[cell])_CondResponse$i")),[".png",".svg"])
end
visResp = uresponsive .| umodulative
display(["uresponsive:", count(uresponsive)])
display(["umodulative:", count(umodulative)])
display(["Responsive cells:", count(visResp)])
result.visResp = visResp
result.modulative = umodulative

## Check which cell is significantly tuning by orientation or direction
oripvalue=[];dirpvalue=[];
for cell in 1:cellNum
    # cell=1
    # Get all trial Id of under maximal sf
    # mcti = @where(condition, :SpatialFreq .== ufm[:SpatialFreq][cell])
    mcti = condition[condition.SpatialFreq .==ufm[:SpatialFreq][cell], :]
    resp=[cellMeanTrial[cell,mcti.i[r][t]] for r in 1:nrow(mcti), t in 1:mcti.n[1]]
    resu= [factorresponsestats(mcti[:dir],resp[:,t],factor=:dir) for t in 1:mcti.n[1]]
    orivec = reduce(vcat,[resu[t].oov for t in 1:mcti.n[1]])
    orip = hotellingt2test([real(orivec) imag(orivec)],[0 0],0.05)
     # check significance of direction selective
    oriang = angle(mean(orivec, dims=1)[1])  # ori angle in ori space
    orivecdir = exp(im*oriang/2)   # ori vector in dir space
    dirvec = reduce(vcat,[resu[t].odv for t in 1:mcti.n[1]])
    dirp = dirsigtest(orivecdir, dirvec)
    push!(oripvalue,orip);push!(dirpvalue,dirp);
end
result.orip = oripvalue
result.dirp = dirpvalue

temp=DataFrame(temp)
cellMeanTrial[cell, collect(Iterators.flatten(mcti.i[1,:]))]'
collect(Iterators.flatten(mcti.i[:,:]))



## Get the optimal factor level using Circular Variance for each cell
ufs = Dict(k=>[] for k in keys(fa))
for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...] # Replace missing with -Inf, then find the x-y coordinates of max value.
    fd = findfirst(f.==keys(fa))   # facotr dimention
    fdn = length(fa[f])  # dimention length/number of factor level
    p[fd]=1:fdn   # pick up a slice for plotting tuning curve
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(cellId[u],fdn),ug=fill(parse(Int, recordPlane), fdn))  # make DataFrame for plotting
    mseuc[f]=fa[f]
    # The optimal dir, ori (based on circular variance) and sf (based on log10 fitting)
    push!(ufs[f],factorresponsestats(dropmissing(mseuc)[f],dropmissing(mseuc)[:m],factor=f))

    plotcondresponse(dropmissing(mseuc),colors=[:black],projection=[],responseline=[], responsetype=:ResponseF)
    # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_$(f)_Tuning$i")),[".png"]#,".svg"])
end
result.optsf = ufs[:SpatialFreq]
tempDF=DataFrame(ufs[:dir])
result.optdir = tempDF.od
result.dirmag = tempDF.odr
result.dircv = tempDF.dcv
result.optori = tempDF.oo
result.orimag = tempDF.oor
result.oricv = tempDF.ocv

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

# Fitting direction and orientation tuning


#Save results
CSV.write(joinpath(resultFolder,join([subject,"_",siteId,"_result.csv"])), result)
save(joinpath(resultFolder,join([subject,"_",siteId,"_result.jld2"])), "result",result)



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

# Condition Test
# tempdict = Dict(cellId[i]=>cellMeanTrial[i,:] for i in 1:cellNum)
# mseuc = condresponse(tempdict,condition)
# @df mseuc plot(:dir,:m,group=:u)

# factors = finalfactor(ctc)
# ci = ctc[factors[1]].!="blank"
# cctc = ctc[ci,factors]   # Remove blank
# ccond = condin(cctc)
# dFMean = Array{Float64}(undef, size(ccond,1))
# dFSem = Array{Float64}(undef, size(ccond,1))
# sigCellResp = DataFrame(dFMean=dFMean, dFSem=dFSem)
# cellResp = OrderedDict("ccond"=>ccond)

# for cell in 1:size(cellMeanTrial,1)
#     for cond in 1:size(ccond,1)
#         sigCellResp.dFMean[cond] = mean(cellMeanTrial[cell,:][ccond.i[cond]])
#         sigCellResp.dFSem[cond] = sem(cellMeanTrial[cell,:][ccond.i[cond]])
#     end
#     push!(cellResp, string(cell)=>sigCellResp)
# end
# insertcols!(cconF, 5, :meanResp => [1, 2])


# foreach(i->ctcli[condition[i,:i]] .= i, 1:condNum)
# vcat(condition[1:end-1, :i]...)


# Tuning map
# plotunitposition(unitposition,color=map(i->HSV(2*i.oo,1,1-i.ocv),ufs[:Ori]),alpha=1)
# foreach(i->savefig(joinpath(resultdir,"UnitPosition_OriTuning$i")),[".png",".svg"])
# plotunitposition(unitposition,color=map(i->HSV(i.od,1,1-i.dcv),ufs[:Ori]),alpha=1)
# foreach(i->savefig(joinpath(resultdir,"UnitPosition_DirTuning$i")),[".png",".svg"])
# save(joinpath(resultdir,"factorresponse.jld2"),"factorstats",ufs,"fms",fms,"fses",fses,"fa",fa)
#
#
# @df DataFrame(ufs[:Ori]) corrplot([:oo :od :ocv :dcv],nbins=30)
