
# Notes: 1. Test with splited data, so load splitted signal and mat files

using NeuroAnalysis,Statistics,DataFrames,StatsPlots,Mmap,Images,StatsBase,Interact, MAT, DataStructures

# Expt info
disk = "O:"
subject = "AF4"  # Animal
recordSession = "005" # Unit
recordSite = "001" # plane/depth
testId = "008"  # Stimulus test

preOffset = 0.1
responseDelay = 0.1

## Prepare data path
siteId = join(filter(!isempty,[recordSession, testId, recordSite]),"_")
dataFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), recordSite, siteId)
metaFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), "metaFiles")
dataExportFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), recordSite, siteId, "DataExport")
resultFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), recordSite, siteId, "Results")
isdir(dataExportFolder) || mkpath(dataExportFolder)
isdir(resultFolder) || mkpath(resultFolder)

## load data, expt, scanning parameters
metaFile=matchfile(Regex("[A-Za-z0-9]*$testId[A-Za-z0-9]*_[A-Za-z0-9]*_meta.mat"),dir=metaFolder,adddir=true)[1]
dataset = prepare(metaFile)
ex = dataset["ex"]
envparam = ex["EnvParam"]
sbx = dataset["sbx"]["info"]

segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.segment"),dir=dataFolder,adddir=true)[1]
segment = prepare(segmentFile)
cellRoi = transpose(segment["vert"])   # ???? Note: this vert structure was saved for Python, col and row are reversed.
cellId = collect(range(1, step=1, stop=length(cellRoi)))
roiMap = segment["mask"]

signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.signals"),dir=dataFolder,adddir=true)[1]
signal = prepare(signalFile)
rawF = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
spike = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

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
minTrialDur = minimum(trialOffFrame-trialOnFrame)
# histogram(trialOffFrame-trialOnFrame,nbins=20,title="Trial Duration(Set to $minTrialDur)")

# On/off frame indces of condations/stimuli
preStim = ex["PreICI"]; stim = ex["CondDur"]; postStim = ex["SufICI"]
trialOnTime = fill(0, size(trialEpoch,1))
condOfftime = preStim + stim
preEpoch = [0 preStim-preOffset]
condEpoch = [preStim condOfftime]
preFrame=epoch2samplerange(preEpoch, sbxfs)
condFrame=epoch2samplerange(condEpoch, sbxfs)
preOn = fill(preFrame.start, size(trialEpoch,1))
preOff = fill(preFrame.stop, size(trialEpoch,1))
condOn = fill(condFrame.start, size(trialEpoch,1))
condOff = fill(condFrame.stop, size(trialEpoch,1))
# trialTime =  sum([preStim, stim, postStim])
# endPre = floor(minTrialDur * (preStim/trialTime))  # in frame numbers
# startPre = floor(preStimSample * preOffsetRatio)  # offset of start of pre stim, in frames
# startStim = ceil(minTrialDur * (preStim + responseDelay)/trialTime)  # Add a response delay, so there are several frames dropped after start stim.
# endStim = ceil(minTrialDur * (preStim+stim) / trialTime)  # in frames

# Plot dF/F traces of all trials for all cells
# Cut raw fluorescence traces according to trial on/off time and calculate dF/F
ys = sbxsubrm(rawF,trialEpoch,cellId;fun=dPrime(preFrame))
# Plot
@manipulate for e in 1:size(ys,1)
    plotanalog(transpose(ys[e,:,:]), fs=sbxfs, timeline=condEpoch.-preStim, xunit=:s, ystep=1,cunit=:p, color=:fire,xext=preStim)
end

# Plot Spike Train for all trials of all cells
# epochext = preicidur
@manipulate for cell in 1:size(spike,1)
ys,ns,ws,is = subrv(spike[cell,:],condOn,condOff,isminzero=true,shift=0)
plotspiketrain(ys,timeline=[0,minCondDur],title="Unit_$(unitid[u])")
end

# for u in 1:length(unitspike)
# ys,ns,ws,is = subrv(unitspike[u],condOn.-epochext,condOff.+epochext,isminzero=true,shift=epochext)
# plotspiketrain(ys,timeline=[0,minCondDur],title="Unit_$(unitid[u])")
# foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_SpikeTrian$i")),[".png",".svg"])
# end
DataFrame(ex["CondTestCond"])
# Condition Test
# Prepare stimulus conditions
meanF = dropdims(mean(ys[:,condFrame,:], dims=2), dims=2)   # Mean response within stim time of all cells for all trials
ctc = DataFrame(ex["CondTestCond"])
factors = finalfactor(ctc)
ci = ctc[factors[1]].!="blank"
cctc = ctc[ci,factors]   # Remove blank
# t = [ctc DataFrame(i=1:nrow(ctc))]
# t = by(t, names(ctc),g->DataFrame(n=nrow(g), i=[g[:,:i]]))
# aa=sort!(t[1:3,:]);
# sort!(t[1:48,:]);
ccond = condin(cctc)
dFMean = Array{Float64}(undef, size(ccond,1))
dFSem = Array{Float64}(undef, size(ccond,1))
sigCellResp = DataFrame(dFMean=dFMean, dFSem=dFSem)
cellResp = OrderedDict("ccond"=>ccond)
test = Array{Int64}(undef, size(ccond,1))
for con in 1:48
    test[con]= findmax(ccond.i[con])[1]
end

tempdict = Dict(i=>meanF[i,ci] for i in 1:size(meanF,1))
mseuc = condresponse(tempdict,condin(cctc))
@df mseuc plot(:SpatialFreq,:m,group=:u)

# for cell in 1:size(meanF,1)
#     for cond in 1:size(ccond,1)
#         sigCellResp.dFMean[cond] = mean(meanF[cell,:][ccond.i[cond]])
#         sigCellResp.dFSem[cond] = sem(meanF[cell,:][ccond.i[cond]])
#     end
#     push!(cellResp, string(cell)=>sigCellResp)
# end


# insertcols!(cconF, 5, :meanResp => [1, 2])
factors = condfactor(ccond)
aa=flin(cctc)
fl = flin(ccond[factors])
fa = OrderedDict(f=>fl[f][f] for f in keys(fl))
ccondOn=condOn[ci]
ccondOff=condOff[ci]

# Prepare blank conditions
bi = .!ci    # Choose blank
bctc = ctc[bi,factors]
bcond = condin(bctc)
bcondOn=condOn[bi]
bcondOff=condOff[bi]

# eval.([:($(Symbol(k))=spike[$k]) for k in keys(spike)])
# condOn = ex["CondTest"]["condOn"]
# condOff = ex["CondTest"]["condOff"]
# CSV.write("dataframe.csv", df)

# Condition Response
responsedelay=0.015
@manipulate for u in 1:length(unitspike)
rs = subrvr(unitspike[u],ccondOn.+responsedelay,ccondOff.+responsedelay)
plotcondresponse(rs,cctc)
end

for u in 1:length(unitspike)
rs = subrvr(unitspike[u],ccondOn.+responsedelay,ccondOff.+responsedelay)
plotcondresponse(rs,cctc)
foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_CondResponse$i")),[".png",".svg"])
end

if !isempty(bcondOn)
    ubr=[]
    for u in 1:length(unitspike)
        br = subrvr(unitspike[u],bcondOn.+responsedelay,bcondOff.+responsedelay)
        mseuc = condresponse(br,condin(bctc))+
        push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
    end
end

# Condition Response in Factor Space
fms,fses,fa=factorresponse(unitspike,cctc,ccondOn,ccondOff,responsedelay=responsedelay)

@manipulate for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
    fd = findfirst(f.==keys(fa))
    fdn = length(fa[f])
    p[fd]=1:fdn
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(unitid[u],fdn))
    mseuc[f]=fa[f]
    plotcondresponse(dropmissing(mseuc),colors=[:black],responseline=ubr[u:u])
end

ufs = Dict(k=>[] for k in keys(fa))
for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
    fd = findfirst(f.==keys(fa))
    fdn = length(fa[f])
    p[fd]=1:fdn
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(unitid[u],fdn))
    mseuc[f]=fa[f]

    push!(ufs[f],factorresponsestats(mseuc[f],mseuc[:m],factor=f))

    # plotcondresponse(dropmissing(mseuc),colors=[:black],responseline=ubr[u:u])
    # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_$(f)_Tuning$i")),[".png",".svg"])
end


# Tuning map
plotunitposition(unitposition,color=map(i->HSV(2*i.oo,1,1-i.ocv),ufs[:Ori]),alpha=1)
foreach(i->savefig(joinpath(resultdir,"UnitPosition_OriTuning$i")),[".png",".svg"])
plotunitposition(unitposition,color=map(i->HSV(i.od,1,1-i.dcv),ufs[:Ori]),alpha=1)
foreach(i->savefig(joinpath(resultdir,"UnitPosition_DirTuning$i")),[".png",".svg"])
save(joinpath(resultdir,"factorresponse.jld2"),"factorstats",ufs,"fms",fms,"fses",fses,"fa",fa)


@df DataFrame(ufs[:Ori]) corrplot([:oo :od :ocv :dcv],nbins=30)
