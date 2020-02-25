
# Peichao's Notes:
# 1. Code was written for 2P data (Direction-Spatial frequency-Hue test) from Scanbox. Will export results (dataframe and csv) for plotting.
# 2. If you have multiple planes, it works with splited & interpolated dat. Note results are slightly different.
# 3. If you have single plane, need to change the code (signal and segmentation) a little bit to make it work.

using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact,CSV,MAT,DataStructures,HypothesisTests,StatsFuns,Random

# User input and Expt info
disk = "O:"
subject = "AF4"  # Animal
recordSession = "005" # Unit
testId = "007"  # Stimulus test
hueSpace = "DKL"   # Color space used? DKL or HSL
interpolatedData = true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.
preOffset = 0.1  # in sec
responseOffset = 0.05  # in sec
α = 0.05   # p value
diraucThres = 0.8   # if passed, calculate hue direction, otherwise calculate hue axis
oriaucThres = 0.5
# Respthres = 0.1  # Set a response threshold to filter out low response cells?
sampnum = 100   # random sampling 100 times
blankId = 36  # Blank Id
excId = [27,28,blankId]  # Exclude some condition?
isplot = false  # Plot figures to investigate?

## Prepare data & result path
exptId = join(filter(!isempty,[recordSession, testId]),"_")
dataFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), exptId)
metaFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), "metaFiles")

## load expt, scanning parameters
metaFile=matchfile(Regex("[A-Za-z0-9]*$testId[A-Za-z0-9]*_[A-Za-z0-9]*_meta.mat"),dir=metaFolder,adddir=true)[1]
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
planeStart = 1

## Use for loop process each plane seperately
for pn in 1:planeNum
    # pn=1  # for test
    global planeStart
    display(planeStart)
    # Initialize DataFrame for saving results
    recordPlane = string("00",pn-1)  # plane/depth, this notation only works for expt has less than 10 planes
    siteId = join(filter(!isempty,[recordSession, testId, recordPlane]),"_")
    dataExportFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "DataExport")
    resultFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "Plots")
    isdir(dataExportFolder) || mkpath(dataExportFolder)
    isdir(resultFolder) || mkpath(resultFolder)
    result = DataFrame()

    cellRoi = segment["seg_ot"]["vert"][pn]   # ???? Note: this vert structure was saved for Python, col and row are reversed.
    cellNum = length(cellRoi)
    display(cellNum)
    cellId = collect(range(1, step=1, stop=cellNum))  # Currently, the cellID is the row index of signal

    if interpolatedData
        rawF = sig[planeStart:planeStart+cellNum-1,:]
        # spike = spks[planeStart:planeStart+cellNum-1,:]
    else
        rawF = transpose(signal["sig_ot"]["sig"][pn])
        # spike = transpose(signal["sig_ot"]["spks"][pn])
    end

    planeStart = planeStart+cellNum   # update, only works when planeStart is globalized
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
    if isplot
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
        # resu= [factorresponsestats(mcti[:dir],resp[:,t],factor=:dir) for t in 1:mcti.n[1]]
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
                resu= [factorresponsestats(mcti[:Dir],resp[:,t],factor=:Dir) for t in 1:mcti[1,:n]]
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

                resu= [factorresponsestats(mcti[:HueAngle],resp[:,t],factor=:HueAngle) for t in 1:mcti[1,:n]]
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
        push!(ufs[f],factorresponsestats(dropmissing(mseuc)[f],dropmissing(mseuc)[:m],factor=f, thres=max(oriAUC[cell], dirAUC[cell], hueaxAUC[cell], huedirAUC[cell])))
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
    # isplot = true
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
    save(joinpath(dataExportFolder,join([subject,"_",siteId,"_tuning.jld2"])), "tuning",tempDF)
end
planeStart = 1  # no clear function like Matab, reset it mannually
