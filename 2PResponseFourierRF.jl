
# Peichao's Notes:
# 1. Code was written for 2P data (Hartley) from Scanbox. Will export results (dataframe and csv) for plotting.
# 2. If you have multiple planes, it works with splited & interpolated dat. Note results are slightly different.
# 3. If you have single plane, need to change the code (signal and segmentation) a little bit to make it work.

using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact, CSV,MAT, DataStructures, HypothesisTests, StatsFuns, Random, Plots

# Expt info
disk = "O:"
subject = "AF4"  # Animal
recordSession = "004" # Unit
testId = "002"  # Stimulus test

interpolatedData = true   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.

delays = 0:0.05:0.5
ntau = length(collect(delays))
print(collect(delays))
isplot = false

## Prepare data & result path
exptId = join(filter(!isempty,[recordSession, testId]),"_")
dataFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), exptId)
metaFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), "metaFiles")

## load expt, scanning parameters
metaFile=matchfile(Regex("[A-Za-z0-9]*$testId[A-Za-z0-9]*_[A-Za-z0-9]*_meta.mat"),dir=metaFolder,adddir=true)[1]
dataset = prepare(metaFile)
ex = dataset["ex"]
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
blkidx = condidx.>5641  # blanks start from 5641
cidx = .!blkidx
condidx2 = condidx.*cidx + blkidx.* 5641
conduniq = unique(condidx2)
## Load data
segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.segment"),dir=dataFolder,adddir=true)[1]
segment = prepare(segmentFile)
signalFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_merged.signals"),dir=dataFolder,adddir=true)[1]
signal = prepare(signalFile)
# sig = transpose(signal["sig"])   # 1st dimention is cell roi, 2nd is fluorescence trace
spks = transpose(signal["spks"])  # 1st dimention is cell roi, 2nd is spike train

## Load data
planeNum = size(segment["mask"],3)  # how many planes
planeStart = vcat(1, length.(segment["seg_ot"]["vert"]).+1)

## Use for loop process each plane seperately
for pn in 1:planeNum
    # pn=2  # for test
    # Initialize DataFrame for saving results
    recordPlane = string("00",pn-1)  # plane/depth, this notation only works for expt has less than 10 planes
    siteId = join(filter(!isempty,[recordSession, testId, recordPlane]),"_")
    dataExportFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "DataExport")
    resultFolder = joinpath(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "Plots")
    isdir(dataExportFolder) || mkpath(dataExportFolder)
    isdir(resultFolder) || mkpath(resultFolder)
    result = DataFrame()

    cellRoi = segment["seg_ot"]["vert"][pn]
    cellNum = length(cellRoi)
    display("plane: $pn")
    display("Cell Number: $cellNum")

    if interpolatedData
        # rawF = sig[planeStart[pn]:planeStart[pn]+cellNum-1,:]
        spike = spks[planeStart[pn]:planeStart[pn]+cellNum-1,:]
    else
        # rawF = transpose(signal["sig_ot"]["sig"][pn])
        spike = transpose(signal["sig_ot"]["spks"][pn])
    end
    result.py = 0:cellNum-1
    result.ani = fill(subject, cellNum)
    result.dataId = fill(siteId, cellNum)
    result.cellId = 1:cellNum
    ## Chop spk trains according delays
    spk=zeros(nstim,ntau,cellNum)
    for d in eachindex(delays)
        y,num,wind,idx = subrv(sbxft,condon.+delays[d], condoff.+delays[d],isminzero=false,ismaxzero=false,shift=0,israte=false)
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

    taumax=[];kur=[];kurmax=[];kernraw=[];kernnor=[];
    signif=[];sfest=[];oriest=[];
    for i = 1:cellNum
        # i=15
        z = r[:,:,:,i]
        q = reshape(z,szhtly^2,:)   # in this case, there are 61^2 pixels in the stimulus.
        k = dropdims(mapslices(kurtosis,q;dims=1).-3, dims=1) # The kurtosis of any univariate normal distribution is 3. It is common to compare the kurtosis of a distribution to this value.
        tmax = findall(x->x==max(k...),k)[1]
        kmax = max(k...)
        sig = kmax>7
        kernRaw = z[:,:,tmax]  # raw kernel without
        kern = log10.(z[:,:,tmax] ./ z[max_k+1,max_k+1,tmax])

        # estimate ori/sf
        # bw = kern .> (max(kern...) .* 0.95)
        bw = kern .== max(kern...)
        idx = findall(x->x==1,bw)
        foreach(x->if x[1]>31 bw[x]=0 end,idx)
        zzm = sum(sum(zz.*bw,dims=1),dims=2)[1] / (length(idx)/2)

        sf = abs(zzm)/szhtly_visangle  # cyc/deg
        ori = rad2deg(angle(zzm)) # deg

        push!(taumax,tmax);push!(kur,k);push!(kurmax, kmax); push!(kernraw,kernRaw);push!(kernnor,kern);
        push!(signif,sig);push!(oriest,ori); push!(sfest,sf);

        # if sig == true
        #     heatmap(kern,yflip=true, aspect_ratio=:equal,color=:coolwarm)
        #     # plot([0 real(zzm)],[0 imag(zzm)],'wo-','linewidth',3,'markersize',14);
        # end
    end
    result.sig = signif
    result.oriest = oriest
    result.sfest = sfest
    result.taumax = taumax
    result.kernnor = kernnor
    result.kernraw = kernraw
    result.kurmax=kurmax
    result.kur = kur

    #Save results
    CSV.write(joinpath(resultFolder,join([subject,"_",siteId,"_",coneType,"_tuning_result.csv"])), result)
    save(joinpath(dataExportFolder,join([subject,"_",siteId,"_",coneType,"_tuning_result.jld2"])), "result",result)

end
