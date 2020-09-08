
# Peichao's Notes:
# 1. Code was written for 2P data (Hartley) from Scanbox. Will export results (dataframe and csv) for plotting.
# 2. If you have multiple planes, it works with splited & interpolated dat. Note results are slightly different.
# 3. If you have single plane, set interpolatedData as false.

using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,LinearAlgebra,Images,StatsBase,Interact, CSV,MAT, DataStructures, HypothesisTests, StatsFuns, Random, Plots

# Expt info
disk = "K:"
subject = "AE6"  # Animal
recordSession = "002" # Unit
testId = "005"  # Stimulus test

interpolatedData = false   # If you have multiplanes. True: use interpolated data; false: use uniterpolated data. Results are slightly different.

delays = -0.066:0.066:0.4
ntau = length(collect(delays))
print(collect(delays))
isplot = false

## Prepare data & result path
exptId = join(filter(!isempty,[recordSession, testId]),"_")
dataFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), exptId)
metaFolder = joinpath(disk,subject, "2P_data", join(["U",recordSession]), "metaFiles")

## load expt, scanning parameters
metaFile=matchfile(Regex("[A-Za-z0-9]*_[A-Za-z0-9]*_$testId*_meta.mat"),dir=metaFolder,join=true)[1]
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
    ## Chop spk trains according to delays
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
            rf = rf + rot180(rf) # average over phases, PL
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
        # i=43
        z = r[:,:,:,i]
        q = reshape(z,szhtly^2,:)   # in this case, there are 61^2 pixels in the stimulus.
        # k = dropdims(mapslices(kurtosis,q;dims=1).-3, dims=1) # The kurtosis of any univariate normal distribution is 3. It is common to compare the kurtosis of a distribution to this value.
        k = [std(q[:,j]) for j in 1:size(q,2)]
        tmax = findall(x->x==max(k...),k)[1]
        kmax = max(k...)
        # sig = kmax>7
        kernRaw = z[:,:,tmax]  # raw kernel without blank normalization
        kern = log10.(z[:,:,tmax] ./ z[max_k+1,max_k+1,tmax])  # kernal normalized by blank
        replace!(kern, -Inf=>0)

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
            # heatmap(kmask,yflip=true, aspect_ratio=:equal,color=:coolwarm)
            # plot([0 real(zzm)],[0 imag(zzm)],'wo-','linewidth',3,'markersize',14);
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
    save(joinpath(dataExportFolder,join([subject,"_",siteId,"_",coneType,"_tuning_result.jld2"])), "result",result)

end
