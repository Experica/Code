# PL: This code is to organize 2P data collected from one location (maybe multiple planes) for following plotting code
# PL: Some experiments used single-plane scanning at mulitple depths (no z scanning), so data from multiplane may have different session and/or expt ids.
#      This code is writtern for this situation.
# PL: Data will be organized at two levels: 1. within plane and  2. across planes
# PL: This code assumes you did ori, hue and L, M, S, Achromatic hartley in the same session.
# Remember: fitting is done when auc>0.5
#

using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact,CSV,MAT,Query,DataStructures, HypothesisTests, StatsFuns, Random
# Expt info
disk = "O:"
subject = "AF4"  # Animal
# recordPlane = ["001", "000"]  # I hard-code the plane Id here. The Id here is matched with the naming rule of multiple plane scanning
recordPlane = ["001"]
oriaucThres = 0.7
diraucThres = 0.7
hueaucThres = 0.8

oricvThres = 0.7
dircvThres = 0.7
huecvThres = 0.7

osiThres = 0.5
dsiThres = 0.5
huesiThres = 0.5
cpiThres = 0.33

numberofColor = 12
colorSpace = "HSL"   # DKL, # HSL

addSTA = false
addFourier = false
staThres = 0.25

## Make folders and path
mainpath = joinpath(disk,subject, "2P_analysis")
dataExportFolder1 = joinpath(mainpath, "Summary", "DataExport")
isdir(dataExportFolder1) || mkpath(dataExportFolder1)

meta = readmeta(joinpath(dataExportFolder1,"metadata.mat"))
for t in 1:nrow(meta)
    meta.files[t]=meta.files[t][3:end]
    meta.Subject_ID[t]=uppercase(meta.Subject_ID[t])
end

## Query Tests
tests = @from i in meta begin
    @where i.Subject_ID == subject
    @where i.sourceformat == "Scanbox"
    @select {i.ID,i.RecordSite,i.filename}
    @collect DataFrame
    end
sort!(tests, [:RecordSite, :filename])
oritests = @from i in tests begin
    # @where i.RecordSite == "u005" #|| i.RecordSite == "u012"
    # @where i.RecordSite != "u001"
    # @where i.RecordSite == "u012"
    @where i.ID == "DirSF"
    @select {i.ID,i.RecordSite,i.filename}
    @collect DataFrame
    end

huetests = @from i in tests begin
    # @where i.RecordSite == "u012"
    @where i.RecordSite != "u008"
    # @where i.RecordSite == "u005" #|| i.RecordSite == "u012"
    @where i.ID == "DirSFColor"
    @select {i.ID,i.RecordSite,i.filename}
    @collect DataFrame
    end

hartleytests = @from i in tests begin
    # @where i.RecordSite != "u001"
    # @where i.RecordSite != "u005"
    @where i.ID == "Hartley"
    @select {i.ID,i.RecordSite,i.filename}
    @collect DataFrame
    end

##
recordSession = unique!(oritests.RecordSite) # Unit
recordSession = map(i->i[2:4],recordSession)
oriExptId = oritests.filename
hueExptId = huetests.filename # Stimulus test
huerefOri = [3 8]
hartleyExptId = hartleytests.filename    # Unit Ids have L, M, S, and Achromatic hartleys.
exptOriNum = length(oriExptId)
exptHueNum = length(hueExptId)
planeNum = Int64(length(oriExptId) / length(recordSession))
sessNum = length(recordSession)*planeNum
if length(recordPlane) != planeNum
    @warn "Warning: More plane acquired than you think!!!!"
    exit()
end

paraName = string("oa",string(oriaucThres),"da",string(diraucThres),"ha",string(hueaucThres),"oc",string(oricvThres),
           "dc",string(dsiThres),"hc",string(huecvThres),"osi",string(osiThres),"dsi",string(dsiThres),"hsi",string(huesiThres),"cpi",string(cpiThres))

## Direction/Orientation data
rois=Any[]; bkgs=Any[]; oriData=DataFrame(); oriSum=DataFrame(); hueData=DataFrame(); hueSum=DataFrame();
staData=DataFrame(); staSum=DataFrame();
for i = 1:exptOriNum
    # oriData=DataFrame();
    display("Processing Oriexpt No: $i")
    j=1#map(j->isodd(i) ? 1 : 2, i)
    local unitId
    local siteId
    unitId = oriExptId[i][5:7]
    siteId = join([unitId, "_",recordPlane[j]])
    dataFolder = joinpath(mainpath, join(["U", unitId]), join([oriExptId[i][5:11], "_", recordPlane[j]]), "DataExport")
    dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,join=true)[1]

    segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
    alignFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*.align"),dir=dataFolder,join=true)[1]

    result = load(dataFile)["result"]
    segment = prepare(segmentFile)
    align = prepare(alignFile)
    cellRoi = segment["vert"]
    bkgImg = align["m"]
    locID = repeat([siteId],length(cellRoi))
    insertcols!(result,4,locId=locID)
    append!(oriData, result)
    push!(rois, cellRoi)
    push!(bkgs, bkgImg)
    summ=DataFrame(id=result.dataId[1], exptid=oriExptId[i], planeid=recordPlane[j], cellNum=length(cellRoi), visResp=sum(result.visResp),
         oriResp=sum((result.oriauc.>oriaucThres).&result.visResp), oriRespfit=sum((result.oriauc.>oriaucThres).&result.visResp.&(.~isnan.(result.fitori))), oriSelcv=sum((result.oricv.<oricvThres).&(result.oriauc.>oriaucThres).&result.visResp), oriSelosi=sum((result.osi.>osiThres).&(result.oriauc.>oriaucThres).&result.visResp),
         dirResp=sum((result.dirauc.>diraucThres).&result.visResp), dirRespfit=sum((result.dirauc.>diraucThres).&result.visResp.&(.~isnan.(result.fitdir))), dirSelcv=sum((result.dircv.<dircvThres).&(result.dirauc.>diraucThres).&result.visResp), dirSeldsi=sum((result.dsi.>dsiThres).&(result.dirauc.>diraucThres).&result.visResp))
    append!(oriSum, summ)
    dataExportFolder2 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
    isdir(dataExportFolder2) || mkpath(dataExportFolder2)
    # save(joinpath(dataExportFolder2,join([subject,"_",siteId,"_roibkg.jld2"])), Dict("roi"=>cellRoi, "bkg"=>bkgImg))   # only need save one for each plane in each unit
    CSV.write(joinpath(dataExportFolder2,join([subject,"_",siteId,"_oriData.csv"])), oriData)
    save(joinpath(dataExportFolder2,join([subject,"_",siteId,"_roibkg.jld2"])), "roi",cellRoi, "bkg",bkgImg)
end

## Hue data
for i = 1:exptHueNum
    # hueData=DataFrame();
    j=1#map(j->isodd(i) ? 1 : 2, i)
    local unitId
    local siteId
    unitId = hueExptId[i][5:7]
    siteId = join([unitId, "_",recordPlane[j]])
    # dataFolder = joinpath(mainpath, join(["U", unitId]), join([oriExptId[huerefOri[i]][5:11], "_", recordPlane[j]]), "DataExport")  # load ori data for cpi calculation
    dataFolder = joinpath(mainpath, join(["U", unitId]), join([oriExptId[i][5:11], "_", recordPlane[j]]), "DataExport")  # load ori data for cpi calculation
    dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,join=true)[1]
    resultori = load(dataFile)["result"]
    dataFolder = joinpath(mainpath, join(["U", unitId]), join([hueExptId[i][5:11], "_", recordPlane[j]]), "DataExport")
    dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,join=true)[1]
    result = load(dataFile)["result"]
    locID = repeat([join([unitId, "_", recordPlane[j]])],size(result)[1])
    insertcols!(result,4,locId=locID)
    hueResp = 1 .- minimum(hcat(result.hueaxcv, result.huedicv),dims=2)
    achroResp = 1 .- resultori.oricv
    cpi = dropdims((hueResp .- achroResp) ./ (hueResp .+ achroResp),dims=2)
    insertcols!(result,11,huecpi=cpi)
    append!(hueData, result)
    summ=DataFrame(id=result.dataId[1], exptid=hueExptId[i], planeid=recordPlane[j], cellNum=size(result)[1], visResp=sum(result.visResp),
         hueaxResp=sum((result.hueaxauc.>hueaucThres).&result.visResp), hueaxRespfit=sum((result.hueaxauc.>hueaucThres).&result.visResp.&(.~isnan.(result.fithueax))),hueaxSelcv=sum((result.hueaxcv.<huecvThres).&(result.hueaxauc.>hueaucThres).&result.visResp), hueaxSelsi=sum((result.hueaxsi.>huesiThres).&(result.hueaxauc.>hueaucThres).&result.visResp),
         huediResp=sum((result.huediauc.>hueaucThres).&result.visResp), huediRespfit=sum((result.huediauc.>hueaucThres).&result.visResp.&(.~isnan.(result.fithuedi))),huediSelcv=sum((result.huedicv.<huecvThres).&(result.huediauc.>hueaucThres).&result.visResp), huediSelsi=sum((result.huedisi.>huesiThres).&(result.huediauc.>hueaucThres).&result.visResp),
         hueoriResp=sum((result.hueoriauc.>oriaucThres).&result.visResp), hueoriRespfit=sum((result.hueoriauc.>oriaucThres).&result.visResp.&(.~isnan.(result.fithueori))),hueoriSelcv=sum((result.hueoricv.<oricvThres).&(result.hueoriauc.>oriaucThres).&result.visResp), hueoriSeldsi=sum((result.hueosi.>osiThres).&(result.hueoriauc.>oriaucThres).&result.visResp),
         huedirResp=sum((result.huedirauc.>diraucThres).&result.visResp), huedirRespfit=sum((result.huedirauc.>diraucThres).&result.visResp.&(.~isnan.(result.fithuedir))),huedirSelcv=sum((result.huedircv.<dircvThres).&(result.huedirauc.>diraucThres).&result.visResp), huedirSeldsi=sum((result.huedsi.>dsiThres).&(result.huedirauc.>diraucThres).&result.visResp),
         hueaucResp=sum(((result.hueaxauc.>hueaucThres).|(result.huediauc.>hueaucThres)).&result.visResp), hueSelcpi=sum(((result.hueaxauc.>hueaucThres).|(result.huediauc.>hueaucThres)).&result.visResp.&(cpi.>cpiThres)))
    append!(hueSum, summ)
    dataExportFolder2 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
    isdir(dataExportFolder2) || mkpath(dataExportFolder2)
    CSV.write(joinpath(dataExportFolder2,join([subject,"_",siteId,"_hueData.csv"])), hueData)
end

## Load Hartely (STA) data
if addSTA
    for i = 1:sessNum
        display("Plane Num: $i")
        local unitId
        j=1#map(j->isodd(i) ? 1 : 2, i)
        # unitId = hartleyExptId[i][5:7]
        unitId = recordSession[i]
        plId = recordPlane[j]
        dataFolder = joinpath(mainpath, join(["U",unitId]),"_Summary","DataExport")  # load ori data for cpi calculation
        dataFile=matchfile(Regex("$subject*_*$unitId*_$plId*_thres$staThres*_sta_datasetFinal.jld2"),dir=dataFolder,join=true)[1]
        display(dataFile)
        sta = load(dataFile)["datasetFinal"]
        cellNum = length(sta["ulsta"])
        ulsta = sort!(OrderedDict(sta["ulsta"]))
        if isequal(j,1)
            global planeData, dataExportFolder2
            planeData = DataFrame()
            dataExportFolder2 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
            isdir(dataExportFolder2) || mkpath(dataExportFolder2)
        end
        result = DataFrame()
        result.py = 0:cellNum-1
        result.ani = fill(subject, cellNum)
        result.unitId = fill(unitId, cellNum)
        result.planeId = fill(plId, cellNum)
        result.cellId = keys(ulsta)
        iscone=[];isachro=[];domicone=[];delycone=[];conemaxExt=[];conemaxMag=[];lcwmn=[];mcwmn=[];scwmn=[];lcwmg=[];mcwmg=[];scwmg=[];
        lcwmgall=[];mcwmgall=[];scwmgall=[];lmg=[];mmg=[];smg=[];amg=[];achResp=[];ls=[];ms=[];ss=[];as=[];isl=[];ism=[];iss=[];isa=[];
        ubcone=sta["ubcone"]
        cowg=sta["coneweight"]
        achroResp=sta["achroResp"]
        ustasign=sta["ustasign"]
        ustamag=sta["ustamag"]
        uconeresponsive = sta["uconeresponsive"]

        cellId = result.cellId
        for k=1:cellNum
        # k=1
            cone=haskey(ubcone,cellId[k]) ? true : false
            push!(isl,uconeresponsive[cellId[k]][1])
            push!(ism,uconeresponsive[cellId[k]][2])
            push!(iss,uconeresponsive[cellId[k]][3])
            push!(isa,uconeresponsive[cellId[k]][4])
            push!(iscone,cone)
            push!(isachro,in(cellId[k],achroResp[:cellId]))
            push!(domicone,map(i->cone==true ? ubcone[i].bstidx : NaN, cellId[k]))
            push!(delycone,map(i->cone==true ? ubcone[i].bstdly : NaN, cellId[k]))
            push!(conemaxExt,map(i->cone==true ? ubcone[i].bstex : NaN, cellId[k]))
            push!(conemaxMag,map(i->cone==true ? ubcone[i].bstmag : NaN, cellId[k]))
            push!(lcwmn,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Lconemn]) ? NaN : cowg[isequal(i).(cowg.cellId),:Lconemn], cellId[k])[1])
            push!(mcwmn,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Mconemn]) ? NaN : cowg[isequal(i).(cowg.cellId),:Mconemn], cellId[k])[1])
            push!(scwmn,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Sconemn]) ? NaN : cowg[isequal(i).(cowg.cellId),:Sconemn], cellId[k])[1])

            push!(lcwmg,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Lconemg]) ? NaN : cowg[isequal(i).(cowg.cellId),:Lconemg], cellId[k])[1])
            push!(mcwmg,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Mconemg]) ? NaN : cowg[isequal(i).(cowg.cellId),:Mconemg], cellId[k])[1])
            push!(scwmg,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Sconemg]) ? NaN : cowg[isequal(i).(cowg.cellId),:Sconemg], cellId[k])[1])

            push!(lcwmgall,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Lconemgall]) ? NaN : cowg[isequal(i).(cowg.cellId),:Lconemgall], cellId[k])[1])
            push!(mcwmgall,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Mconemgall]) ? NaN : cowg[isequal(i).(cowg.cellId),:Mconemgall], cellId[k])[1])
            push!(scwmgall,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Sconemgall]) ? NaN : cowg[isequal(i).(cowg.cellId),:Sconemgall], cellId[k])[1])

            push!(achResp,map(i->isempty(achroResp[isequal(i).(achroResp.cellId),:uachro]) ? NaN : achroResp[isequal(i).(achroResp.cellId),:uachro], cellId[k])[1])
            push!(ls, ustasign[cellId[k]][1])
            push!(ms, ustasign[cellId[k]][2])
            push!(ss, ustasign[cellId[k]][3])
            push!(as, ustasign[cellId[k]][4])
            push!(lmg, ustamag[cellId[k]][6,1])
            push!(mmg, ustamag[cellId[k]][6,2])
            push!(smg, ustamag[cellId[k]][6,3])
            push!(amg, ustamag[cellId[k]][6,4])
        end
        result.isl=isl
        result.ism=ism
        result.iss=iss
        result.isa=isa
        result.iscone=iscone
        result.isachro=isachro
        result.dominantcone=domicone
        result.conedelay=delycone
        result.ConemaxExt=conemaxExt
        result.ConemaxMag=conemaxMag
        result.lcwmn=lcwmn
        result.mcwmn=mcwmn
        result.scwmn=scwmn
        result.lcwmg=lcwmg
        result.mcwmg=mcwmg
        result.scwmg=scwmg
        result.lcwmgall=lcwmgall
        result.mcwmgall=mcwmgall
        result.scwmgall=scwmgall
        result.achroResp=achResp
        result.lsign=ls
        result.msign=ms
        result.ssign=ss
        result.asign=as
        result.lmg=lmg
        result.mmg=mmg
        result.smg=smg
        result.amg=amg
        append!(staData, result)

        append!(planeData, result)
        # if isequal(j,planeNum)
        save(joinpath(dataExportFolder2,join([subject,"_",unitId, "_thres$staThres", "_sta_dataset.jld2"])),"planeData",planeData)
        CSV.write(joinpath(dataExportFolder2,join([subject,"_",unitId,"_thres$staThres", "_sta_dataset.csv"])), planeData)
        # end
        summ=DataFrame(id=result.unitId[1], planeid=plId[1], staNum=cellNum, staCone=sum(result.iscone), staAchro=sum(result.isachro), meandelay=mean(result.conedelay[BitArray(result.iscone)]), stddelay=std(result.conedelay[BitArray(result.iscone)]))
        append!(staSum, summ)
    end
end

## Load Hartely (Fourier) data
if addFourier
    for i = 1:sessNum
        display("Plane Num: $i")
        local unitId
        j=1#map(j->isodd(i) ? 1 : 2, i)
        # unitId = hartleyExptId[i][5:7]
        unitId = recordSession[i]
        plId = recordPlane[j]
        dataFolder = joinpath(mainpath, join(["U",unitId]),"_Summary","DataExport")  # load ori data for cpi calculation
        dataFile=matchfile(Regex("$subject*_*$unitId*_$plId*_thres$staThres*_fourier_dataset.jld2"),dir=dataFolder,join=true)[1]
        display(dataFile)
        fdata = load(dataFile)["dataset"]
        cellNum = length(fdata["signif"])
        kern = sort!(OrderedDict(fdata["kern"]))
        if isequal(j,1)
            global planeData, dataExportFolder2
            planeData = DataFrame()
            dataExportFolder2 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
            isdir(dataExportFolder2) || mkpath(dataExportFolder2)
        end
        result = DataFrame()
        result.py = 0:cellNum-1
        result.ani = fill(subject, cellNum)
        result.unitId = fill(unitId, cellNum)
        result.planeId = fill(plId, cellNum)
        result.cellId = keys(kern)

        isl=[];ism=[];iss=[];isa=[];iscone=[];isachro=[];   # significance
        ldelta=[];mdelta=[];sdelta=[];adelta=[];conemaxDelta=[]; domicone=[];  # magnitude change
        lstd=[];mstd=[];sstd=[];astd=[];  # std/magnitude
        ltau=[];mtau=[];stau=[];atau=[];   # tau/delay
        llambda=[];mlambda=[];slambda=[];alambda=[];   # lambda
        lori=[];mori=[];sori=[];aori=[];   # orientation
        lsf=[];msf=[];ssf=[];asf=[];   # spatial frequency
        lct=[];mct=[];sct=[];act=[];   # cell type

        signif=fdata["signif"]
        tau=fdata["taumax"]
        delta=fdata["kdelta"]
        kstd=fdata["kstdmax"]
        lambda=fdata["slambda"]
        orimax = fdata["orimax"]
        sfmax = fdata["sfmax"]
        orimean = fdata["orimean"]
        sfmean = fdata["sfmean"]

        cellId = result.cellId
        for k=1:cellNum
        # k=1
            cone=haskey(ubcone,cellId[k]) ? true : false
            push!(isl,signif[cellId[k]][1])
            push!(ism,signif[cellId[k]][2])
            push!(iss,signif[cellId[k]][3])
            push!(isa,signif[cellId[k]][4])
            push!(iscone,cone)
            push!(isachro,in(cellId[k],achroResp[:cellId]))
            push!(domicone,map(i->cone==true ? ubcone[i].bstidx : NaN, cellId[k]))
            push!(delycone,map(i->cone==true ? ubcone[i].bstdly : NaN, cellId[k]))
            push!(conemaxExt,map(i->cone==true ? ubcone[i].bstex : NaN, cellId[k]))
            push!(conemaxMag,map(i->cone==true ? ubcone[i].bstmag : NaN, cellId[k]))
            push!(lcwmn,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Lconemn]) ? NaN : cowg[isequal(i).(cowg.cellId),:Lconemn], cellId[k])[1])
            push!(mcwmn,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Mconemn]) ? NaN : cowg[isequal(i).(cowg.cellId),:Mconemn], cellId[k])[1])
            push!(scwmn,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Sconemn]) ? NaN : cowg[isequal(i).(cowg.cellId),:Sconemn], cellId[k])[1])

            push!(lcwmg,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Lconemg]) ? NaN : cowg[isequal(i).(cowg.cellId),:Lconemg], cellId[k])[1])
            push!(mcwmg,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Mconemg]) ? NaN : cowg[isequal(i).(cowg.cellId),:Mconemg], cellId[k])[1])
            push!(scwmg,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Sconemg]) ? NaN : cowg[isequal(i).(cowg.cellId),:Sconemg], cellId[k])[1])

            push!(lcwmgall,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Lconemgall]) ? NaN : cowg[isequal(i).(cowg.cellId),:Lconemgall], cellId[k])[1])
            push!(mcwmgall,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Mconemgall]) ? NaN : cowg[isequal(i).(cowg.cellId),:Mconemgall], cellId[k])[1])
            push!(scwmgall,map(i->isempty(cowg[isequal(i).(cowg.cellId),:Sconemgall]) ? NaN : cowg[isequal(i).(cowg.cellId),:Sconemgall], cellId[k])[1])

            push!(achResp,map(i->isempty(achroResp[isequal(i).(achroResp.cellId),:uachro]) ? NaN : achroResp[isequal(i).(achroResp.cellId),:uachro], cellId[k])[1])
            push!(ls, ustasign[cellId[k]][1])
            push!(ms, ustasign[cellId[k]][2])
            push!(ss, ustasign[cellId[k]][3])
            push!(as, ustasign[cellId[k]][4])
            push!(ldelta, ustamag[cellId[k]][6,1])
            push!(mdelta, ustamag[cellId[k]][6,2])
            push!(sdelta, ustamag[cellId[k]][6,3])
            push!(adelta, ustamag[cellId[k]][6,4])
        end
        result.isl=isl
        result.ism=ism
        result.iss=iss
        result.isa=isa
        result.iscone=iscone
        result.isachro=isachro
        result.dominantcone=domicone
        result.conedelay=delycone
        result.ConemaxExt=conemaxExt
        result.ConemaxMag=conemaxMag
        result.lcwmn=lcwmn
        result.mcwmn=mcwmn
        result.scwmn=scwmn
        result.lcwmg=lcwmg
        result.mcwmg=mcwmg
        result.scwmg=scwmg
        result.lcwmgall=lcwmgall
        result.mcwmgall=mcwmgall
        result.scwmgall=scwmgall
        result.achroResp=achResp
        result.lsign=ls
        result.msign=ms
        result.ssign=ss
        result.asign=as
        result.mdelta=mdelta
        result.mdelta=mdelta
        result.sdelta=sdelta
        result.adelta=adelta
        append!(staData, result)

        append!(planeData, result)
        # if isequal(j,planeNum)
        save(joinpath(dataExportFolder2,join([subject,"_",unitId, "_thres$staThres", "_sta_dataset.jld2"])),"planeData",planeData)
        CSV.write(joinpath(dataExportFolder2,join([subject,"_",unitId,"_thres$staThres", "_sta_dataset.csv"])), planeData)
        # end
        summ=DataFrame(id=result.unitId[1], planeid=plId[1], staNum=cellNum, staCone=sum(result.iscone), staAchro=sum(result.isachro), meandelay=mean(result.conedelay[BitArray(result.iscone)]), stddelay=std(result.conedelay[BitArray(result.iscone)]))
        append!(staSum, summ)
    end
end

## Save results
siteId = join(["U",join(recordSession)])
save(joinpath(dataExportFolder1,join([subject,"_",siteId,"_roibkg.jld2"])), "roi",rois, "bkg",bkgs)
save(joinpath(dataExportFolder1,join([subject,"_",siteId,"_",paraName,"_hue.jld2"])), "hueData",hueData,"hueSum",hueSum)
CSV.write(joinpath(dataExportFolder1,join([subject,"_",siteId,"_hueData.csv"])), hueData)
CSV.write(joinpath(dataExportFolder1,join([subject,"_",siteId,"_",paraName,"_hueSum.csv"])), hueSum)
save(joinpath(dataExportFolder1,join([subject,"_",siteId,"_",paraName,"_ori.jld2"])), "oriData",oriData, "oriSum",oriSum)
CSV.write(joinpath(dataExportFolder1,join([subject,"_",siteId,"_oriData.csv"])), oriData)
CSV.write(joinpath(dataExportFolder1,join([subject,"_",siteId,"_",paraName,"_oriSum.csv"])), oriSum)
if addSTA
    save(joinpath(dataExportFolder1,join([subject,"_",siteId,"_thres$staThres","_sta.jld2"])), "staData",staData,"staSum",staSum)
    CSV.write(joinpath(dataExportFolder1,join([subject,"_",siteId,"_thres$staThres","_staData.csv"])), staData)
    CSV.write(joinpath(dataExportFolder1,join([subject,"_",siteId,"_thres$staThres","_staSum.csv"])), staSum)
end
## Data is organized accoording plane

for i = 1:exptOriNum
    # display("i: $i")
    j=1#map(j->isodd(i) ? 1 : 2, i)
    unitId = oriExptId[i][5:7]
    plId = recordPlane[j]
    if isequal(j,1)
        global planeSum, dataExportFolder3
        planeSum = DataFrame()
        dataExportFolder3 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
        isdir(dataExportFolder3) || mkpath(dataExportFolder3)
    end

    locId = join([unitId, "_", plId])
    roi=rois[i]
    xloc=[];yloc=[];
    for k = 1:size(roi,1)
        x = round(mean(roi[k][:,1]))
        y = round(mean(roi[k][:,2]))
        push!(xloc,x);push!(yloc,y)
    end
    pldata1 = filter(row -> row[:locId]==locId, oriData)
    pldata2 = filter(row -> row[:locId]==locId, hueData)
    insertcols!(pldata2, 6, :xloc=>xloc)
    insertcols!(pldata2, 7, :yloc=>yloc)
    pldata2.visResp = pldata1.visResp + pldata2.visResp
    pldata2.responsive = pldata1.responsive + pldata2.responsive
    pldata2.modulative = pldata1.modulative + pldata2.modulative
    pldata = hcat(pldata2, select(pldata1, Not([:py,:ani,:dataId,:locId,:cellId,:visResp,:responsive,:modulative])))
    planesum=DataFrame(session=unitId, planeid=plId, cellNum=size(pldata)[1], visResp=sum(pldata.visResp.>0),
        oriResp=sum((pldata.oriauc.>oriaucThres).&pldata.visResp), oriRespfit=sum((pldata.oriauc.>oriaucThres).&pldata.visResp.&(.~isnan.(pldata.fitori))), oriSelcv=sum((pldata.oricv.<oricvThres).&(pldata.oriauc.>oriaucThres).&pldata.visResp), oriSelosi=sum((pldata.osi.>osiThres).&(pldata.oriauc.>oriaucThres).&pldata.visResp),
        dirResp=sum((pldata.dirauc.>diraucThres).&pldata.visResp), dirRespfit=sum((pldata.dirauc.>diraucThres).&pldata.visResp.&(.~isnan.(pldata.fitdir))), dirSelcv=sum((pldata.dircv.<dircvThres).&(pldata.dirauc.>diraucThres).&pldata.visResp), dirSeldsi=sum((pldata.dsi.>dsiThres).&(pldata.dirauc.>diraucThres).&pldata.visResp),
        hueaxResp=sum((pldata.hueaxauc.>hueaucThres).&pldata.visResp), hueaxRespfit=sum((pldata.hueaxauc.>hueaucThres).&pldata.visResp.&(.~isnan.(pldata.fithueax))),hueaxSelcv=sum((pldata.hueaxcv.<huecvThres).&(pldata.hueaxauc.>hueaucThres).&pldata.visResp), hueaxSelsi=sum((pldata.hueaxsi.>huesiThres).&(pldata.hueaxauc.>hueaucThres).&pldata.visResp),
        huediResp=sum((pldata.huediauc.>hueaucThres).&pldata.visResp), huediRespfit=sum((pldata.huediauc.>hueaucThres).&pldata.visResp.&(.~isnan.(pldata.fithuedi))),huediSelcv=sum((pldata.huedicv.<huecvThres).&(pldata.huediauc.>hueaucThres).&pldata.visResp), huediSelsi=sum((pldata.huedisi.>huesiThres).&(pldata.huediauc.>hueaucThres).&pldata.visResp),
        hueoriResp=sum((pldata.hueoriauc.>oriaucThres).&pldata.visResp), hueoriRespfit=sum((pldata.hueoriauc.>oriaucThres).&pldata.visResp.&(.~isnan.(pldata.fithueori))),hueoriSelcv=sum((pldata.hueoricv.<oricvThres).&(pldata.hueoriauc.>oriaucThres).&pldata.visResp), hueoriSeldsi=sum((pldata.hueosi.>osiThres).&(pldata.hueoriauc.>oriaucThres).&pldata.visResp),
        huedirResp=sum((pldata.huedirauc.>diraucThres).&pldata.visResp), huedirRespfit=sum((pldata.huedirauc.>diraucThres).&pldata.visResp.&(.~isnan.(pldata.fithuedir))),huedirSelcv=sum((pldata.huedircv.<dircvThres).&(pldata.huedirauc.>diraucThres).&pldata.visResp), huedirSeldsi=sum((pldata.huedsi.>dsiThres).&(pldata.huedirauc.>diraucThres).&pldata.visResp),
        hueaucResp=sum(((pldata.hueaxauc.>hueaucThres).|(pldata.huediauc.>hueaucThres)).&pldata.visResp), hueSelcpi=sum(((pldata.hueaxauc.>hueaucThres).|(pldata.huediauc.>hueaucThres)).&pldata.visResp.&(pldata.huecpi.>cpiThres)))
    append!(planeSum, planesum)
    save(joinpath(dataExportFolder3,join([subject,"_",locId,"_",paraName,"_sum.jld2"])), locId,pldata,string(locId,"_sum"),planesum)
    CSV.write(joinpath(dataExportFolder3,join([subject,"_",locId,"_sum.csv"])), pldata)

    if isequal(j, planeNum)
        save(joinpath(dataExportFolder3,join([subject,"_",unitId,"_",paraName,"_sum.jld2"])), "PlaneSum", planeSum)
        CSV.write(joinpath(dataExportFolder3,join([subject,"_",unitId,"_",paraName,"_sum.csv"])), planeSum)
    end
end
display("Processing is Done!!!!")
# exptNum = exptOriNum + exptHueNum
# CSV.write(joinpath(dataExportFolder1,join([subject,"_hue.csv"])), hueData)
