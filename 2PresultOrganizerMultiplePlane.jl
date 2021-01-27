# PL: This code is to organize 2P data collected from one location (maybe multiple planes) for following plotting code
# Remember: fitting is done when auc>0.5
#

using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact, CSV,MAT,Query,DataStructures, HypothesisTests, StatsFuns, Random

# Expt info
disk = "O:"
subject = "AF4"  # Animal
# recordSession = ["002","003","004", "005", "006", "007", "008"] # Unit
recordPlane = ["000", "001"]

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
colorSpace = "DKL"   # DKL, # HSL

addSTA = false
addFourier = true
staThres = 0.25

## Make folders and path
mainpath = joinpath(disk,subject, "2P_analysis")
dataExportFolder1 = joinpath(mainpath, "Summary", "DataExport")
isdir(dataExportFolder1) || mkpath(dataExportFolder1)

meta = readmeta(joinpath(dataExportFolder1,"metadata.mat"))
for t in 1:nrow(meta)
    meta.files[t]=meta.files[t][3:end]
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
    @where i.RecordSite == "u004"
    # @where i.RecordSite != "u008"
    @where i.ID == "DirSF"
    @select {i.ID,i.RecordSite,i.filename}
    @collect DataFrame
    end
# deleterows!(oritests, [2,3,4,5,6,8,9,11,12,13,14,15,16,17,19,20,21,22,24,25,26,27])
# deleterows!(oritests, [2,4,5,8,10])
deleterows!(oritests, [2,3,4,5,6])
deleterows!(oritests, [2])
huetests = @from i in tests begin
    # @where i.RecordSite != "u005"
    @where i.RecordSite == "u004"
    @where i.ID == "DirSFColor"
    @select {i.ID,i.RecordSite,i.filename}
    @collect DataFrame
    end
# deleterows!(huetests, [2,3,6,7,9])
deleterows!(huetests, [2])

hartleytests = @from i in tests begin
    # @where i.RecordSite != "u005"
    @where i.RecordSite == "u004"
    @where i.ID == "Hartley"
    @select {i.ID,i.RecordSite,i.filename}
    @collect DataFrame
    end

##
recordSession = unique!(oritests.RecordSite) # Unit
recordSession = map(i->i[2:4],recordSession)
oriExptId = oritests.filename
hueExptId = huetests.filename # Stimulus test
# hartleyExptId = hartleytests.filename    # Unit Ids have L, M, S, and Achromatic hartleys.
exptOriNum = length(oriExptId)
exptHueNum = length(hueExptId)
planeNum = length(recordPlane)
sessNum = length(recordSession)

paraName = string("oa",string(oriaucThres),"da",string(diraucThres),"ha",string(hueaucThres),"oc",string(oricvThres),
           "dc",string(dsiThres),"hc",string(huecvThres),"osi",string(osiThres),"dsi",string(dsiThres),"hsi",string(huesiThres),"cpi",string(cpiThres))

## Direction/Orientation data
rois=Any[]; bkgs=Any[]; oriData=DataFrame(); oriSum=DataFrame(); hueData=DataFrame(); hueSum=DataFrame();
staData=DataFrame(); staSum=DataFrame();fourierData=DataFrame(); fourierSum=DataFrame();
for i = 1:exptOriNum
    display("Processing Oriexpt No: $i")
    for j = 1:planeNum
        local unitId
        local siteId
        local result
        unitId = oriExptId[i][5:7]
        siteId = join([unitId, "_",recordPlane[j]])
        dataFolder = joinpath(mainpath, join(["U", unitId]), join([oriExptId[i][5:11], "_", recordPlane[j]]), "DataExport")
        dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,join=true)[1]
        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
        alignFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.align"),dir=dataFolder,join=true)[1]

        # segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
        # alignFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*.align"),dir=dataFolder,join=true)[1]

        result = load(dataFile)["result"]
        segment = prepare(segmentFile)
        align = prepare(alignFile)
        cellRoi = segment["vert"]
        bkgImg = align["m"]
        locID = repeat([siteId],length(cellRoi))
        insertcols!(result,4,:locId => locID)
        append!(oriData, result)
        push!(rois, cellRoi)
        push!(bkgs, bkgImg)
        summ=DataFrame(id=result.dataId[1], exptid=oriExptId[i], planeid=recordPlane[j], cellNum=length(cellRoi), visResp=sum(result.visResp),
             oriResp=sum((result.oriauc.>oriaucThres).&result.visResp), oriRespfit=sum((result.oriauc.>oriaucThres).&result.visResp.&(.~isnan.(result.fitori))), oriSelcv=sum((result.oricv.<oricvThres).&(result.oriauc.>oriaucThres).&result.visResp), oriSelosi=sum((result.osi1.>osiThres).&(result.oriauc.>oriaucThres).&result.visResp),
             dirResp=sum((result.dirauc.>diraucThres).&result.visResp), dirRespfit=sum((result.dirauc.>diraucThres).&result.visResp.&(.~isnan.(result.fitdir))), dirSelcv=sum((result.dircv.<dircvThres).&(result.dirauc.>diraucThres).&result.visResp), dirSeldsi=sum((result.dsi1.>dsiThres).&(result.dirauc.>diraucThres).&result.visResp))
        append!(oriSum, summ)
        dataExportFolder2 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
        isdir(dataExportFolder2) || mkpath(dataExportFolder2)
        # save(joinpath(dataExportFolder2,join([subject,"_",siteId,"_roibkg.jld2"])), Dict("roi"=>cellRoi, "bkg"=>bkgImg))   # only need save one for each plane in each unit
        save(joinpath(dataExportFolder2,join([subject,"_",siteId,"_roibkg.jld2"])), "roi",cellRoi, "bkg",bkgImg)
    end
end

## Hue data
for i = 1:exptHueNum
    display("Processing Hueexpt No: $i")
    for j = 1:planeNum
        local siteId
        local result
        siteId = join([recordSession[i], "_",recordPlane[j]])
        dataFolder = joinpath(mainpath, join(["U", recordSession[i]]), join([oriExptId[i][5:11], "_", recordPlane[j]]), "DataExport")
        dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,join=true)[1]
        resultori = load(dataFile)["result"]
        dataFolder = joinpath(mainpath, join(["U", recordSession[i]]), join([hueExptId[i][5:11], "_", recordPlane[j]]), "DataExport")
        dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,join=true)[1]
        result = load(dataFile)["result"]
        locID = repeat([siteId],size(result)[1])
        insertcols!(result,4,locId=locID)
        hueResp = 1 .- minimum(hcat(result.hueaxcv, result.huedicv),dims=2)
        achroResp = 1 .- resultori.oricv
        cpi = dropdims((hueResp .- achroResp) ./ (hueResp .+ achroResp),dims=2)
        insertcols!(result,11,huecpi=cpi)
        append!(hueData, result)
        summ=DataFrame(id=result.dataId[1], exptid=hueExptId[i], planeid=recordPlane[j], cellNum=size(result)[1], visResp=sum(result.visResp),
             hueaxResp=sum((result.hueaxauc.>hueaucThres).&result.visResp), hueaxRespfit=sum((result.hueaxauc.>hueaucThres).&result.visResp.&(.~isnan.(result.fithueax))),hueaxSelcv=sum((result.hueaxcv.<huecvThres).&(result.hueaxauc.>hueaucThres).&result.visResp), hueaxSelsi=sum((result.hueaxsi1.>huesiThres).&(result.hueaxauc.>hueaucThres).&result.visResp),
             huediResp=sum((result.huediauc.>hueaucThres).&result.visResp), huediRespfit=sum((result.huediauc.>hueaucThres).&result.visResp.&(.~isnan.(result.fithuedi))),huediSelcv=sum((result.huedicv.<huecvThres).&(result.huediauc.>hueaucThres).&result.visResp), huediSelsi=sum((result.huedisi1.>huesiThres).&(result.huediauc.>hueaucThres).&result.visResp),
             hueoriResp=sum((result.hueoriauc.>oriaucThres).&result.visResp), hueoriRespfit=sum((result.hueoriauc.>oriaucThres).&result.visResp.&(.~isnan.(result.fithueori))),hueoriSelcv=sum((result.hueoricv.<oricvThres).&(result.hueoriauc.>oriaucThres).&result.visResp), hueoriSeldsi=sum((result.hueosi1.>osiThres).&(result.hueoriauc.>oriaucThres).&result.visResp),
             huedirResp=sum((result.huedirauc.>diraucThres).&result.visResp), huedirRespfit=sum((result.huedirauc.>diraucThres).&result.visResp.&(.~isnan.(result.fithuedir))),huedirSelcv=sum((result.huedircv.<dircvThres).&(result.huedirauc.>diraucThres).&result.visResp), huedirSeldsi=sum((result.huedsi1.>dsiThres).&(result.huedirauc.>diraucThres).&result.visResp),
             hueaucResp=sum(((result.hueaxauc.>hueaucThres).|(result.huediauc.>hueaucThres)).&result.visResp))#, hueSelcpi=sum(((result.hueaxauc.>hueaucThres).|(result.huediauc.>hueaucThres)).&result.visResp.&(cpi.>cpiThres)))
        append!(hueSum, summ)
    end
end

## Load Hartely (STA) data
if addSTA
    for i = 1:sessNum
        display("RecordSession Num: $i")
        for j = 1:planeNum
            local unitId
            local dataFolder
            unitId = recordSession[i]
            plId = recordPlane[j]
            dataFolder = joinpath(mainpath, join(["U",unitId]),"_Summary","DataExport")  # load ori data for cpi calculation
            dataFile=matchfile(Regex("$subject*_*$unitId*_$plId*_thres$staThres*_sta_datasetFinal.jld2"),dir=dataFolder,join=true)[1]
            sta = load(dataFile)["datasetFinal"]
            cellNum = length(sta["ulsta"])
            ulsta = sort!(OrderedDict(sta["ulsta"]))
            if isequal(j,1)
                global planeData, dataExportFolder2
                planeData = DataFrame()
                dataExportFolder2 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
                isdir(dataExportFolder2) || mkpath(dataExportFolder2)
            end
            local result
            result = DataFrame()
            result.py = 0:cellNum-1
            result.ani = fill(subject, cellNum)
            result.unitId = fill(unitId, cellNum)
            result.planeId = fill(plId, cellNum)
            result.cellId = collect(keys(ulsta))
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
                push!(iscone,cone)
                push!(isl,uconeresponsive[cellId[k]][1])
                push!(ism,uconeresponsive[cellId[k]][2])
                push!(iss,uconeresponsive[cellId[k]][3])
                push!(isa,uconeresponsive[cellId[k]][4])
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
            if isequal(j,planeNum)
                save(joinpath(dataExportFolder2,join([subject,"_",unitId, "_thres$staThres", "_sta_dataset.jld2"])),"planeData",planeData)
                CSV.write(joinpath(dataExportFolder2,join([subject,"_",unitId,"_thres$staThres", "_sta_dataset.csv"])), planeData)
            end
            summ=DataFrame(id=result.unitId[1], planeid=plId[1], staNum=cellNum, staCone=sum(result.iscone), staAchro=sum(result.isachro), meandelay=mean(delycone[BitArray(iscone)]), stddelay=std(delycone[BitArray(iscone)]))
            append!(staSum, summ)
        end
    end
end

## Load Hartely (Fourier) data
if addFourier
    for i = 1:sessNum
        display("RecordSession Num: $i")
        for j = 1:planeNum

            local unitId
            local dataFolder
            unitId = recordSession[i]
            plId = recordPlane[j]
            display("Plane Id.: $plId")
            dataFolder = joinpath(mainpath, join(["U",unitId]),"_Summary","DataExport")  # load ori data for cpi calculation
            dataFile=matchfile(Regex("$subject*_*$unitId*_$plId*_thres$staThres*_fourier_dataset.jld2"),dir=dataFolder,join=true)[1]
            fdata = load(dataFile)["dataset"]
            cellNum = length(fdata["signif"])
            kern = sort!(OrderedDict(fdata["kern"]))
            if isequal(j,1)
                global planeData, dataExportFolder2
                planeData = DataFrame()
                dataExportFolder2 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
                isdir(dataExportFolder2) || mkpath(dataExportFolder2)
            end
            local result
            result = DataFrame()
            result.py = 0:cellNum-1
            result.ani = fill(subject, cellNum)
            result.unitId = fill(unitId, cellNum)
            result.planeId = fill(plId, cellNum)
            result.cellId = collect(keys(kern))

            isl=[];ism=[];iss=[];isa=[];iscone=[];  # significance
            ldelta=[];mdelta=[];sdelta=[];adelta=[];conemaxDelta=[]; domicone=[];  # magnitude change
            lstd=[];mstd=[];sstd=[];astd=[];  # std/magnitude
            ltau=[];mtau=[];stau=[];atau=[];   # tau/delay
            llambda=[];mlambda=[];slambda=[];alambda=[];   # lambda
            lorimax=[];morimax=[];sorimax=[];aorimax=[];   # orientation based on maximal value
            lorimean=[];morimean=[];sorimean=[];aorimean=[];   # orientation based on mean value
            lpo=[];mpo=[];spo=[];apo=[];   # Preferred orientation from fitting
            lohw=[];mohw=[];sohw=[];aohw=[];   # half-width of orientaiton tuning.
            losi1=[];mosi1=[];sosi1=[];aosi1=[];   # osi of orientaiton tuning.
            losi2=[];mosi2=[];sosi2=[];aosi2=[];   # osi of orientaiton tuning.
            lsfmax=[];msfmax=[];ssfmax=[];asfmax=[];   # spatial frequency based maximal value
            lsfmean=[];msfmean=[];ssfmean=[];asfmean=[];   # spatial frequency based mean value
            lpsf=[];mpsf=[];spsf=[];apsf=[];   # Preferred SF from fitting
            lsfhw=[];msfhw=[];ssfhw=[];asfhw=[];   # half-width of sf tuning.
            lsfbw=[];msfbw=[];ssfbw=[];asfbw=[];   # half-width of sf tuning.
            lsft=[];msft=[];ssft=[];asft=[];   # SF type

            signif=fdata["signif"]
            tau=fdata["taumax"]
            delta=fdata["kdelta"]
            kstd=fdata["kstdmax"]
            lambda=fdata["slambda"]
            orimax = fdata["orimax"]
            sfmax = fdata["sfmax"]
            orimean = fdata["orimean"]
            sfmean = fdata["sfmean"]
            sffit = fdata["sffit"]
            orifit = fdata["orifit"]

            cellId = result.cellId
            for k=1:cellNum
                # k=1
                # print(k,"-")
                cone=isequal(sum(signif[cellId[k]][1:3]),0) ? false : true
                push!(iscone,cone)
                push!(isl,signif[cellId[k]][1])
                push!(ism,signif[cellId[k]][2])
                push!(iss,signif[cellId[k]][3])
                push!(isa,signif[cellId[k]][4])

                push!(ldelta,delta[cellId[k]][1])
                push!(mdelta,delta[cellId[k]][2])
                push!(sdelta,delta[cellId[k]][3])
                push!(adelta,delta[cellId[k]][4])
                # push!(domicone,map(i->cone==true ? ubcone[i].bstidx : NaN, cellId[k]))

                push!(lstd,kstd[cellId[k]][1])
                push!(mstd,kstd[cellId[k]][2])
                push!(sstd,kstd[cellId[k]][3])
                push!(astd,kstd[cellId[k]][4])

                push!(ltau,tau[cellId[k]][1])
                push!(mtau,tau[cellId[k]][2])
                push!(stau,tau[cellId[k]][3])
                push!(atau,tau[cellId[k]][4])

                push!(llambda,lambda[cellId[k]][1])
                push!(mlambda,lambda[cellId[k]][2])
                push!(slambda,lambda[cellId[k]][3])
                push!(alambda,lambda[cellId[k]][4])

                push!(lorimax,orimax[cellId[k]][1])
                push!(morimax,orimax[cellId[k]][2])
                push!(sorimax,orimax[cellId[k]][3])
                push!(aorimax,orimax[cellId[k]][4])

                push!(lorimean,orimean[cellId[k]][1])
                push!(morimean,orimean[cellId[k]][2])
                push!(sorimean,orimean[cellId[k]][3])
                push!(aorimean,orimean[cellId[k]][4])

                push!(lpo,orifit[cellId[k]][1].po)
                push!(mpo,orifit[cellId[k]][2].po)
                push!(spo,orifit[cellId[k]][3].po)
                push!(apo,orifit[cellId[k]][4].po)

                push!(lohw,orifit[cellId[k]][1].ohw)
                push!(mohw,orifit[cellId[k]][2].ohw)
                push!(sohw,orifit[cellId[k]][3].ohw)
                push!(aohw,orifit[cellId[k]][4].ohw)

                push!(losi1,orifit[cellId[k]][1].osi1)
                push!(mosi1,orifit[cellId[k]][2].osi1)
                push!(sosi1,orifit[cellId[k]][3].osi1)
                push!(aosi1,orifit[cellId[k]][4].osi1)

                push!(losi2,orifit[cellId[k]][1].osi2)
                push!(mosi2,orifit[cellId[k]][2].osi2)
                push!(sosi2,orifit[cellId[k]][3].osi2)
                push!(aosi2,orifit[cellId[k]][4].osi2)

                push!(lsfmax,sfmax[cellId[k]][1])
                push!(msfmax,sfmax[cellId[k]][2])
                push!(ssfmax,sfmax[cellId[k]][3])
                push!(asfmax,sfmax[cellId[k]][4])

                push!(lsfmean,sfmean[cellId[k]][1])
                push!(msfmean,sfmean[cellId[k]][2])
                push!(ssfmean,sfmean[cellId[k]][3])
                push!(asfmean,sfmean[cellId[k]][4])

                push!(lpsf,sffit[cellId[k]][1].psf)
                push!(mpsf,sffit[cellId[k]][2].psf)
                push!(spsf,sffit[cellId[k]][3].psf)
                push!(apsf,sffit[cellId[k]][4].psf)

                push!(lsfhw,sffit[cellId[k]][1].sfhw)
                push!(msfhw,sffit[cellId[k]][2].sfhw)
                push!(ssfhw,sffit[cellId[k]][3].sfhw)
                push!(asfhw,sffit[cellId[k]][4].sfhw)

                push!(lsfbw,sffit[cellId[k]][1].sfbw)
                push!(msfbw,sffit[cellId[k]][2].sfbw)
                push!(ssfbw,sffit[cellId[k]][3].sfbw)
                push!(asfbw,sffit[cellId[k]][4].sfbw)

                push!(lsft,sffit[cellId[k]][1].sftype)
                push!(msft,sffit[cellId[k]][2].sftype)
                push!(ssft,sffit[cellId[k]][3].sftype)
                push!(asft,sffit[cellId[k]][4].sftype)

                # push!(domicone,map(i->cone==true ? ubcone[i].bstidx : NaN, cellId[k]))
                # push!(delycone,map(i->cone==true ? ubcone[i].bstdly : NaN, cellId[k]))
                # push!(conemaxExt,map(i->cone==true ? ubcone[i].bstex : NaN, cellId[k]))
                # push!(conemaxMag,map(i->cone==true ? ubcone[i].bstmag : NaN, cellId[k]))

            end
            result.isl=isl
            result.ism=ism
            result.iss=iss
            result.isa=isa
            result.iscone=iscone

            # result.dominantcone=domicone
            # result.conedelay=delycone
            # result.ConemaxExt=conemaxExt
            # result.ConemaxMag=conemaxMag
            result.ldelta=ldelta
            result.mdelta=mdelta
            result.sdelta=sdelta
            result.adelta=adelta

            result.lstd=lstd
            result.mstd=mstd
            result.sstd=sstd
            result.astd=astd

            result.ltau=ltau
            result.mtau=mtau
            result.stau=stau
            result.atau=atau

            result.llambda=llambda
            result.mlambda=mlambda
            result.slambda=slambda
            result.alambda=alambda

            result.lorimax=lorimax
            result.morimax=morimax
            result.sorimax=sorimax
            result.aorimax=aorimax

            result.lorimean=lorimean
            result.morimean=morimean
            result.sorimean=sorimean
            result.aorimean=aorimean

            result.lpo=lpo
            result.mpo=mpo
            result.spo=spo
            result.apo=apo

            result.lohw=lohw
            result.mohw=mohw
            result.sohw=sohw
            result.aohw=aohw

            result.losi1=losi1
            result.mosi1=mosi1
            result.sosi1=sosi1
            result.aosi1=aosi1

            result.losi2=losi2
            result.mosi2=mosi2
            result.sosi2=sosi2
            result.aosi2=aosi2

            result.lsfmax=lsfmax
            result.msfmax=msfmax
            result.ssfmax=ssfmax
            result.asfmax=asfmax

            result.lsfmean=lsfmean
            result.msfmean=msfmean
            result.ssfmean=ssfmean
            result.asfmean=asfmean

            result.lpsf=lpsf
            result.mpsf=mpsf
            result.spsf=spsf
            result.apsf=apsf

            result.lsfhw=lsfhw
            result.msfhw=msfhw
            result.ssfhw=ssfhw
            result.asfhw=asfhw

            result.lsfbw=lsfbw
            result.msfbw=msfbw
            result.ssfbw=ssfbw
            result.asfbw=asfbw

            result.lsft=lsft
            result.msft=msft
            result.ssft=ssft
            result.asft=asft

            append!(fourierData, result)

            append!(planeData, result)
            if isequal(j,planeNum)
                save(joinpath(dataExportFolder2,join([subject,"_",unitId, "_thres$staThres", "_fourier_dataset_organized.jld2"])),"planeData",planeData)
                CSV.write(joinpath(dataExportFolder2,join([subject,"_",unitId,"_thres$staThres", "_fourier_dataset_organized.csv"])), planeData)
            end
            # summ=DataFrame(id=result.unitId[1], planeid=plId[1], fourierNum=cellNum, staCone=sum(result.iscone), staAchro=sum(result.isachro), meandelay=mean(delycone[BitArray(iscone)]), stddelay=std(delycone[BitArray(iscone)]))
            # append!(fourierSum, summ)
        end
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
ct=0
for i = 1:sessNum
    local unitId
    unitId = recordSession[i]
    planeSum = DataFrame()
    dataExportFolder3 = joinpath(mainpath, join(["U", unitId]), "_Summary", "DataExport")
    isdir(dataExportFolder3) || mkpath(dataExportFolder3)
        for j = 1:planeNum
            plId = recordPlane[j]
            locId = join([unitId, "_", plId])

            global ct
            ct=ct+1
            roi=rois[ct]

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
            planesum=DataFrame(session=recordSession[i], planeid=recordPlane[j], cellNum=size(pldata)[1], visResp=sum(pldata.visResp.>0),
                oriResp=sum((pldata.oriauc.>oriaucThres).&pldata.visResp), oriRespfit=sum((pldata.oriauc.>oriaucThres).&pldata.visResp.&(.~isnan.(pldata.fitori))), oriSelcv=sum((pldata.oricv.<oricvThres).&(pldata.oriauc.>oriaucThres).&pldata.visResp), oriSelosi=sum((pldata.osi1.>osiThres).&(pldata.oriauc.>oriaucThres).&pldata.visResp),
                dirResp=sum((pldata.dirauc.>diraucThres).&pldata.visResp), dirRespfit=sum((pldata.dirauc.>diraucThres).&pldata.visResp.&(.~isnan.(pldata.fitdir))), dirSelcv=sum((pldata.dircv.<dircvThres).&(pldata.dirauc.>diraucThres).&pldata.visResp), dirSeldsi=sum((pldata.dsi1.>dsiThres).&(pldata.dirauc.>diraucThres).&pldata.visResp),
                hueaxResp=sum((pldata.hueaxauc.>hueaucThres).&pldata.visResp), hueaxRespfit=sum((pldata.hueaxauc.>hueaucThres).&pldata.visResp.&(.~isnan.(pldata.fithueax))),hueaxSelcv=sum((pldata.hueaxcv.<huecvThres).&(pldata.hueaxauc.>hueaucThres).&pldata.visResp), hueaxSelsi=sum((pldata.hueaxsi1.>huesiThres).&(pldata.hueaxauc.>hueaucThres).&pldata.visResp),
                huediResp=sum((pldata.huediauc.>hueaucThres).&pldata.visResp), huediRespfit=sum((pldata.huediauc.>hueaucThres).&pldata.visResp.&(.~isnan.(pldata.fithuedi))),huediSelcv=sum((pldata.huedicv.<huecvThres).&(pldata.huediauc.>hueaucThres).&pldata.visResp), huediSelsi=sum((pldata.huedisi1.>huesiThres).&(pldata.huediauc.>hueaucThres).&pldata.visResp),
                hueoriResp=sum((pldata.hueoriauc.>oriaucThres).&pldata.visResp), hueoriRespfit=sum((pldata.hueoriauc.>oriaucThres).&pldata.visResp.&(.~isnan.(pldata.fithueori))),hueoriSelcv=sum((pldata.hueoricv.<oricvThres).&(pldata.hueoriauc.>oriaucThres).&pldata.visResp), hueoriSeldsi=sum((pldata.hueosi1.>osiThres).&(pldata.hueoriauc.>oriaucThres).&pldata.visResp),
                huedirResp=sum((pldata.huedirauc.>diraucThres).&pldata.visResp), huedirRespfit=sum((pldata.huedirauc.>diraucThres).&pldata.visResp.&(.~isnan.(pldata.fithuedir))),huedirSelcv=sum((pldata.huedircv.<dircvThres).&(pldata.huedirauc.>diraucThres).&pldata.visResp), huedirSeldsi=sum((pldata.huedsi1.>dsiThres).&(pldata.huedirauc.>diraucThres).&pldata.visResp),
                hueaucResp=sum(((pldata.hueaxauc.>hueaucThres).|(pldata.huediauc.>hueaucThres)).&pldata.visResp), hueSelcpi=sum(((pldata.hueaxauc.>hueaucThres).|(pldata.huediauc.>hueaucThres)).&pldata.visResp.&(pldata.huecpi.>cpiThres)))
            append!(planeSum, planesum)
            save(joinpath(dataExportFolder3,join([subject,"_",locId,"_",paraName,"_sum.jld2"])), locId,pldata,string(locId,"_sum"),planesum)
            CSV.write(joinpath(dataExportFolder3,join([subject,"_",locId,"_sum.csv"])), pldata)
        end
    save(joinpath(dataExportFolder3,join([subject,"_",unitId,"_",paraName,"_sum.jld2"])), "PlaneSum", planeSum)
    CSV.write(joinpath(dataExportFolder3,join([subject,"_",unitId,"_",paraName,"_sum.csv"])), planeSum)
end
display("Processing is Done!!!!")
