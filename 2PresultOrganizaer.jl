
# Rember fitting is done when auc>0.5

using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact, CSV,MAT, DataStructures, HypothesisTests, StatsFuns, Random

# Expt info
disk = "H:"
subject = "AE7"  # Animal
# recordSession = ["002","003","004", "005", "006","007"] # Unit
# recordPlane = ["000", "001"]
# oriExptId = ["002_000","003_008", "004_005", "005_010","006_007","007_009"]
# hueExptId = ["002_006","003_012", "004_008", "005_007","006_005","007_006"]  # Stimulus test

recordSession = ["001"]#,"003","004", "005", "006","007"] # Unit
recordPlane = ["000"]
oriExptId = ["001_001"]
hueExptId = []  #

oriaucThres = 0.7
diraucThres = 0.7
hueaucThres = 0.9

oricvThres = 0.7
dircvThres = 0.7
huecvThres = 0.7

osiThres = 0.5
dsiThres = 0.5
huesiThres = 0.5
cpiThres = 0.33

numberofColor = 12
colorSpace = "DKL"

mainpath = joinpath(disk,subject, "2P_analysis")
dataExportFolder1 = joinpath(mainpath, "Summary", "DataExport")
isdir(dataExportFolder1) || mkpath(dataExportFolder1)
paraName = string("oa",string(oriaucThres),"da",string(diraucThres),"ha",string(hueaucThres),"oc",string(oricvThres),
           "dc",string(dsiThres),"hc",string(huecvThres),"osi",string(osiThres),"dsi",string(dsiThres),"hsi",string(huesiThres),"cpi",string(cpiThres))

exptOriNum = length(oriExptId)
exptHueNum = length(hueExptId)
planeNum = length(recordPlane)

## Direction/Orientation data
rois = Any[]; bkgs = Any[]; oriData = DataFrame(); oriSum =DataFrame(); hueData = DataFrame(); hueSum = DataFrame();

for i = 1:exptOriNum
    for j = 1:planeNum
        siteId = join([oriExptId[i][1:3], "_",recordPlane[j]])
        dataFolder = joinpath(mainpath, join(["U", oriExptId[i][1:3]]), join([oriExptId[i], "_", recordPlane[j]]), "DataExport")
        dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,adddir=true)[1]
        # segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.segment"),dir=dataFolder,adddir=true)[1]
        # alignFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.align"),dir=dataFolder,adddir=true)[1]

        segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*.segment"),dir=dataFolder,adddir=true)[1]
        alignFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*.align"),dir=dataFolder,adddir=true)[1]

        result = load(dataFile)["result"]
        segment = prepare(segmentFile)
        align = prepare(alignFile)
        cellRoi = segment["vert"]
        bkgImg = align["m"]
        locID = repeat([join([oriExptId[i][1:4], recordPlane[j]])],length(cellRoi))
        insertcols!(result,4,locId=locID)
        append!(oriData, result)
        push!(rois, cellRoi)
        push!(bkgs, bkgImg)
        summ=DataFrame(id=result.dataId[1], exptid=oriExptId[i], planeid=recordPlane[j], cellNum=length(cellRoi), visResp=sum(result.visResp),
             oriResp=sum((result.oriauc.>oriaucThres).&result.visResp), oriRespfit=sum((result.oriauc.>oriaucThres).&result.visResp.&(.~isnan.(result.fitori))), oriSelcv=sum((result.oricv.<oricvThres).&(result.oriauc.>oriaucThres).&result.visResp), oriSelosi=sum((result.osi.>osiThres).&(result.oriauc.>oriaucThres).&result.visResp),
             dirResp=sum((result.dirauc.>diraucThres).&result.visResp), dirRespfit=sum((result.dirauc.>diraucThres).&result.visResp.&(.~isnan.(result.fitdir))), dirSelcv=sum((result.dircv.<dircvThres).&(result.dirauc.>diraucThres).&result.visResp), dirSeldsi=sum((result.dsi.>dsiThres).&(result.dirauc.>diraucThres).&result.visResp))
        append!(oriSum, summ)
        dataExportFolder2 = joinpath(mainpath, join(["U", oriExptId[i][1:3]]), "_Summary", "DataExport")
        isdir(dataExportFolder2) || mkpath(dataExportFolder2)
        # save(joinpath(dataExportFolder2,join([subject,"_",siteId,"_roibkg.jld2"])), Dict("roi"=>cellRoi, "bkg"=>bkgImg))   # only need save one for each plane in each unit
        save(joinpath(dataExportFolder2,join([subject,"_",siteId,"_roibkg.jld2"])), "roi",cellRoi, "bkg",bkgImg)
    end
end

## Hue data
for i = 1:exptHueNum
    for j = 1:planeNum
        dataFolder = joinpath(mainpath, join(["U", oriExptId[i][1:3]]), join([oriExptId[i], "_", recordPlane[j]]), "DataExport")
        dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,adddir=true)[1]
        resultori = load(dataFile)["result"]
        dataFolder = joinpath(mainpath, join(["U", hueExptId[i][1:3]]), join([hueExptId[i], "_", recordPlane[j]]), "DataExport")
        dataFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_result.jld2"),dir=dataFolder,adddir=true)[1]
        result = load(dataFile)["result"]
        locID = repeat([join([hueExptId[i][1:4], recordPlane[j]])],size(result)[1])
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

## Data is organized accoording plane
ct=0
for i = 1:exptOriNum
    # display("i: $i")
    planeSum = DataFrame()
    dataExportFolder2 = joinpath(mainpath, join(["U", oriExptId[i][1:3]]), "_Summary", "DataExport")
    isdir(dataExportFolder2) || mkpath(dataExportFolder2)
    for j = 1:planeNum
        # display("j:$j")
        global ct
        ct = ct+1
        locId = join([oriExptId[i][1:4], recordPlane[j]])
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
            oriResp=sum((pldata.oriauc.>oriaucThres).&pldata.visResp), oriRespfit=sum((pldata.oriauc.>oriaucThres).&pldata.visResp.&(.~isnan.(pldata.fitori))), oriSelcv=sum((pldata.oricv.<oricvThres).&(pldata.oriauc.>oriaucThres).&pldata.visResp), oriSelosi=sum((pldata.osi.>osiThres).&(pldata.oriauc.>oriaucThres).&pldata.visResp),
            dirResp=sum((pldata.dirauc.>diraucThres).&pldata.visResp), dirRespfit=sum((pldata.dirauc.>diraucThres).&pldata.visResp.&(.~isnan.(pldata.fitdir))), dirSelcv=sum((pldata.dircv.<dircvThres).&(pldata.dirauc.>diraucThres).&pldata.visResp), dirSeldsi=sum((pldata.dsi.>dsiThres).&(pldata.dirauc.>diraucThres).&pldata.visResp),
            hueaxResp=sum((pldata.hueaxauc.>hueaucThres).&pldata.visResp), hueaxRespfit=sum((pldata.hueaxauc.>hueaucThres).&pldata.visResp.&(.~isnan.(pldata.fithueax))),hueaxSelcv=sum((pldata.hueaxcv.<huecvThres).&(pldata.hueaxauc.>hueaucThres).&pldata.visResp), hueaxSelsi=sum((pldata.hueaxsi.>huesiThres).&(pldata.hueaxauc.>hueaucThres).&pldata.visResp),
            huediResp=sum((pldata.huediauc.>hueaucThres).&pldata.visResp), huediRespfit=sum((pldata.huediauc.>hueaucThres).&pldata.visResp.&(.~isnan.(pldata.fithuedi))),huediSelcv=sum((pldata.huedicv.<huecvThres).&(pldata.huediauc.>hueaucThres).&pldata.visResp), huediSelsi=sum((pldata.huedisi.>huesiThres).&(pldata.huediauc.>hueaucThres).&pldata.visResp),
            hueoriResp=sum((pldata.hueoriauc.>oriaucThres).&pldata.visResp), hueoriRespfit=sum((pldata.hueoriauc.>oriaucThres).&pldata.visResp.&(.~isnan.(pldata.fithueori))),hueoriSelcv=sum((pldata.hueoricv.<oricvThres).&(pldata.hueoriauc.>oriaucThres).&pldata.visResp), hueoriSeldsi=sum((pldata.hueosi.>osiThres).&(pldata.hueoriauc.>oriaucThres).&pldata.visResp),
            huedirResp=sum((pldata.huedirauc.>diraucThres).&pldata.visResp), huedirRespfit=sum((pldata.huedirauc.>diraucThres).&pldata.visResp.&(.~isnan.(pldata.fithuedir))),huedirSelcv=sum((pldata.huedircv.<dircvThres).&(pldata.huedirauc.>diraucThres).&pldata.visResp), huedirSeldsi=sum((pldata.huedsi.>dsiThres).&(pldata.huedirauc.>diraucThres).&pldata.visResp),
            hueaucResp=sum(((pldata.hueaxauc.>hueaucThres).|(pldata.huediauc.>hueaucThres)).&pldata.visResp), hueSelcpi=sum(((pldata.hueaxauc.>hueaucThres).|(pldata.huediauc.>hueaucThres)).&pldata.visResp.&(pldata.huecpi.>cpiThres)))
        append!(planeSum, planesum)
        save(joinpath(dataExportFolder2,join([subject,"_",locId,"_",paraName,"_sum.jld2"])), locId,pldata,string(locId,"_sum"),planesum)
        CSV.write(joinpath(dataExportFolder2,join([subject,"_",locId,"_sum.csv"])), pldata)
    end
    CSV.write(joinpath(dataExportFolder2,join([subject,"_",recordSession[i],"_",paraName,"_sum.csv"])), planeSum)
end

exptNum = exptOriNum + exptHueNum
# CSV.write(joinpath(dataExportFolder1,join([subject,"_hue.csv"])), hueData)
