using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact,CSV,MAT,Query,DataStructures, HypothesisTests, StatsFuns, Random
# Expt info
dataExportFolder = "O:\\AF4\\COFD_manuscript"
animalId = ["AE6","AE7","AF4"]
Manu_title = "Functional Organization for Color Appearance Mechanisms in Primary Visual Cortex"
oriaucThres = 0.7
hueaucThres = 0.8

oripath=matchfile(Regex("oriData.csv"), dir=dataExportFolder,join=true)
huepath=matchfile(Regex("hueData.csv"), dir=dataExportFolder,join=true)
stapath=matchfile(Regex("staData.csv"), dir=dataExportFolder,join=true)

## covert String to Bitarray
# aa=CSV.read(huepath[1])
# aa.visResp =  aa.visResp .== "TRUE"
# aa.responsive =  aa.responsive .== "TRUE"
# aa.modulative =  aa.modulative .== "TRUE"
# aa.isl = aa.isl .== "TRUE"
# aa.ism = aa.ism .== "TRUE"
# aa.iss = aa.iss .== "TRUE"
# aa.isa = aa.isa .== "TRUE"
# aa.iscone = aa.iscone .== "TRUE"
# aa.isachro = aa.isachro .== "TRUE"
# CSV.write(joinpath(dataExportFolder,join(["AF6_hueData.csv"])), aa)

## Load data
oriData=DataFrame(); hueData=DataFrame(); staData=DataFrame();result=Dict();
for i=1:size(oripath,1)
    append!(oriData, CSV.read(oripath[i]))
end
for i=1:size(huepath,1)
    append!(hueData, CSV.read(huepath[i]))
end
for i=1:size(stapath,1)
    append!(staData, CSV.read(stapath[i]))
end

## Summary of cell Num in all cases
visResp = oriData.visResp .| hueData.visResp
oriSlec = (oriData.oriauc .> oriaucThres) .& visResp
hueSlec = ((hueData.hueaxauc.>hueaucThres) .| (hueData.huediauc.>hueaucThres)) .& visResp
staSig = staData.iscone .| staData.isachro

result["title"] = Manu_title
result["animalId"] = animalId
result["oriaucThres"] = oriaucThres
result["hueaucThres"] = hueaucThres
result["oripath"] = oripath
result["huepath"] = huepath
result["stapath"] = stapath

result["oriData"] = oriData
result["hueData"] = hueData
result["staData"] = staData

result["visRespNum"] = sum(visResp)
result["oriNum"] = sum(oriSlec)
result["hueNum"] = sum(hueSlec)
result["AllstaNum"] = sum(staSig)
result["LstaNum"] = sum(staData.isl)
result["MstaNum"] = sum(staData.ism)
result["SstaNum"] = sum(staData.iss)
result["AstaNum"] = sum(staData.isachro)
result["staProportion"] = sum(staSig) / sum(visResp)

## Cell Summary in DKL cases
DKLhuepath = matchfile(Regex("AF4_[A-Za-z0-9]*_hueData.csv"), dir=dataExportFolder,join=true)
Oridklpath = matchfile(Regex("AF4_[A-Za-z0-9]*_oriData.csv"), dir=dataExportFolder,join=true)



DKLhueData = DataFrame();OridklData = DataFrame();DKLstaData = DataFrame();DKLSData = DataFrame();
for i=1:size(DKLhuepath,1)
    append!(DKLhueData, CSV.read(DKLhuepath[i]))
end
for i=1:size(Oridklpath,1)
    append!(OridklData, CSV.read(Oridklpath[i]))
end


DKLvisResp =  DKLhueData.visResp .| OridklData.visResp
DKLhueInd = ((DKLhueData.hueaxauc.>hueaucThres) .| (DKLhueData.huediauc.>hueaucThres)) .& DKLvisResp
DKLhuecell = DKLhueData.maxhue[DKLhueInd]
SaxisDKL = ((DKLhueData.hueaxauc.>hueaucThres) .| (DKLhueData.huediauc.>hueaucThres)) .& DKLvisResp .& ((DKLhueData.maxhue.== 90) .| (DKLhueData.maxhue.==270))

result["DKLhueData"] = DKLhueData
result["OridklData"] = OridklData
result["DKLcellNum"] = sum(DKLhueInd)
result["SaxisNum"] = sum(SaxisDKL)
result["SaxiscellPro"] = sum(SaxisDKL) / sum(DKLhueInd)

## AF5 Cone correlation bin proportion
AF5matpath = matchfile(Regex("AF5_Left_V1_OnOff Relation_result.mat"), dir=dataExportFolder,join=true)
AF5 = matread(AF5matpath[1])
# totalPix = length(AF5["result"]["L"])
pLM = AF5["result"]["pLM"]
# pLS = AF5["result"]["pLS"]
# pMS = AF5["result"]["pMS"]

# plot()
# heatmap!(pLM)
# proportion of pixels in LM opponent quadrants.
result["AF5V1OnOffResltion"] = AF5
result["AF5LMpixelprop"] = sum(pLM[26:50, 1:25]) + sum(pLM[1:25,26:50])

## Diamond plot: Cell proportion
STADKLpath = matchfile(Regex("STADKL_List.csv"), dir=dataExportFolder,join=true)

STADKLData = DataFrame();
for i=1:size(STADKLpath,1)
    append!(STADKLData, CSV.read(STADKLpath[i]))
end

LonMoff_tlt_Ind = (STADKLData.lcw .> 0) .& (STADKLData.mcw .< 0)
LonMoff_hue_Ind = (STADKLData.lcw .> 0) .& (STADKLData.mcw .< 0) .& (STADKLData.HueSelec .== "True")
LonMoff_Lhue_Ind = (STADKLData.lcw .> 0) .& (STADKLData.mcw .< 0) .& (STADKLData.HueSelec .== "True") .& ((0 .<= STADKLData.maxhue .<= 60) .| (300 .<= STADKLData.maxhue .<= 330))

result["STAcellDKLhuePrefer"] = STADKLData
result["LhueCellNum"] = sum(LonMoff_Lhue_Ind)
result["LonMoffQuaTotalCellum"] = sum(LonMoff_tlt_Ind)
sum(LonMoff_Lhue_Ind)
result["LonMoffQuaHueCellum"] = sum(LonMoff_hue_Ind)


MonLoff_tlt_Ind = (STADKLData.lcw .< 0) .& (STADKLData.mcw .> 0)
MonLoff_hue_Ind = (STADKLData.lcw .< 0) .& (STADKLData.mcw .> 0) .& (STADKLData.HueSelec .== "True")
MonLoff_Lhue_Ind = (STADKLData.lcw .< 0) .& (STADKLData.mcw .> 0) .& (STADKLData.HueSelec .== "True") .& ((0 .<= STADKLData.maxhue .<= 60) .| (300 .<= STADKLData.maxhue .<= 330))

result["MhueCellNum"] = sum(MonLoff_Lhue_Ind)
result["MonLoffQuaTotalCellum"] = sum(MonLoff_tlt_Ind)

result["MonLoffQuaHueCellum"] = sum(MonLoff_hue_Ind)

sum(LonMoff_Lhue_Ind)/sum(LonMoff_tlt_Ind)
sum(MonLoff_Lhue_Ind)/sum(MonLoff_tlt_Ind)

sum(LonMoff_Lhue_Ind)/sum(LonMoff_hue_Ind)
sum(MonLoff_Lhue_Ind)/sum(MonLoff_hue_Ind)

## S cone dominant cells' hue preference
DKLstapath = matchfile(Regex("AF4_[A-Za-z0-9]*_thres0.25_staData.csv"), dir=dataExportFolder,join=true)
DKLSpath = matchfile(Regex("Sonoff_List.csv"), dir=dataExportFolder,join=true)
for i=1:size(DKLstapath,1)
    append!(DKLstaData, CSV.read(DKLstapath[i]))
end
for i=1:size(DKLSpath,1)
    append!(DKLSData, CSV.read(DKLSpath[i]))
end

SonId = DKLstaData.cellId[DKLstaData.iss .& (DKLstaData.scwmn.>= 0)]
SoffId = DKLstaData.cellId[DKLstaData.iss .& (DKLstaData.scwmn.< 0)]
SonDKLcell = DKLSData.maxhue[DKLSData.sign .== 1]
SoffDKLcell = DKLSData.maxhue[DKLSData.sign .== -1]
sum(30 .<= SonDKLcell .<= 150)/length(SonDKLcell)
sum(30 .<= SoffDKLcell .<= 150)/length(SoffDKLcell)

sum(210 .<= SonDKLcell .<= 330)/length(SonDKLcell)
sum(210 .<= SoffDKLcell .<= 330)/length(SoffDKLcell)
a=60
b=120
c=240
d=300
SonSon = sum(a .<= SonDKLcell .<= b)
SonSoff = sum(a .<= SoffDKLcell .<= b)
FisherExactTest(SonSon, length(SonDKLcell)-SonSon, SonSoff, length(SoffDKLcell)-SonSoff)

SoffSon = sum(c .<= SonDKLcell .<= d)
SoffSoff = sum(c .<= SoffDKLcell .<= d)
FisherExactTest(SonSon,SoffSon,SonSoff,SoffSoff)
pyplot()
histogram(SonDKLcell,bins=range(-15,stop=345,length=13), normalized=:probability)
histogram(SoffDKLcell,bins=range(-15,stop=345,length=13), normalized=:probability)

DKLhueDis = filter(!iszero,proportions(Int64.(DKLhuecell)))
SonDis = filter(!iszero,proportions(Int64.(SonDKLcell)))
SoffDis = filter(!iszero,proportions(Int64.(SoffDKLcell)))

SonDisNorm = SonDis ./ DKLhueDis
SoffDisNorm = SoffDis ./ DKLhueDis

FisherExactTest(Int64.(round.((sum(SonDisNorm[2:6]),sum(SonDisNorm[8:12]),sum(SoffDisNorm[2:6]),sum(SoffDisNorm[8:12])))))
FisherExactTest(9,4,4,9)

result["DKLstaData"] = DKLstaData
result["DKLSData"] = DKLSData
result["SonDKLNum"] = length(SonDKLcell)
result["SoffDKLNum"] = length(SoffDKLcell)



save(joinpath(dataExportFolder,join([Manu_title,"_cellNum_Summary.jld2"])),"result",result)
