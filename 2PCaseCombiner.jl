using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact,CSV,MAT,Query,DataStructures, HypothesisTests, StatsFuns, Random
# Expt info
dataExportFolder = "O:\\AF4\\COFD_manuscript"
animalId = ["AE6","AE7","AE8","AE9","AF3","AF4","AF5"]
Manu_title = "Functional Organization for Color Appearance Mechanisms in Primary Visual Cortex"
oriaucThres = 0.7
hueaucThres = 0.8

resultpath=matchfile(Regex("_Summary.jld2"), dir=dataExportFolder,join=true)
if length(resultpath) > 0
    result=load(resultpath[1], "result")
else
    result=Dict()
end
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

## Summary of cell Num in all cases
oriData=DataFrame(); hueData=DataFrame(); staData=DataFrame();
for i=1:size(oripath,1)
    append!(oriData, CSV.read(oripath[i]))
end
for i=1:size(huepath,1)
    append!(hueData, CSV.read(huepath[i]))
end
for i=1:size(stapath,1)
    append!(staData, CSV.read(stapath[i]))
end


visResp = oriData.visResp .| hueData.visResp
oriSlec = (oriData.oriauc .> oriaucThres) .& visResp
hueSlec = ((hueData.hueaxauc.>hueaucThres) .| (hueData.huediauc.>hueaucThres)) .& visResp
staSig = staData.iscone .| staData.isachro

result["title"] = Manu_title
result["animalId"] = animalId
result["oriaucThres"] = oriaucThres
result["hueaucThres"] = hueaucThres
result["rawData"]=Dict()
result["rawData"]["oripath"] = oripath
result["rawData"]["huepath"] = huepath
result["rawData"]["stapath"] = stapath
result["rawData"]["oriData"] = oriData
result["rawData"]["hueData"] = hueData
result["rawData"]["staData"] = staData
result["cellNumSum"]=Dict()
result["cellNumSum"]["visRespNum"] = sum(visResp)
result["cellNumSum"]["oriNum"] = sum(oriSlec)
result["cellNumSum"]["hueNum"] = sum(hueSlec)
result["cellNumSum"]["AllstaNum"] = sum(staSig)
result["cellNumSum"]["LstaNum"] = sum(staData.isl)
result["cellNumSum"]["MstaNum"] = sum(staData.ism)
result["cellNumSum"]["SstaNum"] = sum(staData.iss)
result["cellNumSum"]["AstaNum"] = sum(staData.isachro)
result["cellNumSum"]["staProportion"] = sum(staSig) / sum(visResp)

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

result["rawData"]["DKLhueData"] = DKLhueData
result["rawData"]["OridklData"] = OridklData
result["DKLcellNumSum"]=Dict()
result["DKLcellNumSum"]["DKLcellNum"] = sum(DKLhueInd)
result["DKLcellNumSum"]["SaxisNum"] = sum(SaxisDKL)
result["DKLcellNumSum"]["SaxiscellPro"] = sum(SaxisDKL) / sum(DKLhueInd)

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
result["rawData"]["AF5V1OnOffResltion"] = AF5
result["AF5LMpixelprop"] = sum(pLM[26:50, 1:25]) + sum(pLM[1:25,26:50])

## AF4 2P cells within ISI COFD distribution
AF4matpath = matchfile(Regex("AF4_U002U003U004U005U006_hue_ISI Stength_struct.mat"), dir=dataExportFolder,join=true)
cell_COFD = matread(AF4matpath[1])["pltstruct"]
inCOFD = cell_COFD["COFD"]
hueselective = cell_COFD["hueselect"]
Lcoord = cell_COFD["Lscale"]
Mcoord = cell_COFD["Mscale"]
Scoord = cell_COFD["Sscale"]
maxHue_cellCOFD = cell_COFD["maxhue"]
hueCellInCOFDInd = (inCOFD .== 1) .& (hueselective .== 1)
LonMoffQuaInd = ((Lcoord .> 0) .& (Mcoord .< 0)) .& hueCellInCOFDInd
MonLoffQuaInd = ((Lcoord .< 0) .& (Mcoord .> 0)) .& hueCellInCOFDInd
LMoppoQuaInd = LonMoffQuaInd .| MonLoffQuaInd

LonMoffCell_hue = maxHue_cellCOFD[LonMoffQuaInd]
MonLoffCell_hue = maxHue_cellCOFD[MonLoffQuaInd]

LonMoffCell_Lhue = sum(((0 .<= LonMoffCell_hue.<= 60) .| (300 .<= LonMoffCell_hue .<= 330)))
LonMoffCell_Mhue = sum((120 .<= LonMoffCell_hue .<= 240))

MonLoffCell_Lhue = sum(((0 .<= MonLoffCell_hue.<= 60) .| (300 .<= MonLoffCell_hue .<= 330)))
MonLoffCell_Mhue = sum((120 .<= MonLoffCell_hue .<= 240))
#
# sum(LMoppoQuaInd)
# sum(hueselective)
# sum(30 .<= LonMoffCell .<= 150)

result["rawData"]["AF4_2Pcell_ISICOFD_distribution"] = cell_COFD
result["AF4_2Pcell_ISICOFD"]=Dict()
result["AF4_2Pcell_ISICOFD"]["AF4_2Pcell_ISICOFD_LMquaCellProportion"] = sum(LMoppoQuaInd) / sum(hueCellInCOFDInd)
result["AF4_2Pcell_ISICOFD"]["LonMoffCell_hue"] = length(LonMoffCell_hue)
result["AF4_2Pcell_ISICOFD"]["MonLoffCell_hue"] = length(MonLoffCell_hue)
result["AF4_2Pcell_ISICOFD"]["LonMoffCell_Lhue"] = LonMoffCell_Lhue
result["AF4_2Pcell_ISICOFD"]["LonMoffCell_Mhue"] = LonMoffCell_Mhue
result["AF4_2Pcell_ISICOFD"]["MonLoffCell_Mhue"] = MonLoffCell_Mhue
result["AF4_2Pcell_ISICOFD"]["MonLoffCell_Lhue"] = MonLoffCell_Lhue

# result["LonMoffCell_Lhue"] / result["LonMoffCell_hue"]
LonMoffCell_Lhue/length(LonMoffCell_hue)
LonMoffCell_Mhue/length(LonMoffCell_hue)

MonLoffCell_Mhue/length(MonLoffCell_hue)
MonLoffCell_Lhue/length(MonLoffCell_hue)

FisherExactTest(LonMoffCell_Lhue, length(LonMoffCell_hue), LonMoffCell_Mhue, length(LonMoffCell_hue))
FisherExactTest(MonLoffCell_Mhue, length(MonLoffCell_hue), MonLoffCell_Lhue, length(MonLoffCell_hue))
## AF4 Diamond plot (Only cells have L, M, or S STA): Cell proportion
STADKLpath = matchfile(Regex("STADKL_List.csv"), dir=dataExportFolder,join=true)

STADKLData = DataFrame();
for i=1:size(STADKLpath,1)
    append!(STADKLData, CSV.read(STADKLpath[i]))
end

LonMoff_tlt_Ind = (STADKLData.lcw .> 0) .& (STADKLData.mcw .< 0)
LonMoff_hue_Ind = LonMoff_tlt_Ind .& (STADKLData.HueSelec .== true)
LonMoff_Lhue_Ind = LonMoff_hue_Ind .& ((0 .<= STADKLData.maxhue .<= 60) .| (300 .<= STADKLData.maxhue .<= 330))

result["rawData"]["STAcellDKLhuePrefer"] = STADKLData
result["AF4_STA_Diamond"]=Dict()
result["AF4_STA_Diamond"]["LhueCellNum_inLonMoff"] = sum(LonMoff_Lhue_Ind)
result["AF4_STA_Diamond"]["LonMoffQuaTotalCellNum"] = sum(LonMoff_tlt_Ind)
result["AF4_STA_Diamond"]["LonMoffQuaHueCellNum"] = sum(LonMoff_hue_Ind)

MonLoff_tlt_Ind = (STADKLData.lcw .< 0) .& (STADKLData.mcw .> 0)
MonLoff_hue_Ind = MonLoff_tlt_Ind .& (STADKLData.HueSelec .== true)
MonLoff_Mhue_Ind = MonLoff_hue_Ind .& (120 .<= STADKLData.maxhue .<= 240)

result["AF4_STA_Diamond"]["MhueCellNum_inMonLoff"] = sum(MonLoff_Mhue_Ind)
result["AF4_STA_Diamond"]["MonLoffQuaTotalCellNum"] = sum(MonLoff_tlt_Ind)
result["AF4_STA_Diamond"]["MonLoffQuaHueCellNum"] = sum(MonLoff_hue_Ind)

# result["LhueCellNum_inLonMoff"]/result["LonMoffQuaTotalCellNum"]
# result["MhueCellNum_inMonLoff"]/result["MonLoffQuaTotalCellNum"]

# result["LhueCellNum_inLonMoff"]/result["LonMoffQuaHueCellNum"]
# result["MhueCellNum_inMonLoff"]/result["MonLoffQuaHueCellNum"]

## S cone dominant cells' hue preference
DKLstapath = matchfile(Regex("AF4_[A-Za-z0-9]*_thres0.25_staData.csv"), dir=dataExportFolder,join=true)
DKLSpath = matchfile(Regex("Sonoff_List.csv"), dir=dataExportFolder,join=true)
DKLstaData = DataFrame(); DKLSData = DataFrame()
for i=1:size(DKLstapath,1)
    append!(DKLstaData, CSV.read(DKLstapath[i]))
end
for i=1:size(DKLSpath,1)
    append!(DKLSData, CSV.read(DKLSpath[i]))
end

result["rawData"]["DKLstaData"] = DKLstaData
result["rawData"]["DKLSData"] = DKLSData

DKLhueData = result["rawData"]["DKLhueData"]
# SonId = DKLstaData.cellId[DKLstaData.iss .& (DKLstaData.scwmn.>= 0)]
# SoffId = DKLstaData.cellId[DKLstaData.iss .& (DKLstaData.scwmn.< 0)]
sum(DKLstaData.iscone)
aa=[parse(Int,DKLhueData.locId[ii][3]) for ii in 1:size(DKLhueData,1)]
bb=[parse(Int,DKLhueData.locId[ii][7]) for ii in 1:size(DKLhueData,1)]
allDKLId = hcat(aa,bb,DKLhueData.cellId)
SstaId = hcat(DKLstaData.unitId[DKLstaData.iss], DKLstaData.planeId[DKLstaData.iss], DKLstaData.cellId[DKLstaData.iss])
# Find non-S cone cells from all cells
nonSstaInd = trues(size(DKLhueData,1))
for ii in 1:size(SstaId,1)
    for jj in 1:size(DKLhueData,1)
        if allDKLId[jj,:] == SstaId[ii,:]
            nonSstaInd[jj] = false
        end
    end
end
nonSstaInd = nonSstaInd .& DKLvisResp
SstaInd = (.!nonSstaInd) .& DKLvisResp
sum(nonSstaInd)
sum(SstaInd)
# DKL hue selective cells
DKLhueInd = ((DKLhueData.hueaxauc.>hueaucThres) .| (DKLhueData.huediauc.>hueaucThres)) .& DKLvisResp
sum(DKLhueInd)

# non-S cone STAs and DKL hue selective cells
nonSDKLInd = DKLhueInd .& nonSstaInd
# sum(nonSDKLInd)
SDKLInd = DKLhueInd .& (.!nonSstaInd)
# sum(SDKLInd)
nonSDKLmaxHue = DKLhueData.maxhue[nonSDKLInd]
SDKLmaxHue = DKLhueData.maxhue[SDKLInd]
SonDKLcell = DKLSData.maxhue[DKLSData.sign .== 1]
SoffDKLcell = DKLSData.maxhue[DKLSData.sign .== -1]
result["AF4_Scone_DKL"]=Dict()
result["AF4_Scone_DKL"]["SonDKLNum"] = length(SonDKLcell)
result["AF4_Scone_DKL"]["SoffDKLNum"] = length(SoffDKLcell)
result["AF4_Scone_DKL"]["nonSDKLmaxHue"] = nonSDKLmaxHue
result["AF4_Scone_DKL"]["SDKLmaxHue"] = SDKLmaxHue
# length(nonSDKLmaxHue)
# a=60
# b=120
# c=240
# d=300

## 90 in Son vs. 90 in nonS cells

# SdirSon = sum(a .<= SonDKLcell .<= b)
# SdirNonS = sum(a .<= nonSDKLmaxHue .<= b)
SdirSon = sum(SonDKLcell .== 90)
SdirNonS = sum(nonSDKLmaxHue .== 90)
FisherExactTest(SdirSon, length(SonDKLcell), SdirNonS, length(nonSDKLmaxHue))
# SdirSon/length(SonDKLcell)
# SdirNonS/length(nonSDKLmaxHue)
result["AF4_Scone_DKL"]["SdirSon"] = SdirSon
result["AF4_Scone_DKL"]["SdirNonS"] = SdirNonS

result["AF4_Scone_DKL"]["nonSDKLmaxHue"] = length(nonSDKLmaxHue)
## 270 in Soff vs. 270 in nonS cells

# SoffdirSoff = sum(c .<= SoffDKLcell .<= d)
# SoffdirNonS = sum(c .<= nonSDKLmaxHue .<= d)
SoffdirSoff = sum(SoffDKLcell .== 270)
SoffdirNonS = sum(nonSDKLmaxHue .== 270)
FisherExactTest(SoffdirSoff, length(SoffDKLcell), SoffdirNonS, length(nonSDKLmaxHue))
# SoffdirSoff/length(SoffDKLcell)
# SoffdirNonS/length(nonSDKLmaxHue)
#

# nonSconeId = DKLstaData.cellId[DKLstaData.iss .== false
#

# sum(30 .<= SonDKLcell .<= 150)/length(SonDKLcell)
# sum(30 .<= SoffDKLcell .<= 150)/length(SoffDKLcell)
#
# sum(210 .<= SonDKLcell .<= 330)/length(SonDKLcell)
# sum(210 .<= SoffDKLcell .<= 330)/length(SoffDKLcell)
# a=60
# b=120
# c=240
# d=300
result["AF4_Scone_DKL"]["SoffdirSoff"] = SoffdirSoff
result["AF4_Scone_DKL"]["SoffdirNonS"] = SoffdirNonS

## 90 in Son vs. 90 in Soff cells

# SonSon = sum(a .<= SonDKLcell .<= b)
# SonSoff = sum(a .<= SoffDKLcell .<= b)
SdirSon = sum(SonDKLcell .== 90)
SdirSoff = sum(SoffDKLcell .== 90)
# FisherExactTest(SonSon, length(SonDKLcell)-SonSon, SonSoff, length(SoffDKLcell)-SonSoff)
FisherExactTest(SdirSon, length(SonDKLcell), SdirSoff, length(SoffDKLcell))

# SoffSon = sum(c .<= SonDKLcell .<= d)
# SoffSoff = sum(c .<= SoffDKLcell .<= d)
result["AF4_Scone_DKL"]["SdirSoff"] = SdirSoff
## 270 in Son vs. 270 in Soff cells

SoffdirSon = sum(SonDKLcell .== 270)
SoffdirSoff = sum(SoffDKLcell .== 270)
FisherExactTest(SoffdirSon, length(SonDKLcell), SoffdirSoff, length(SoffDKLcell))

# FisherExactTest(SonSon,SoffSon,SonSoff,SoffSoff)
# pyplot()
# histogram(SonDKLcell,bins=range(-15,stop=345,length=13), normalized=:probability)
# histogram(SoffDKLcell,bins=range(-15,stop=345,length=13), normalized=:probability)
#
# DKLhueDis = filter(!iszero,proportions(Int64.(DKLhuecell)))
# SonDis = filter(!iszero,proportions(Int64.(SonDKLcell)))
# SoffDis = filter(!iszero,proportions(Int64.(SoffDKLcell)))
#
# SonDisNorm = SonDis ./ DKLhueDis
# SoffDisNorm = SoffDis ./ DKLhueDis
#
# FisherExactTest(Int64.(round.((sum(SonDisNorm[2:6]),sum(SonDisNorm[8:12]),sum(SoffDisNorm[2:6]),sum(SoffDisNorm[8:12])))))
# FisherExactTest(9,4,4,9)
result["AF4_Scone_DKL"]["SoffdirSon"] = SoffdirSon
result["AF4_Scone_DKL"]["SoffdirSoff"] = SoffdirSoff

## CO intensity in COFD vs nonCOFD

COpath = matchfile(Regex("V1_ConeNoncone_result.mat"), dir=dataExportFolder,join=true)
COpath2 = matchfile(Regex("V1_COintensity_stats.mat"), dir=dataExportFolder,join=true)
COintensity=Dict(); COFDmedDistr=[]; COFDmedInt=[]; nonCOFDmedInt=[];
allCOFDstats=DataFrame(); ani=Any[]; tests=Any[];pv=[];

for i = 1:length(COpath)
    # i=1
    case_result = matread(COpath[i])["result"]
    case_result2 = matread(COpath2[i])["result"]
    cofd = filter!(!isnan, dropdims(case_result2["LMS"],dims=2))
    noncofd = median(dropdims(case_result["noncone"],dims=1)) .* ones(size(cofd))
    push!(ani,split(COpath[i],"\\")[4])
    push!(tests,SignedRankTest(cofd,noncofd))
    # SignTest(cofd,noncofd)
    # push!(tests,MannWhitneyUTest(cofd,noncofd))
    push!(pv, pvalue(SignedRankTest(cofd,noncofd)))
    push!(COFDmedDistr,cofd)
    push!(COFDmedInt,median(case_result["cone"],dims=2))
    push!(nonCOFDmedInt, median(case_result["noncone"],dims=2))
end

allCOFDstats.animal=ani
allCOFDstats.COFDnonCOFDtest=tests
allCOFDstats.pvalue=pv
allCOFDstats.COFDmedDistr = COFDmedDistr
allCOFDstats.COFDIntensityMedian = COFDmedInt
allCOFDstats.nonCOFDIntensityMedian = nonCOFDmedInt
COintensity["individualCaseAllCOFD"]=allCOFDstats
cofd=dropdims(reduce(hcat,COFDmedInt),dims=1)
noncofd=dropdims(reduce(hcat,nonCOFDmedInt),dims=1)
COintensity["allcase"]=SignTest(cofd,noncofd)

# show(allCOFDstats, allrows = true, allcols = true)
## CO intensity in COFD subregions

COpath2 = matchfile(Regex("V1_COintensity_stats.mat"), dir=dataExportFolder,join=true)
Lon=[]; Loff=[]; Mon=[]; Moff=[]; Son=[]; Soff=[]; Aon=[]; Aoff=[];
COFDsubRegstats=DataFrame(); ani=[]; testname=Any[]; tests=Any[]; pv=[];

for i = 2:length(COpath2)
    # i=2
    case_result = matread(COpath2[i])["result"]
    lon = filter!(!isnan,dropdims(case_result["Lon"],dims=2))
    loff = filter!(!isnan,dropdims(case_result["Loff"],dims=2))
    mon = filter!(!isnan,dropdims(case_result["Mon"],dims=2))
    moff = filter!(!isnan,dropdims(case_result["Moff"],dims=2))
    son = filter!(!isnan,dropdims(case_result["Son"],dims=2))
    soff = filter!(!isnan,dropdims(case_result["Soff"],dims=2))
    # aon = filter!(!isnan,dropdims(case_result["WB_on"],dims=2))
    # aoff = filter!(!isnan,dropdims(case_result["WB_off"],dims=2))
    L = filter!(!isnan,dropdims(case_result["Liso"],dims=2))
    M = filter!(!isnan,dropdims(case_result["Miso"],dims=2))
    S = filter!(!isnan,dropdims(case_result["Siso"],dims=2))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"LM")
    push!(pv,pvalue(MannWhitneyUTest(L,M)))
    push!(tests,MannWhitneyUTest(L,M))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"SL")
    push!(pv,pvalue(MannWhitneyUTest(S,L)))
    push!(tests,MannWhitneyUTest(S,L))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"SM")
    push!(pv,pvalue(MannWhitneyUTest(S,M)))
    push!(tests,MannWhitneyUTest(S,M))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"LonLoff")
    push!(pv,pvalue(MannWhitneyUTest(lon,loff)))
    push!(tests,MannWhitneyUTest(lon,loff))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"MonMoff")
    push!(pv,pvalue(MannWhitneyUTest(mon,moff)))
    push!(tests,MannWhitneyUTest(mon,moff))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"SonSoff")
    push!(pv,pvalue(MannWhitneyUTest(son,soff)))
    push!(tests,MannWhitneyUTest(son,soff))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"SonLon")
    push!(pv,pvalue(MannWhitneyUTest(son,lon)))
    push!(tests,MannWhitneyUTest(son,lon))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"SonMon")
    push!(pv,pvalue(MannWhitneyUTest(son,mon)))
    push!(tests,MannWhitneyUTest(son,mon))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"SoffLon")
    push!(pv,pvalue(MannWhitneyUTest(soff,lon)))
    push!(tests,MannWhitneyUTest(soff,lon))
    push!(ani,split(COpath[i],"\\")[4])
    push!(testname,"SoffMon")
    push!(pv,pvalue(MannWhitneyUTest(soff,mon)))
    push!(tests,MannWhitneyUTest(soff,mon))
end
COFDsubRegstats.ani = ani
COFDsubRegstats.COFDsubRegtestName = testname
COFDsubRegstats.COFDsubRegtestPv = pv
COFDsubRegstats.COFDsubRegtest = tests
CSV.write(joinpath(dataExportFolder,join([Manu_title,"COFDsubReg_stats_result.csv"])), COFDsubRegstats)

COintensity["individualCaseCOFDsubReg"]=COFDsubRegstats
result["COintensity"] = COintensity


save(joinpath(dataExportFolder,join([Manu_title,"_stats_Summary.jld2"])),"result",result)
