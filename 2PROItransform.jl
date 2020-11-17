using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,StatsPlots,Mmap,Images,StatsBase,Interact, CSV,MAT,Query,DataStructures, HypothesisTests, StatsFuns, Random

dataFolder = "J:\\AF3\\2P_analysis\\U005\\005_007_000\\DataExport"

subject = split(dataFolder, "\\")[2]
siteId = split(dataFolder, "\\")[5]

segmentFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.segment"),dir=dataFolder,join=true)[1]
alignFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_ot_[A-Za-z0-9]*.align"),dir=dataFolder,join=true)[1]

segment = prepare(segmentFile)
align = prepare(alignFile)
cellRoi = segment["vert"]
bkgImg = align["m"]

save(joinpath(dataFolder,join([subject,"_",siteId,"_roibkg.jld2"])), "roi",cellRoi, "bkg",bkgImg)
