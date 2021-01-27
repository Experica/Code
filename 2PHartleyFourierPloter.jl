using NeuroAnalysis,Statistics,StatsBase,FileIO,Images,Plots,LsqFit,FFTW

disk = "O:"
subject = "AF4"

#-----------------------
# recordSession = "002"  #AF2
# testId = ["004", "005","006","007"]

# recordSession = "003"  #AF2
# testId = ["004", "005","006","003"]

recordSession = "004"  #AF2
testId = [ "005","006","007","004"]

# recordSession = "005"  #AF2
# testId = ["003","004", "005","002"]

# recordSession = "006"  #AF2
# testId = ["005","007", "008","004"]

#-------------
# recordSession = "003"  #AF3
# testId = ["004","005","007", "003"]

#----------------
# recordSession = "002"  #AF4
# testId = ["008","009","010","011"]  # In the order of L, M, S, and achromatic

# recordSession = "003"  #AF4
# testId = ["000","001","002", "004"]

# recordSession = "004"  #AF4
# testId = ["001","002","003","004"]

# recordSession = "005"  #AF4
# testId = ["000","001","003","004"]

# recordSession = "006"  #AF4
# testId = ["000","001","002","003"]

# --------- AE7
# recordSession = "003"
# testId = ["002","003","004","001"]

# recordSession = "005"
# testId = ["004","005","006","003"]

# recordSession = "006"
# testId = ["004","006","007","003"]

# recordSession = "008"
# testId = ["003","004","005","002"]

# recordSession = "009"
# testId = ["002","003","004","001"]

# recordSession = "011"
# testId = ["002","004","005","001"]

# recordSession = "012"
# testId = ["002","003","004","001"]

# recordSession = "013"
# testId = ["003","004","005","002"]

#----------------------
# recordSession = "002"  #AE6
# testId = ["006", "007", "008","005"]    # 001
# testId = ["012", "013", "014","011"]    # 000

# recordSession = "003"  #AE6
# testId = ["005", "006","007","002"]   # 001
# testId = ["011", "012","013","010"]  # 000

# recordSession = "004"  #AE6
# testId = ["003", "007", "008","002"]   # 001 # In the order of L, M, S, and achromatic
# testId = ["013", "014", "015","012"]  # 000

recordPlane = "001"
delays = collect(-0.066:0.033:0.4)
print(collect(delays))

lbTime = 0.198
ubTime = 0.330
blkTime = 0.099
respThres = 0.25

# cw = datasetFinal["coneweight"]
# achroResp = datasetFinal["achroResp"]
# CSV.write(joinpath(resultFolder,join([subject,"_",recordSession,"_",recordPlane,"_thres",respThres,"_sta_dataset.csv"])), cw)
# CSV.write(joinpath(resultFolder,join([subject,"_",recordSession,"_",recordPlane,"_thres",respThres,"_achrosta_dataset.csv"])), achroResp)
## Prepare data & result path
siteId=[]
for i =1:size(testId,1)
    siteid = join(["$(recordSession)_", testId[i], "_$(recordPlane)"])
    push!(siteId,siteid)
end

dataFolder = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]))
dataExportFolder = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]), siteId, "DataExport")
resultFolder = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]), "_Summary", "DataExport")
resultFolderPlot = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]), "_Summary", join(["plane_",recordPlane]),"0. Original maps","Fourier")
isdir(resultFolder) || mkpath(resultFolder)
isdir(resultFolderPlot) || mkpath(resultFolderPlot)

## Load all stas
# testids = ["$(siteId)_HartleySubspace_$i" for i in 1:4]
dataFile=[]
for i=1:size(testId,1)
    datafile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*_tuning_result.jld2"), dir=dataExportFolder[i],join=true)[1]
    push!(dataFile, datafile)
end

## Join Results from L, M, S, and achromatic expts.
dataset = sbxjoinhartleyFourier(load.(dataFile))
save(joinpath(resultFolder,join([subject,"_",recordSession,"_",recordPlane,"_thres",respThres,"_fourier_dataset.jld2"])),"dataset",dataset)
## To load saved data
dataset = load(joinpath(resultFolder,join([subject,"_",recordSession,"_",recordPlane,"_thres",respThres,"_fourier_dataset.jld2"])),"dataset")
## Plot Hartley Fourier image

for u in sort(collect(keys(dataset["kern"])))
    # u=10
    colors =["L", "M", "S", "A"]
    kern = dataset["kern"][u]
    replace!(kern, -Inf=>0)
    maxmag = maximum(abs.(kern))
    kern = kern ./ maxmag   # Kernal is normalized by maximum of all kernals (L, M, S, A) for this cell
    # kern= normalized(kern;equal=false)
    delta = dataset["kdelta"]
    imagesize = size(kern)[1]
    truestimsz = 12   # Arbitrary
    colorbar_hack = zeros(size(kern))
    xylim = [0,round(truestimsz,digits=1)]
    xy = range(xylim...,length=imagesize)

    p = Plots.plot(layout=(1,5),legend=false,size=(1650,600))
    # Plots.bar!(p,subplot=1,dataset["color"],ucex,frame=:zerolines)
    foreach(c->Plots.heatmap!(p,subplot=c,xy,xy,kern[:,:,c],aspect_ratio=:equal,frame=:grid,
    color=:bwr,clims=(-1,1),xlims=xylim,ylims=xylim,xticks=[],yticks=[],yflip=true,xlabel=string(colors[c]),title="Cell_$(u)_Fourier"),1:4)
    heatmap!(p,subplot=5,xy,xy,colorbar_hack[:,:,1],aspect_ratio=:equal,frame=:grid,xticks=[],yticks=[],color=:bwr,clims=(-maxmag,maxmag),colorbar=:left) # add colorbar
    foreach(c->Plots.plot!(p,subplot=c,[6], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=3, linealpha=1, xticks=([6],["0"]), label=""),1:4)
    foreach(c->Plots.plot!(p,subplot=c,[6], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=3, linealpha=1, label=""),1:4)
    foreach(c->Plots.plot!(p,subplot=c,yticks=([6],["0"])),1)
    # :diverging_bwr_40_95_c42_n256   RdBu_5
    # foreach(c->Plots.plot!(p,subplot=c+1,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"]), label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"])),1)

    # foreach(c->Plots.plot!(p,subplot=c+1,[0.4,0.8,1.2,1.6,2.0], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.6,1.2,1.8,2.4],["-1.2","-0.6","0","0.6","1.2"]), label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,[0.4,0.8,1.2,1.6,2.0], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.6,1.2,1.8,2.4],["-1.2","-0.6","0","0.6","1.2"])),1)
    p
    savefig(joinpath(resultFolderPlot,join([subject,"_U",recordSession,"_Plane",recordPlane, "_Cell",u,".png"])))
end

## Plot example cell

##
