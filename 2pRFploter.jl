using NeuroAnalysis,Statistics,StatsBase,FileIO,Images,Plots,LsqFit,FFTW

disk = "F:"
subject = "AF2"

# recordSession = "003"  #AF2
# testId = ["004", "005","006","003"]

recordSession = "004"  #AF2
testId = [ "005","006","007","004"]

# recordSession = "005"  #AF2
# testId = ["003","004", "005","002"]

# recordSession = "006"  #AF2
# testId = ["005","007", "008","004"]

# recordSession = "003"  #AF3
# testId = ["004","005","007", "003"]

# recordSession = "002"  #AF4
# testId = ["008","009","010","011"]  # In the order of L, M, S, and achromatic

# recordSession = "003"  #AF4
# testId = ["000","001","002", "004"]
# #
# recordSession = "004"  #AF4
# testId = ["001","002","003","004"]

# recordSession = "005"  #AF4
# testId = ["000","001","003","004"]

# recordSession = "006"  #AF4
# testId = ["000","001","002","003"]

# recordSession = "004"  #AE6
# testId = ["006", "007", "008","005"]
# testId = ["012", "013", "014","011"]
# testId = ["005", "006","007","002"]
# testId = ["011", "012","013","010"]
# testId = ["003", "007", "008","002"]   # In the order of L, M, S, and achromatic
# testId = ["013", "014", "015","012"]

recordPlane = "000"
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
# resultFolderPlot = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]), "_Summary", join(["plane_",recordPlane]),"0. Original maps","STA")
resultFolderPlot = joinpath.(disk,subject, "2P_analysis", join(["U",recordSession]), "_Summary", join(["plane_",recordPlane]),"0. Original maps","Fourier")
isdir(resultFolder) || mkpath(resultFolder)
isdir(resultFolderPlot) || mkpath(resultFolderPlot)

## Load all stas
# testids = ["$(siteId)_HartleySubspace_$i" for i in 1:4]
dataFile=[]
for i=1:size(testId,1)
    # datafile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*_sta.jld2"), dir=dataExportFolder[i],join=true)[1]
    datafile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*_tuning_result.jld2"), dir=dataExportFolder[i],join=true)[1]
    push!(dataFile, datafile)
end

## Join STAs from L, M, S, and achromatic expts. Also find the strongest response withe corresponding delay.
# dataset = sbxjoinsta(load.(dataFile),lbTime,ubTime,blkTime)
dataset = sbxjoinhartleyFourier(load.(dataFile))
save(joinpath(resultFolder,join([subject,"_",recordSession,"_",recordPlane,"_thres",respThres,"_fourier_dataset.jld2"])),"dataset",dataset)
## Check responsiveness of cell based on threshold and delay range, filter out low-response and cell 'response' too early or too late

dataset = sbxresponsivesta!(dataset,lbTime,ubTime,respThres)
dataset = sbxrffit!(dataset)
dataset = sbxgetbestconesta(dataset)
datasetFinal=sbxgetconeweight(dataset)
save(joinpath(resultFolder,join([subject,"_",recordSession,"_",recordPlane,"_thres",respThres,"_sta_datasetFinal.jld2"])),"datasetFinal",datasetFinal)

# stop
dataset = load(joinpath(resultFolder,join([subject,"_",recordSession,"_",recordPlane,"_thres",respThres,"_sta_datasetFinal.jld2"])),"datasetFinal")

## Plot Hartley Fourier image

for u in sort(collect(keys(dataset["kern"])))
    # u=53
    colors =["L", "M", "S", "A"]
    kern = dataset["kern"][u]
    replace!(kern, -Inf=>0)
    kern= normalized(kern;equal=false)
    delta = dataset["delta"]
    imagesize = size(kern)[1]
    # stisize = dataset["stisize"]
    # maskradius = dataset["maskradius"]
    # truestimsz = stisize*maskradius*2
    truestimsz = 12
    # ucex = map(i->i.ex,dataset["ulcex"][u])
    # ucd = map(i->i.pd,dataset["ulcex"][u])
    # clim = maximum(abs.(kern))    clims=(-clim,clim),
    xylim = [0,round(truestimsz,digits=1)]
    xy = range(xylim...,length=imagesize)

    p = Plots.plot(layout=(1,4),legend=false,size=(1650,600))
    # Plots.bar!(p,subplot=1,dataset["color"],ucex,frame=:zerolines)
    foreach(c->Plots.heatmap!(p,subplot=c,xy,xy,kern[:,:,c],aspect_ratio=:equal,frame=:grid,
    color=:bwr,clims=(-1,1),xlims=xylim,ylims=xylim,xticks=[],yticks=[],yflip=true,xlabel=string(colors[c]),title="Cell_$(u)_Fourier"),1:4)
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

## stas
# @manipulate for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
@manipulate for u in sort(collect(keys(dataset["ulsta"])))#,d in eachindex(dataset["delays"])
    usta=dataset["usta"][u]
    # delays = dataset["delays"]
    # d=8
    imagesize = dataset["imagesize"]
    stisize = dataset["stisize"]
    uresponsive = dataset["uresponsive"]
    clim = maximum(map(i->abs(i.ex),dataset["ucex"][u]))
    xlim = [0,5]#[0,stisize[2]]
    ylim = [0,5]#stisize[1]]
    x = range(xlim...,length=imagesize[2])
    y = range(ylim...,length=imagesize[1])

    p = Plots.plot(layout=(1,4),legend=false,size=(1600,600),title="Cell_$(u)_STA_$(delays[d])_Responsive_$(uresponsive[u])")
    foreach(c->Plots.heatmap!(p,subplot=c,x,y,usta[:,:,d,c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xlim,ylims=ylim,xticks=xlim,yticks=ylim,yflip=true,xlabel=dataset["color"][c]),1:4)
    p
    savefig("./temp/Unit_$u.png")
end
## best roi stas
 # for u in sort(collect(keys(dataset["ulsta"])))
for u in sort(collect(keys(dataset["ulsta"])))
    # u=11
    local delays
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    imagesize = size(ulsta)[1]
    stisize = dataset["stisize"]
    maskradius = dataset["maskradius"]
    truestimsz = stisize*maskradius*2
    ucex = map(i->i.ex,dataset["ulcex"][u])
    ucd = map(i->i.pd,dataset["ulcex"][u])
    clim = maximum(abs.(ucex))
    xylim = [0,round(truestimsz,digits=1)]
    xy = range(xylim...,length=imagesize)

    p = Plots.plot(layout=(1,5),legend=false,size=(1650,600))
    Plots.bar!(p,subplot=1,dataset["color"],ucex,frame=:zerolines)
    foreach(c->Plots.heatmap!(p,subplot=c+1,xy,xy,ulsta[:,:,ucd[c],c],aspect_ratio=:equal,frame=:grid,
    color=:bwr,clims=(-clim,clim),xlims=xylim,ylims=xylim,xticks=[],yticks=[],
    yflip=true,xlabel=string(dataset["color"][c]),title="Cell_$(u)_STA_$(delays[ucd[c]])"),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,[0.3,0.6,0.9,1.2,1.5], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.3,0.6,0.9,1.2,1.5,1.8],["-0.9","-0.6","-0.3","0","0.3","0.6","0.9"]), label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,[0.3,0.6,0.9,1.2,1.5], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.3,0.6,0.9,1.2,1.5,1.8],["-0.9","-0.6","-0.3","0","0.3","0.6","0.9"])),1)

    foreach(c->Plots.plot!(p,subplot=c+1,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"]), label=""),1:4)
    foreach(c->Plots.plot!(p,subplot=c+1,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
    foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"])),1)

    # foreach(c->Plots.plot!(p,subplot=c+1,[0.4,0.8,1.2,1.6,2.0], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, xticks=([0,0.6,1.2,1.8,2.4],["-1.2","-0.6","0","0.6","1.2"]), label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,[0.4,0.8,1.2,1.6,2.0], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label=""),1:4)
    # foreach(c->Plots.plot!(p,subplot=c+1,yticks=([0,0.6,1.2,1.8,2.4],["-1.2","-0.6","0","0.6","1.2"])),1)

    p

    savefig(joinpath(resultFolderPlot,join([subject,"_U",recordSession,"_Plane",recordPlane, "_Cell",u,".png"])))
end

sort(collect(keys(dataset["urf"])))
collect(keys(first(values(dataset["urf"]))))

## Plot example STA
u=399
coneIdex = 2
ulsta = dataset["ulsta"][u][:,:,10,coneIdex]

maskradius = dataset["maskradius"]
ppd=dataset["ppd"]
xi = dataset["xi"]
stisize = dataset["stisize"]

imagesize = size(ulsta)[1]
truestimsz = stisize*maskradius*2
ucex = map(i->i.ex,dataset["ulcex"][u])[coneIdex]
ucd = map(i->i.pd,dataset["ulcex"][u])[coneIdex]
clim = maximum(abs.(ucex))
xylim = [0,round(truestimsz,digits=1)]
xy = range(xylim...,length=imagesize)

p = Plots.plot(legend=false,size=(600,600))

Plots.heatmap!(p,xy,xy,ulsta,aspect_ratio=:equal,frame=:grid,
color=:bwr,clims=(-clim,clim),xlims=xylim,ylims=xylim,xticks=[],yticks=[],yflip=true)
Plots.plot!(p,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="vline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label="")
Plots.plot!(p,[0.2,0.4,0.6,0.8,1.0,1.2,1.4], seriestype="hline", linecolor=:gray, linestyle=:dot, linewidth=1, linealpha=0.5, label="")
# Plots.plot!(p,yticks=([0,0.4,0.8,1.2,1.6],["-0.8","-0.4","0","0.4","0.8"]))
p
save(joinpath(resultFolder,"AF4_U004_plane001_cell399_M.svg"),p)
a

# r = rand(-10:10,10,10)
# p=Plots.heatmap(r, color=:bwr)
# p
# save(joinpath(resultFolder,"bwr_colorbar.svg"),p)
## rf fit
# @manipulate
for u in sort(collect(keys(dataset["urf"])))#,m in collect(keys(first(values(dataset["urf"]))))
    nc = length(dataset["color"])
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    imagesize = size(ulsta)[1]
    ppd = dataset["ppd"]
    stisize = dataset["stisize"]
    maskradius = dataset["maskradius"]
    truestimsz = stisize*maskradius*2
    ucex = map(i->i.ex,dataset["ucex"][u])
    ucd = map(i->i.d,dataset["ucex"][u])
    clim = maximum(abs.(ucex))
    xylim = [0,round(truestimsz,digits=1)]
    xy = range(xylim...,length=imagesize)

    p = Plots.plot(layout=(2,nc),legend=false,size=(1650,900))
    foreach(c->Plots.heatmap!(p,subplot=c,xy,xy,ulsta[:,:,ucd[c],c],aspect_ratio=:equal,frame=:grid,color=:coolwarm,clims=(-clim,clim),
    xlims=xylim,ylims=xylim,xticks=[],yticks=[],yflip=true,title="Unit_$(u)_STA_$(delays[ucd[c]])"),1:nc)


    fit=dataset["urf"][u]
    pr = (imagesize[1]-1)/2
    x=y=(-pr:pr)./ppd
    sr = round(x[end],digits=1)
    xylim=[-sr,sr]
    # if m == :dog
    m = :dog
    rfs = map(f->isnothing(f) ? nothing : [rfdog(i,j,f.param...) for j in reverse(y), i in x],fit[m])
    foreach(c->isnothing(rfs[c]) ? Plots.plot!(p,subplot=c+nc,frame=:none) :
    Plots.heatmap!(p,subplot=c+nc,x,y,rfs[c],aspect_ratio=:equal,frame=:grid,color=:coolwarm,clims=(-clim,clim),
    xlims=xylim,ylims=xylim,xticks=[],yticks=[],yflip=true,xlabel=string(dataset["color"][c])),1:nc)
    # else
    #     rfs = map(f->isnothing(f) ? nothing : [rfgabor(i,j,f.param...) for j in reverse(y), i in x],fit[m])
    # end


    # foreach(c->isnothing(rfs[c]) ? Plots.plot!(p,subplot=c+2nc,frame=:none) :
    # Plots.histogram!(p,subplot=c+2nc,fit[m][c].resid,frame=:grid,color=:coolwarm,linecolor=:match,bar_width=1,xlims=[-abs(maximum(fit[m][c].resid)),abs(maximum(fit[m][c].resid))],
    # xlabel="Residual",title="",grid=false),1:nc)
    #
    # foreach(c->isnothing(rfs[c]) ? Plots.plot!(p,subplot=c+3nc,frame=:none) :
    # Plots.scatter!(p,subplot=c+3nc,vec(ulsta[:,:,ucd[c],c]),vec(rfs[c]),frame=:grid,color=:coolwarm,grid=false,
    # xlabel="y",ylabel="predict y",title="r = $(round(fit[m][c].r,digits=3))",markerstrokewidth=0,markersize=1),1:nc)
    p
    savefig(joinpath(resultFolderPlot,join([subject,"_U",recordSession,"_Plane",recordPlane, "_Cell",u,".png"])))
end









function getbestrf(dataset;rt=0.65)
    urf = dataset["urf"]
    nc = length(dataset["color"])
    ubrf=Dict()
    for u in keys(urf)
        # u!=96 && continue
        mfs=[]
        for m in keys(urf[u])
            push!(mfs, map(f->begin
             if isnothing(f)
                 nothing
             else
                 (param=f.param,r=f.r,m=m)
             end
         end,urf[u][m]))
         end
         # display(mfs)
         crf = []
         for i in 1:nc
             push!(crf,mapreduce(mf->mf[i],(m1,m2)->begin
             if isnothing(m1)
                 if isnothing(m2)
                     nothing
                 else
                     m2.r > rt ? m2 : nothing
                 end
             else
                 if isnothing(m2)
                     m1.r > rt ? m1 : nothing
                 else
                     t = m1.r > m2.r ? m1 : m2
                     t.r > rt ? t : nothing
                 end
             end
         end,mfs))
         end
         ubrf[u]=crf
    end
    return ubrf
end


usrf = getbestrf(dataset)


urf = map(i->i[:gabor][2],values(dataset["urf"]))
urft = map(i->i.converged,urf)
vurf = urf[urft]

vups = mapreduce(i->i.param,hcat,vurf)

using Clustering

r = kmeans(vups[[3,5],:],4,display=:final)

r = kmeans([abs.(vups[3:3,:]./vups[5:5,:].-1);maximum(vups[[3,5],:],dims=1).*vups[7:7,:]],3,display=:final)

r = kmeans(maximum(vups[[3,5],:],dims=1).*vups[7:7,:],2,display=:final)


ci = assignments(r)


cu = map(i->us[findall(i.==ci)],1:3)


sort(cu[1])

Plots.scatter(vups[3,:],vups[7,:],group=ci)



# spike = load(joinpath(resultFolder,testids[1],"spike.jld2"),"spike")
# unitid = spike["unitid"];unitgood = spike["unitgood"];unitposition = spike["unitposition"];unitlayer = assignlayer(unitposition[:,2],layer)

us = collect(keys(dataset["ulsta"]))
# up = unitposition[indexin(us,unitid),:]
# ul = unitlayer[indexin(us,unitid)]



uc = [argmax(map(j->abs(j.ex),dataset["ucex"][u])) for u in keys(dataset["ulsta"])]
ud = map((u,c)->dataset["ucex"][u][c].d,keys(dataset["ulsta"]),uc)

ucsta = map((u,d,c)->mapcolor(dataset["ulsta"][u][:,:,d,c],dataset["minmaxcg"][c]),keys(dataset["ulsta"]),ud,uc)
ucmsta = map(i->alphamask(i,radius=0.5,sigma=10,masktype="Disk")[1],ucsta)


p=plotunitpositionimage(up,ucmsta,layer=layer)
save("UnitPosition_STA.svg",p)

p=plotunitlayerimage(ul,ucmsta,unitid=us)
save("UnitLayer_STA.svg",p)


uexs = mapreduce(u->map(j->j.ex,dataset["ucex"][u]),hcat,us)
uds = mapreduce(u->map(j->j.d,dataset["ucex"][u]),hcat,us)

function rftype(cws,cpd)
    a=cws[1];l=cws[2];m=cws[3];s=cws[4];c=1;d = cpd[argmax(abs.(exs))]
    t = join(filter(k->!isnothing(k),map((i,j)->i==0 ? nothing : "$j$(i > 0 ? "+" : "-")",cws,["A","L","M","S"])),"_")
    if l != 0
        if m != 0 # both significent l,m
            if l > 0 && m < 0
                t = "L+M-"
                c = 2
                d = cpd[argmax(abs.(cws[2:3]))+1]
            elseif l < 0 && m > 0
                t = "M+L-"
                c = 3
                d = cpd[argmax(abs.(cws[2:3]))+1]
            elseif l > 0 && m > 0
                if s < 0 && abs(s) > maximum(abs.(cws[2:3]))
                    t = "S-L+M+"
                    c = 4
                    d = cpd[4]
                else
                    t = "+-"
                    c = 1
                    d = cpd[1]
                end
            elseif l < 0 && m < 0
                if s > 0 && abs(s) > maximum(abs.(cws[2:3]))
                    t = "S+L-M-"
                    c = 4
                    d = cpd[4]
                else
                    t = "+-"
                    c = 1
                    d = cpd[1]
                end
            end
        else
            if a ==0
                t = "L+"
                c = 4
                d = cpd[4]
            else
            end
        end
    end

        if l>0 && m<0
            t = "L+M-"
            c = 2
            d = ds[argmax(abs.(exs[2:3]))+1]
        elseif l<0 && m>0
            t="M+L-"
            c=3
            d = ds[argmax(abs.(exs[2:3]))+1]
        elseif l>0 && m>0
            if s<0
                t="S-L+M+"
                c=4
                d = ds[argmax(abs.(exs[2:end]))+1]
            else
                t = "+-"
                c=1
                d = ds[argmax(abs.(exs))]
            end
        elseif l<0 && m<0
            if s>0
                t="S+L-M-"
                c=4
                d = ds[argmax(abs.(exs[2:end]))+1]
            else
                t = "-+"
                c=1
                d = ds[argmax(abs.(exs))]
            end
        end
    return (t=t,c=c,d=d)
end

tcd = [rftype(uexs[:,i],uds[:,i]) for i in 1:size(uexs,2)]
uc = map(i->i.c,tcd)
ud = map(i->i.d,tcd)

_,nsg,_,_=histrv(unitposition[unitgood,2],0,3500,binwidth=60)
_,ns,_,_=histrv(up[:,2],0,3500,binwidth=60)

bar(0:60:3500,[nsg ns],orientation=:h,bar_position=:stack,xlims=[-0.3,10])

plot([nsg ns],range(0,step=40,length=length(ns)),fillalpha=0.4,fill=true,labels=["SU" "SU with STA"],linealpha=0,xlims=[-0.3,7])
hline!([layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(-0.3,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)

savefig("unit_layer_dist.svg")



import Base:sin
sin(x,y;fx=1,fy=1) = sin(2π*(fx*x+fy*y))



t = dataset["ulsta"][92][:,:,10,2]

Plots.heatmap(t)

ft=abs.(fft(t))
Plots.heatmap(ft)

freqs = fftfreq(25,30)


lmpi = [Tuple(argmax(clc))...]

lmi = clc[:,:,lmpi[3:end]...]

heatmap(lmi,aspect_ratio=:equal,frame=:all,color=:coolwarm,yflip=true)
scatter!([lmpi[2]],[lmpi[1]],markersize=10)
seg = seeded_region_growing(lmi,[(CartesianIndex(lmpi[1:2]...),2),(CartesianIndex(1,1),1)])

lr = labels_map(seg)
heatmap(lr,aspect_ratio=:equal,frame=:all,color=:grays)

ri = mapreduce(i->[i...],hcat, Tuple.(findall(lr.>1)))


t = extrema(ri,dims=2)
tr = map(i->i[2]-i[1],t)
rsize = floor(Int,maximum(tr)/2)

lp = floor.(Int,mean.(t))[:]

hspan!([lp[1]-rsize,lp[1]+rsize],alpha=0.2)
vspan!([lp[2]-rsize,lp[2]+rsize],alpha=0.2)


@manipulate for u in uids
t1 = stas[1]["usta"][u]
t2=stas[2]["usta"][u]
t3=stas[3]["usta"][u]
t4=stas[4]["usta"][u]

sds = Array{Float64}(undef,imagesize...,41)
for d in 1:41
t = t1[d,:]
tt = fill(mean(t),imagesize)
tt[xi]=t
sds[:,:,d]= mapwindow(std,tt,(5,5))
end

mi = [Tuple(argmax(sds))...][1:2]
display(mi)
mir= map(i->filter(j->j>0,(-5:5) .+i),mi)
display(mir)
sds1 = [sum(sds[mir...,d]) for d in 1:41]
# bm1 = mean(t1[1,:])
# display(bm1)
# bsd1 = std(t1[1,:])
# display(bsd1)
# bm2 = mean(t2[1,:])
# bsd2 = std(t2[1,:])
# bm3 = mean(t3[1,:])
# bsd3 = std(t3[1,:])
# bm4 = mean(t4[1,:])
# bsd4 = std(t4[1,:])
#
#
# sds1=sum((t1.-bm1)./bsd1.^2,dims=2)[:]
# sds2=sum((t2.-bm2)./bsd2.^2,dims=2)[:]
# sds3=sum((t3.-bm3)./bsd3.^2,dims=2)[:]
# sds4=sum((t4.-bm4)./bsd4.^2,dims=2)[:]

plot(sds1)


# c = [:black :red :lightgreen :blue]
# plot([sds1 sds2 sds3 sds4],color=c,lw=2,grid=false)
# a=[sds1[[1,2,end,end-1]],sds2[[1,2,end,end-1]],sds3[[1,2,end,end-1]],sds4[[1,2,end,end-1]]]'
# hline!(mean.(a).+5*std.(a),color=c)
end


skm = pyimport("skimage.metrics")

@manipulate for u in uids
    t1 = stas[1]["usta"][u][1,:]
    t2 = stas[1]["usta"][u][25,:]
    ti1=fill(mean(t1),imagesize)
    ti1[xi]=t1
    ti2=fill(mean(t2),imagesize)
    ti2[xi]=t2
    ssim,di=skm.structural_similarity(ti1,ti2,full=true)

    heatmap(di,clims=(0.4,1))
end
