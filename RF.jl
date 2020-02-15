using NeuroAnalysis,Statistics,StatsBase,FileIO,Images,Plots,Interact,ImageSegmentation,LsqFit,Makie,FFTW,ProgressMeter

dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL3"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
resultsitedir = joinpath(resultroot,subject,siteid)
layer = load(joinpath(resultsitedir,"layer.jld2"),"layer")



function exd(stas;sig=true)
    exi = [argmin(stas),argmax(stas)]
    ex = stas[exi]
    absexi = argmax(abs.(ex))
    (ex=sig ? ex[absexi] : 0,d=exi[absexi][3])
end
function joinsta(stas)
    cn = length(stas)
    imagesize = stas[1]["imagesize"]
    delays = stas[1]["delays"]
    stisize = stas[1]["stisize"]
    ppd = imagesize[1]/stisize[1]
    xi = stas[1]["xi"]
    nxi = setdiff(1:prod(imagesize),xi)
    cii=CartesianIndices(imagesize)
    bdi = filter(i->!isnothing(i),indexin([-20,-15,-10,-5,0],delays))

    dataset = Dict("stisize"=>stisize,"imagesize"=>imagesize,"delays"=>delays,"ppd"=>ppd,"bdi"=>bdi)
    dataset["log"] = [i["log"] for i in stas]
    dataset["color"]=[i["color"] for i in stas]
    dataset["minmaxcolor"] = [(i["mincolor"],i["maxcolor"]) for i in stas]
    dataset["minmaxcg"] = map(i->minmaxcolorgradient(i...,n=256),dataset["minmaxcolor"])
    usta = Dict();ucex=Dict()

    uids = mapreduce(c->keys(stas[c]["usta"]),intersect,1:cn)
    for u in uids
        csta = Array{Float64}(undef,imagesize...,length(delays),cn)
        bsta = mapreduce(c->stas[c]["usta"][u][bdi,:],vcat,1:cn)
        bm = mean(bsta);bsd = std(bsta)
        cex = []
        for c in 1:cn
            zsta = zscore(stas[c]["usta"][u],bm,bsd)
            for d in eachindex(delays)
                csta[cii[nxi],d,c].=mean(zsta[d,:])
                csta[cii[xi],d,c] = zsta[d,:]
            end
            push!(cex,exd(csta[:,:,:,c]))
        end
        usta[u] = csta
        ucex[u] = cex
    end

    dataset["usta"] = usta
    dataset["ucex"] = ucex
    return dataset
end
"Local RMS Contrast of each image"
function localcontrast(csta;ws=0.5,ppd=50)
    dims = size(csta)
    clc = Array{Float64}(undef,dims)
    w = floor(Int,ws*ppd)
    w = iseven(w) ? w+1 : w
    for c in 1:dims[4], d in 1:dims[3]
        clc[:,:,d,c] = mapwindow(std,csta[:,:,d,c],(w,w))
    end
    return clc
end
function peakroi(clc)
    pi = [Tuple(argmax(clc))...]
    plc = clc[:,:,pi[3:end]...]
    segs = seeded_region_growing(plc,[(CartesianIndex(pi[1:2]...),2),(CartesianIndex(1,1),1)])
    idx = findall(labels_map(segs).==2)
    idxlim = dropdims(extrema(mapreduce(i->[Tuple(i)...],hcat,idx),dims=2),dims=2)
    hw = map(i->i[2]-i[1],idxlim)
    radius = floor(Int,maximum(hw)/2)
    center = round.(Int,mean.(idxlim))
    cpd = [Tuple(argmax(clc[:,:,:,c]))[3] for c in 1:size(clc,4)]
    return idx,center,radius,pi[4],cpd
end
import NeuroAnalysis: isresponsive
function isresponsive(csta,idx,d,bdi;msdfactor=3,csdfactor=3)
    d in bdi && return false
    blc = [std(csta[idx,j]) for j in bdi]
    bcm = mean(blc);bcsd=std(blc)
    blm = [mean(csta[idx,j]) for j in bdi]
    bmm = mean(blm);bmsd=std(blm)
    (std(csta[idx,d]) > bcm+csdfactor*bcsd) || (mean(csta[idx,d]) > bmm+msdfactor*bmsd)
end
function responsivesta!(dataset;ws=0.5,msdfactor=3.5,csdfactor=3.5,roimargin=0.2,roirlim=1.5)
    imagesize = dataset["imagesize"]
    ppd=dataset["ppd"]
    bdi = dataset["bdi"]
    usta = dataset["usta"]
    ucex = dataset["ucex"]

    uresponsive=Dict();ulsta=Dict();ulroi=Dict()
    p = ProgressMeter.Progress(length(usta),desc="Test STAs ... ")
    for u in keys(usta)
        clc = localcontrast(usta[u],ws=ws,ppd=ppd)
        idx,center,radius,c,cpd = peakroi(clc)
        if radius > roirlim*ppd
            uresponsive[u] = false
            next!(p)
            continue
        end
        rs= [isresponsive(usta[u][:,:,:,i],idx,cpd[i],bdi,msdfactor=msdfactor,csdfactor=csdfactor) for i in 1:size(clc,4)]
        if any(rs)
            uresponsive[u] = true
            radius = ceil(Int,radius*(1+roimargin))
            vr = map(i->intersect(i.+(-radius:radius),1:imagesize[1]),center)
            radius = floor(Int,minimum(map(r->(r[end]-r[1])/2,vr)))
            ulsta[u] = usta[u][map(i->i.+(-radius:radius),center)...,:,:]
            ulroi[u] = (center=center,stisize=(2radius+1)/ppd,d=cpd[c],c=c)

            ucex[u]=map(i->exd(ulsta[u][:,:,:,i],sig=rs[i]),1:size(ulsta[u],4))
        end
        next!(p)
    end
    dataset["uresponsive"]=uresponsive
    dataset["ulsta"]=ulsta
    dataset["ulroi"]=ulroi

    return dataset
end
rfgabor(x,y,p...) = gaborf(x,y,a=p[1],μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[3]*p[5],θ=p[6],f=p[7],phase=p[8])
rfdog(x,y,p...) = dogf(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3]*p[5],θₑ=p[6],aᵢ=p[7],μᵢ₁=p[8],σᵢ₁=p[9],μᵢ₂=p[10],σᵢ₂=p[9]*p[11],θᵢ=p[12])
function modelfit(data,ppd;model=:gabor)
    dmin,dmax = abs.(extrema(data))
    dlim = maximum([dmin,dmax])
    pr = (size(data)[1]-1)/2
    sr = pr/ppd

    x = (mapreduce(i->[i[2] -i[1]],vcat,CartesianIndices(data)) .+ [-pr pr])/ppd
    y = vec(data)
    if model == :dog
        lb=[0.5dmax,   -0.4sr,    0.1sr,   -0.4sr,    0.5,    0,     0.5dmin,     -0.4sr,     0.2sr,    -0.4sr,    0.5,    0]
        ub=[1.5dmax,    0.4sr,    0.4sr,    0.4sr,    2,      π,     1.5dmin,      0.4sr,     0.5sr,     0.4sr,    2,      π]
        p0=[dmax,       0,        0.2sr,    0,        1,      0,     dmin,         0,         0.3sr,     0,        1,      0]

        mfit = curve_fit((x,p)->rfdog.(x[:,1],x[:,2],p...),x,y,
        p0,lower=lb,upper=ub,maxIter=2000,x_tol=1e-10,g_tol=1e-14)
    else
        lb=[0.5dlim,  -0.4sr,   0.1sr,   -0.4sr,   0.2,    0,    0.4/sr,   0]
        ub=[1.5dlim,   0.4sr,   0.6sr,    0.4sr,   3,      π,    3/sr,     1]
        p0=[dlim,      0,       0.2sr,    0,       1,      0,    1/sr,     0]

        mfit = curve_fit((x,p)->rfgabor.(x[:,1],x[:,2],p...),x,y,
        p0,lower=lb,upper=ub,maxIter=2000,x_tol=1e-10,g_tol=1e-14,min_step_quality=1e-4,good_step_quality=0.5,lambda_increase=5,lambda_decrease=0.5)
    end
end
function rffit!(dataset;model=[:gabor,:dog])
    ulsta = dataset["ulsta"]
    ulroi = dataset["ulroi"]
    ppd = dataset["ppd"]

    urf=Dict()
    p = ProgressMeter.Progress(length(ulsta),desc="Fit RFs ... ")
    for u in keys(ulsta)
        exs=map(i->i.ex,dataset["ucex"][u])
        ds=map(i->i.d,dataset["ucex"][u])
        urf[u] = Dict(m=>map(i->exs[i]==0 ? missing : modelfit(ulsta[u][:,:,ds[i],i],ppd,model=m),1:length(exs)) for m in model)
        next!(p)
    end

    dataset["urf"]=urf
    return dataset
end



# Load all stas
testids = ["$(siteid)_HartleySubspace_$i" for i in 1:4]
dataset = joinsta(load.(joinpath.(resultsitedir,testids,"sta.jld2")))
dataset = responsivesta!(dataset)
dataset = rffit!(dataset)

save(joinpath(resultsitedir,"stadataset.jld2"),"dataset",dataset)
dataset = load(joinpath(resultsitedir,"stadataset.jld2"),"dataset")

@manipulate for u in sort(collect(keys(dataset["usta"]))),d in eachindex(dataset["delays"])
    usta=dataset["usta"][u]
    delays = dataset["delays"]
    imagesize = dataset["imagesize"]
    stisize = dataset["stisize"]
    uresponsive = dataset["uresponsive"]
    clim = maximum(map(i->abs(i.ex),dataset["ucex"][u]))
    xlim = [0,stisize[2]]
    ylim = [0,stisize[1]]
    x = range(xlim...,length=imagesize[2])
    y = range(ylim...,length=imagesize[1])

    p = Plots.plot(layout=(1,4),legend=false,size=(1600,600),title="Unit_$(u)_STA_$(delays[d])_Responsive_$(uresponsive[u])")
    foreach(c->Plots.heatmap!(p,subplot=c,x,y,usta[:,:,d,c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xlim,ylims=ylim,xticks=xlim,yticks=ylim,yflip=true,xlabel=dataset["color"][c]),1:4)
    p
end

@manipulate for u in sort(collect(keys(dataset["ulsta"])))
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    imagesize = size(ulsta)
    stisize = dataset["ulroi"][u].stisize
    ucex = map(i->i.ex,dataset["ucex"][u])
    ucd = map(i->i.d,dataset["ucex"][u])
    clim = maximum(abs.(ucex))
    xylim = [0,round(stisize,digits=1)]
    xy = range(xylim...,length=imagesize[2])

    p = Plots.plot(layout=(1,5),legend=false,size=(1650,600))
    Plots.bar!(p,subplot=1,dataset["color"],ucex,frame=:zerolines)
    foreach(c->Plots.heatmap!(p,subplot=c+1,xy,xy,ulsta[:,:,ucd[c],c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xylim,ylims=xylim,xticks=xylim,yticks=[],yflip=true,xlabel=dataset["color"][c],title="Unit_$(u)_STA_$(delays[ucd[c]])"),1:4)
    p
end

@manipulate for u in sort(collect(keys(dataset["urf"]))),m in collect(keys(first(values(dataset["urf"]))))
    ulsta = dataset["ulsta"][u]
    delays = dataset["delays"]
    imagesize = size(ulsta)
    stisize = dataset["ulroi"][u].stisize
    ucex = map(i->i.ex,dataset["ucex"][u])
    ucd = map(i->i.d,dataset["ucex"][u])
    clim = maximum(abs.(ucex))
    xylim = [0,round(stisize,digits=1)]
    xy = range(xylim...,length=imagesize[2])

    p = Plots.plot(layout=(2,4),legend=false,size=(1500,650))
    foreach(c->Plots.heatmap!(p,subplot=c,xy,xy,ulsta[:,:,ucd[c],c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xylim,ylims=xylim,xticks=xylim,yticks=[],yflip=true,title="Unit_$(u)_STA_$(delays[ucd[c]])"),1:4)

    fit=dataset["urf"][u]
    sr = round(stisize/2,digits=1)
    x=y=-sr:0.01:sr
    xylim=[-sr,sr]
    if m == :dog
        rfs = map(f->ismissing(f) ? missing : [rfdog(i,j,f.param...) for j in reverse(y), i in x],fit[m])
    else
        rfs = map(f->ismissing(f) ? missing : [rfgabor(i,j,f.param...) for j in reverse(y), i in x],fit[m])
    end
    foreach(c->ismissing(rfs[c]) ? Plots.plot!(p,subplot=c+4,frame=:none) :
    Plots.heatmap!(p,subplot=c+4,x,y,rfs[c],aspect_ratio=:equal,frame=:semi,color=:coolwarm,clims=(-clim,clim),
    xlims=xylim,ylims=xylim,xticks=xylim,yticks=xylim,yformatter=i->-i,yflip=true,xlabel=dataset["color"][c],
    title="Unit_$(u)_STA_$(delays[ucd[c]])_FIT_$(m)"),1:4)
    p
end


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



spike = load(joinpath(resultsitedir,testids[1],"spike.jld2"),"spike")
unitid = spike["unitid"];unitposition = spike["unitposition"];unitlayer = assignlayer(unitposition[:,2],layer)

us = collect(keys(dataset["ulsta"]))
up = unitposition[indexin(us,unitid),:]
ul = unitlayer[indexin(us,unitid)]



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
