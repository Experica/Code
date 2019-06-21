using NeuroAnalysis,Statistics,DataFrames,StatsPlots,Images,StatsBase,Combinatorics,Interact


# Prepare Dataset
subject="AE9";recordsession="u004";test="005"
sessionid = join([subject,recordsession],"_")
testid = join([subject,recordsession,test],"_")

dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"
resultdir = joinpath(resultroot,testid)
isdir(resultdir) || mkdir(resultdir)
lfregex = Regex("^$(uppercase(testid))[A-Za-z0-9_]*.imec.lf.bin")
lffile = joinpath(dataroot,filter(i->occursin(lfregex,i),readdir(dataroot))[1])
dataset = prepare(joinpath(dataexportroot,"$testid.mat"))

## Condition Tests
ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
spike = dataset["spike"]
eval.([:($(Symbol(k))=spike[$k]) for k in keys(spike)])
condon = ex["CondTest"]["CondOn"]
condoff = ex["CondTest"]["CondOff"]
condidx = ex["CondTest"]["CondIndex"]
minconddur=minimum(condoff-condon)
histogram(condoff-condon,bins=20,title="Condition Duration(Set to $conddur)",legend=false)


# Unit Spike Trian
epochext = preicidur
@manipulate for u in 1:length(unitspike)
ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
plotspiketrain(ys,timeline=[0,minconddur],title="Unit $(unitid[u])")
end

for u in 1:length(unitspike)
ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
plotspiketrain(ys,timeline=[0,minconddur],title="Unit $(unitid[u])")
foreach(i->savefig(joinpath(resultdir,"Unit$(unitid[u])_SpikeTrian$i")),[".png",".svg"])
end

# Prepare Conditions
ctc = condtestcond(ex["CondTestCond"])
factors = finalfactor(ctc)
ci = .!isnan.(ctc[factors[1]])
cctc = ctc[ci,factors]
ccondidx = condidx[ci]
ccondon=condon[ci]
ccondoff=condoff[ci]

# Condition Response
delay=0.02
@manipulate for u in 1:length(unitspike)
rs = subrvr(unitspike[u],ccondon.+delay,ccondoff.+delay)
plotcondresponse(rs,cctc,factors)
end

for u in 1:length(unitspike)
rs = subrvr(unitspike[u],ccondon.+delay,ccondoff.+delay)
plotcondresponse(rs,cctc,factors)
foreach(i->savefig(joinpath(resultdir,"Unit$(unitid[u])_CondResponse$i")),[".png",".svg"])
end

# Factor Response
fms,fses,fa=factorresponse(unitspike,cctc,ccondon,ccondoff,delay=delay)

@manipulate for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
    fd = findfirst(f.==keys(fa))
    fdn = length(fa[f])
    p[fd]=1:fdn
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(unitid[u],fdn))
    mseuc[f]=fa[f]
    plotcondresponse(dropmissing(mseuc),colors=[:black])
end

for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
    fd = findfirst(f.==keys(fa))
    fdn = length(fa[f])
    p[fd]=1:fdn
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(unitid[u],fdn))
    mseuc[f]=fa[f]
    plotcondresponse(dropmissing(mseuc),colors=[:black])
    foreach(i->savefig(joinpath(resultdir,"Unit$(unitid[u])_$(f)_Tuning$i")),[".png",".svg"])
end



oris=float.(cctc[:Ori])
udf=[]
for u in 1:length(unitspike),f in factors
r = subrvr(unitspike[u],ccondon.+delay,ccondoff.+delay)
stats = statsori(oris,float.(r))
push!(udf,DataFrame(stats))
# plotcondresponse(r,cctc,f)
# foreach(i->savefig(joinpath(resultdir,"Unit$(u)_$(f)_Tuning")),[".png",".svg"])
end

udf = vcat(udf...)
@df udf corrplot([:dcv :ocv :pdir :pori],nbins=30)
foreach(i->savefig(joinpath(resultdir,"ori$i")),[".png",".svg"])
save(joinpath(resultdir,"ori.jld2"),"ori",udf)






# Prepare Images
imgfile = joinpath(dataroot,"$(testid)_image.mat")
imgs = readmat(imgfile)["imageArray"]
imgs=mat2julia!(imgs)
#px=Int(envparam["x_pos"]);py=Int(envparam["y_pos"])
px=613;py=135;rp=24


imageset = map(i->dropdims(mean(i[py-rp:py+rp,px-rp:px+rp,:],dims=3),dims=3),imgs)
imagesize = size(imageset[1])

@manipulate for i in 1:length(imageset)
    heatmap(imageset[i],yflip=true,aspectratio=:equal,color=:grays)
end


x = Array{Float64}(undef,length(ccondidx),prod(imagesize))
foreach(i->x[i,:]=vec(imageset[ccondidx[i]]),1:size(x,1))

@manipulate for u in 1:length(unitspike), d in 0.02:0.01:0.14
    _,y,_,_ = subrv(unitspike[u],ccondon.+d,ccondoff.+d,israte=false)
    plotsta(sta(x,y),imagesize,delay=d,filter=Kernel.gaussian(s))
end


for u in 1:length(unitspike), d in 0.04:0.01:0.12
    _,y,_,_ = subrv(unitspike[u],ccondon.+d,ccondoff.+d,israte=false)
    plotsta(sta(x,y),imagesize,delay=d)
    foreach(i->savefig(joinpath(resultdir,"Unit$(unitid[u])_STA_$d$i")),[".png",".svg"])
end
