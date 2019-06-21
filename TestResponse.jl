using NeuroAnalysis,Statistics,DataFrames,StatsPlots,Mmap,Images,StatsBase,Combinatorics,Interact


# Prepare Dataset
subject="AE9";recordsession="u004";test="010"
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
ci = ctc[factors[1]].!="blank"
cctc = ctc[ci,factors]
ccondon=condon[ci]
ccondoff=condoff[ci]

bi = .!ci
bctc = ctc[bi,factors]
bcondon=condon[bi]
bcondoff=condoff[bi]

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

if !isempty(bcondon)
    ubr=[]
    for u in 1:length(unitspike)
        br = subrvr(unitspike[u],bcondon.+delay,bcondoff.+delay)
        mseuc = condresponse(br,condin(bctc))
        push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
    end
end

# Factor Response
fms,fses,fa=factorresponse(unitspike,cctc,ccondon,ccondoff,delay=delay)

@manipulate for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
    #opt = map((i,j)->j[i],c,values(fa))
    fd = findfirst(f.==keys(fa))
    fdn = length(fa[f])
    p[fd]=1:fdn
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(unitid[u],fdn))
    mseuc[f]=fa[f]
    plotcondresponse(dropmissing(mseuc),colors=[:black],responseline=ubr[u:u])
end

for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
    fd = findfirst(f.==keys(fa))
    fdn = length(fa[f])
    p[fd]=1:fdn
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(unitid[u],fdn))
    mseuc[f]=fa[f]
    plotcondresponse(dropmissing(mseuc),colors=[:black],responseline=ubr[u:u])
    foreach(i->savefig(joinpath(resultdir,"Unit$(unitid[u])_$(f)_Tuning$i")),[".png",".svg"])
end



a=DataFrame(a=[1,2,1,2],b=[1,1,2,2])
b=(a=1,b=1)

a.==[b]

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


RGBA(0.1,0.1,0.3,0.8)

OrderedDict(map((i,j)->i=>j,1:3,4:6))
