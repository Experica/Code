using NeuroAnalysis,Statistics,DataFrames,StatsPlots,Mmap,Images,StatsBase,Interact


# Prepare Dataset
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"
settimeunit(1)

subject = "AE9";recordsession = "";recordsite = "u004"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
sitedir = joinpath(resultroot,subject,siteid)
testid = "$(siteid)_013"
resultdir = joinpath(sitedir,testid)
isdir(resultdir) || mkpath(resultdir)
dataset = prepare(joinpath(dataexportroot,"$testid.mat"))

# Condition Test
ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
spike = dataset["spike"]
eval.([:($(Symbol(k))=spike[$k]) for k in keys(spike)])
condon = ex["CondTest"]["CondOn"]
condoff = ex["CondTest"]["CondOff"]
minconddur=minimum(condoff-condon)
histogram(condoff-condon,nbins=20,title="Condition Duration(Set to $conddur)")


# Unit Spike Trian
epochext = preicidur
@manipulate for u in 1:length(unitspike)
ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
plotspiketrain(ys,timeline=[0,minconddur],title="Unit_$(unitid[u])")
end

for u in 1:length(unitspike)
ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
plotspiketrain(ys,timeline=[0,minconddur],title="Unit_$(unitid[u])")
foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_SpikeTrian$i")),[".png",".svg"])
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
responsedelay=0.015
@manipulate for u in 1:length(unitspike)
rs = subrvr(unitspike[u],ccondon.+responsedelay,ccondoff.+responsedelay)
plotcondresponse(rs,cctc)
end

for u in 1:length(unitspike)
rs = subrvr(unitspike[u],ccondon.+responsedelay,ccondoff.+responsedelay)
plotcondresponse(rs,cctc)
foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_CondResponse$i")),[".png",".svg"])
end

if !isempty(bcondon)
    ubr=[]
    for u in 1:length(unitspike)
        br = subrvr(unitspike[u],bcondon.+responsedelay,bcondoff.+responsedelay)
        mseuc = condresponse(br,condin(bctc))
        push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
    end
end

# Condition Response in Factor Space
fms,fses,fa=factorresponse(unitspike,cctc,ccondon,ccondoff,responsedelay=responsedelay)

@manipulate for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
    fd = findfirst(f.==keys(fa))
    fdn = length(fa[f])
    p[fd]=1:fdn
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(unitid[u],fdn))
    mseuc[f]=fa[f]
    plotcondresponse(dropmissing(mseuc),colors=[:black],responseline=ubr[u:u])
end

ufs = Dict(k=>[] for k in keys(fa))
for u in 1:length(fms), f in collect(keys(fa))
    p = Any[Tuple(argmax(coalesce.(fms[u],-Inf)))...]
    fd = findfirst(f.==keys(fa))
    fdn = length(fa[f])
    p[fd]=1:fdn
    mseuc=DataFrame(m=fms[u][p...],se=fses[u][p...],u=fill(unitid[u],fdn))
    mseuc[f]=fa[f]

    push!(ufs[f],factorresponsestats(mseuc[f],mseuc[:m],factor=f))

    # plotcondresponse(dropmissing(mseuc),colors=[:black],responseline=ubr[u:u])
    # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_$(f)_Tuning$i")),[".png",".svg"])
end


# Tuning map
plotunitposition(unitposition,color=map(i->HSV(2*i.oo,1,1-i.ocv),ufs[:Ori]),alpha=1)
foreach(i->savefig(joinpath(resultdir,"UnitPosition_OriTuning$i")),[".png",".svg"])
plotunitposition(unitposition,color=map(i->HSV(i.od,1,1-i.dcv),ufs[:Ori]),alpha=1)
foreach(i->savefig(joinpath(resultdir,"UnitPosition_DirTuning$i")),[".png",".svg"])
save(joinpath(resultdir,"factorresponse.jld2"),"factorstats",ufs,"fms",fms,"fses",fses,"fa",fa)


@df DataFrame(ufs[:Ori]) corrplot([:oo :od :ocv :dcv],nbins=30)
