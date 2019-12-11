using NeuroAnalysis,Statistics,DataFrames,StatsPlots,Mmap,Images,StatsBase,Interact


# Prepare Dataset
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL1";test = "Color_2"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
datadir = joinpath(dataroot,subject,siteid)
resultsitedir = joinpath(resultroot,subject,siteid)
testid = join(filter(!isempty,[siteid,test]),"_")
resultdir = joinpath(resultsitedir,testid)
isdir(resultdir) || mkpath(resultdir)
dataset = prepare(joinpath(dataexportroot,subject,"$testid.mat"))

# Condition Test
ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
condon = ex["CondTest"]["CondOn"]
condoff = ex["CondTest"]["CondOff"]
minconddur=minimum(condoff-condon)
histogram(condoff-condon,nbins=20,title="Condition Duration(Set to $conddur)")


# Unit Spike Trian segments for all condtests
spike = dataset["spike"]
eval.([:($(Symbol(k))=spike[$k]) for k in keys(spike)])

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

## Prepare Conditions
ctc = condtestcond(ex["CondTestCond"])
cond = condin(ctc)
factors = finalfactor(ctc)
# for color test
factors = [:HueAngle]
blank = (:Color,36)

ci = ctc[!,blank[1]].!=blank[2]
cctc = ctc[ci,factors]
ccond = condin(cctc)
ccondon=condon[ci]
ccondoff=condoff[ci]

bi = .!ci
bctc = ctc[bi,factors]
bcondon=condon[bi]
bcondoff=condoff[bi]
isblank = !isempty(bcondon)

## Condition Response
responsedelay=15
@manipulate for u in 1:length(unitspike)
rs = subrvr(unitspike[u],ccondon.+responsedelay,ccondoff.+responsedelay)
plotcondresponse(rs,cctc)
end

ubr=[];uresponsive=[];umodulative=[]
for u in 1:length(unitspike)
    rs = subrvr(unitspike[u],ccondon.+responsedelay,ccondoff.+responsedelay)
    basers = subrvr(unitspike[u],ccondon.+(responsedelay-preicidur),ccondon.+responsedelay)
    push!(uresponsive,isresponsive(basers,rs,ccond.i))
    push!(umodulative,ismodulative([DataFrame(Y=rs) cctc]))
    if isblank
        br = subrvr(unitspike[u],bcondon.+responsedelay,bcondoff.+responsedelay)
        mseuc = condresponse(br,condin(bctc))
        push!(ubr,(m=mseuc[1,:m],se=mseuc[1,:se]))
    end
    # plotcondresponse(rs,cctc)
    # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_CondResponse$i")),[".png",".svg"])
end

## Condition Response in Factor Space
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
    mseuc[!,f]=fa[f]

    push!(ufs[f],factorresponsestats(mseuc[!,f],mseuc[!,:m],factor=f))

    plotcondresponse(dropmissing(mseuc),colors=[:black],responseline=isblank ? ubr[u:u] : [])
    foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_$(f)_Tuning$i")),[".png",".svg"])
end
save(joinpath(resultdir,"factorresponse.jld2"),"factorstats",ufs,"fms",fms,"fses",fses,"fa",fa,"responsive",uresponsive,"modulative",umodulative)

## Tuning map
plotunitposition(unitposition,color=map((i,j)->j ? HSV(1*i.oo,1-i.ocv,1) : HSVA(0,0,0.5,0.2),ufs[:Ori],umodulative))
foreach(i->savefig(joinpath(resultdir,"UnitPosition_OriTuning$i")),[".png",".svg"])
plotunitposition(unitposition,color=map((i,j)->j ? HSV(i.od,1-i.dcv,1) : HSVA(0,0,0.5,0.2),ufs[:Ori],umodulative))
foreach(i->savefig(joinpath(resultdir,"UnitPosition_DirTuning$i")),[".png",".svg"])
plotunitposition(unitposition,color=map((i,j)->j ? HSV(i.oh,1-i.hcv,1) : HSVA(0,0,0.5,0.2),ufs[:ColorID],umodulative))
foreach(i->savefig(joinpath(resultdir,"UnitPosition_HueTuning$i")),[".png",".svg"])

@df DataFrame(ufs[:Ori]) corrplot([:oo :od :ocv :dcv],nbins=30)

# for color test
vi = uresponsive.&umodulative
colorspace = ex["Param"]["ColorSpace"]
hues = ex["Param"]["Color"]
ivc = colorspace == "HSL" ? RGBA(0.5,0.5,0.5,1) : RGBA(cond[1,:BGColor]...)
title = "UnitPosition_$(colorspace)_$(hues)_OptimalHue"
plotunitposition(unitposition,color=map((i,j)->j ? RGBA(cond[cond[:HueAngle].==i.oh,:Color][1]...) : ivc,ufs[:HueAngle],vi),title=title)
foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])

ys,ns,ws,is = histrv(map(i->i.oh,ufs[:HueAngle][vi]),0,360,nbins=10)
title = "$(colorspace)_$(hues)_OptimalHue_Distribution"
plot(deg2rad.(mean.(ws)),ns,title=title,projection=:polar)
foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])

scond = sort(cond,:HueAngle)
map(i->RGB(i[1:3]...),scond[:Color])
