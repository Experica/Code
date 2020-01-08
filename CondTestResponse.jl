using NeuroAnalysis,Statistics,DataFrames,StatsPlots,Mmap,Images,StatsBase,Interact


# Prepare Dataset
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"
imageroot = "../NaturalStimuli"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL3";test = "HartleySubspace_1"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
datadir = joinpath(dataroot,subject,siteid)
resultsitedir = joinpath(resultroot,subject,siteid)
testid = join(filter(!isempty,[siteid,test]),"_")
resultdir = joinpath(resultsitedir,testid)
isdir(resultdir) || mkpath(resultdir)
dataset = prepare(joinpath(dataexportroot,subject,"$testid.mat"))

# Condition Test
ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
spike = dataset["spike"]
eval.([:($(Symbol(k))=spike[$k]) for k in keys(spike)])
condon = ex["CondTest"]["CondOn"]
condoff = ex["CondTest"]["CondOff"]
condidx = ex["CondTest"]["CondIndex"]
minconddur=minimum(condoff-condon)
histogram(condoff-condon,nbins=20,title="Condition Duration(Set to $conddur)")

ugs = map(i->i ? "Single-" : "Multi-",unitgood)
# Unit Position
plotunitposition(spike,layer=layer)
foreach(i->savefig(joinpath(resultdir,"Unit_Position$i")),[".png",".svg"])

# Unit Spike Trian segments for all condtests
epochext = preicidur
@manipulate for u in 1:length(unitspike)
ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
plotspiketrain(ys,timeline=[0,minconddur],title="Unit_$(unitid[u])")
end

for u in 1:length(unitspike)
    ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
    plotspiketrain(ys,timeline=[0,minconddur],title="$(ugs[u])Unit_$(unitid[u])")
    foreach(i->savefig(joinpath(resultdir,"$(ugs[u])Unit_$(unitid[u])_SpikeTrian$i")),[".png",".svg"])
end

## Prepare Conditions
ctc = condtestcond(ex["CondTestCond"])
cond = condin(ctc)
factors = finalfactor(ctc)
# for color test
# factors = [:HueAngle]
# blank = (:Color,36)

# for orisf test
# blank = (:SpatialFreq,0)

# for hartley test
blank = (:Ori_Final,NaN)

ci = ctc[!,blank[1]].!==blank[2]
cctc = ctc[ci,factors]
ccond = condin(cctc)
ccondon=condon[ci]
ccondoff=condoff[ci]
ccondidx = condidx[ci]

bi = .!ci
bctc = ctc[bi,factors]
bcondon=condon[bi]
bcondoff=condoff[bi]
bcondidx = condidx[bi]
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

    # plotcondresponse(dropmissing(mseuc),colors=[:black],projection=f==:Ori_Final ? :polar : :cartesian,responseline=isblank ? ubr[u:u] : [])
    # foreach(i->savefig(joinpath(resultdir,"Unit_$(unitid[u])_$(f)_Tuning$i")),[".png",".svg"])
end
save(joinpath(resultdir,"factorresponse.jld2"),"factorstats",ufs,"fms",fms,"fses",fses,"fa",fa,"responsive",uresponsive,"modulative",umodulative)

## Tuning map
vi = uresponsive.&umodulative.&unitgood
# for orisf test
title = "UnitPosition_OptimalOri"
plotunitposition(unitposition,color=map((i,j)->j ? HSV(i.oo,1-i.ocv/2,1) : RGBA(0.5,0.5,0.5,0.1),ufs[:Ori_Final],vi),title=title)
foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
title = "UnitPosition_OptimalDir"
plotunitposition(unitposition,color=map((i,j)->j ? HSV(i.od,1-i.dcv/2,1) : RGBA(0.5,0.5,0.5,0.1),ufs[:Ori_Final],vi),title=title)
foreach(i->savefig(joinpath(resultdir,"UnitPosition_DirTuning$i")),[".png",".svg"])
plotunitposition(unitposition,color=map((i,j)->j ? HSV(i.oh,1-i.hcv,1) : HSVA(0,0,0.5,0.2),ufs[:ColorID],umodulative))
foreach(i->savefig(joinpath(resultdir,"UnitPosition_HueTuning$i")),[".png",".svg"])

@df DataFrame(ufs[:Ori]) corrplot([:oo :od :ocv :dcv],nbins=30)

# for color test
colorspace = ex["Param"]["ColorSpace"]
hues = ex["Param"]["Color"]
title = "UnitPosition_$(colorspace)_$(hues)_OptimalHue"
plotunitposition(unitposition,color=map((i,j)->j ? RGBA(clamp.(cond[findclosestangle(deg2rad(i.oh),deg2rad.(cond[!,:HueAngle])),:Color],0,1)...) : RGBA(0.5,0.5,0.5,0.1),ufs[:HueAngle],vi),layer=layer,title=title)
foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])

ys,ns,ws,is = histrv(map(i->i.oh,ufs[:HueAngle][vi]),0,360,nbins=10)
title = "$(colorspace)_$(hues)_OptimalHue_Distribution"
plot(deg2rad.(mean.(ws)),ns,title=title,projection=:polar)
foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])



## Prepare Images
imagesetname = splitext(splitdir(ex["CondPath"])[2])[1]
imageset = map(i->GrayA.(i),prepare(joinpath(imageroot,"$(imagesetname).mat"))["imgs"])
@manipulate for i in 1:length(imageset)
    imageset[i]
end

ppd = size(imageset[1],1)/5


getparam(envparam,"Diameter")





bgcolor = oftype(imageset[1][1],RGBA(getparam(envparam,"BGColor")...))
masktype = getparam(envparam,"MaskType")
maskradius = getparam(envparam,"MaskRadius")
masksigma = getparam(envparam,"Sigma")

nscale = 2;downsample = 2;sigma = 1.5
imagestimuli = map(i->gaussian_pyramid(i, nscale-1, downsample, sigma),imageset)
imagesize = map(i->size(i),imagestimuli[1])
unmaskindex = map(i->alphamask(i,radius=maskradius,sigma=masksigma,masktype=masktype)[2],imagestimuli[1])
imagestimuli = map(s->map(i->alphablend.(alphamask(i[s],radius=maskradius,sigma=masksigma,masktype=masktype)[1],[bgcolor]),imagestimuli),1:nscale)
@manipulate for i in 1:length(imageset),s in 1:length(imagestimuli)
    imagestimuli[s][i]
end


scaleindex=1
uii = unique(ccondidx)
uci = map(i->findall(ccondidx.==i),uii)

x = Array{Float64}(undef,length(uii),prod(imagesize[scaleindex]))
foreach(i->x[i,:]=vec(gray.(imagestimuli[scaleindex][uii[i]])),1:size(x,1))

@manipulate for u in 1:length(unitspike), d in 0:20:40
    y = subrvr(unitspike[u],ccondon.+d,ccondoff.+d,israte=false)
    y = map(i->sum(y[i]),uccii)
    plotsta(sta(x,y),imagesize[imagescale],delay=d)
end

u=155
d=60
y = subrvr(unitspike[u],ccondon.+d,ccondoff.+d,israte=false)
y = map(i->sum(y[i]),uci)
plotsta(sta(x,y),imagesize[scaleindex],delay=d,title="$(ugs[u])Unit_$(unitid[u])_STA_$(d)")


for u in 1:length(unitspike), d in 0.04:0.01:0.12
    _,y,_,_ = subrv(unitspike[u],ccondon.+d,ccondoff.+d,israte=false)
    plotsta(sta(x,y),imagesize,delay=d)
    foreach(i->savefig(joinpath(resultdir,"Unit$(unitid[u])_STA_$d$i")),[".png",".svg"])
end




om=[(1,0.5),(-1,0.5),(-1,-0.5),(1,-0.5)]
polygon = om

scatter(rand(10),rand(10),marker=(Shape(polygon),10,:red))
