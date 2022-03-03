using NeuroAnalysis,Statistics,DataFrames,StatsPlots,StatsBase

function collectcondtest(indir;unit=DataFrame(id=[]),ccode=ccode,datafile="factorresponse.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            fr = load(joinpath(root,datafile))
            siteid = fr["siteid"]
            exid = fr["exid"]
            id = map((u,g)->g ? "$(siteid)_SU$u" : "$(siteid)_MU$u",fr["unitid"],fr["unitgood"])
            c = ccode[fr["color"]]
            ks = ["responsive","modulative","enoughresponse","frf","fa"]
            isempty(fr["f1f0"]) || push!(ks,"f1f0")
            df1 = ("$(k)!$(exid)!$c"=>fr[k] for k in ks)
            # df2 = ("frf_$(k)!$c"=>fr["factorresponsefeature"][k] for k in keys(fr["factorresponsefeature"]))
            # df3 = ("pzfr_$(k)!$c"=>map((i,p,m)->(m[i...].-mean(p[i...]))/std(p[i...]),fr["optfri"][k],fr["pfms"],fr["fms"]) for k in keys(fr["fa"]))
            # df4 = ("f_$(k)!$c"=>fill(fr["fa"][k],length(id)) for k in keys(fr["fa"]))
            # df = DataFrame(df1...,df2...,df3...,df4...)
            df = DataFrame(df1...)
            df.siteid .= siteid
            df.id = id
            unit = outerjoin(unit,df,on=:id)
        end
    end
    return unit
end


## CondTests
resultroot = "Z:/"
figfmt = [".png"]

staunit = collectcondtest(resultroot)

## Tuning map
vi = uresponsive.&umodulative.&unitgood
colorspace = ex["Param"]["ColorSpace"]
hues = ex["Param"]["Color"]
# for orisf test
title = "UnitPosition_OptimalOriDirSF"
p=plotunitpositionproperty(unitposition[vi,:],title=title,ori=map(i->i.oo,ufs[:Ori_Final][vi]),os=map(i->1-i.ocv,ufs[:Ori_Final][vi]),
dir = map(i->i.od,ufs[:Ori_Final][vi]),ds=map(i->1-i.dcv,ufs[:Ori_Final][vi]),sf=map(i->i.osf,ufs[:SpatialFreq][vi]),layer=layer)

save("testup.svg",p)
pfactor=intersect([:Ori,:Ori_Final],keys(ufs))

oos = map(i->i.oo,ufs[:Ori_Final][vi])
oos = [oos;oos.+180]
ys,ns,ws,is = histrv(oos,0,360,nbins=20)
oa = deg2rad.(mean.([ws;ws[1]]))
oh = [ns;ns[1]]./2

ys,ns,ws,is = histrv(map(i->i.od,ufs[:Ori_Final][vi]),0,360,nbins=20)
da = deg2rad.(mean.([ws;ws[1]]))
dh = [ns;ns[1]]

title = "$(colorspace)_$(hues)_OptOriDir_Distribution"
plot([oa da],[oh dh],title=title,projection=:polar,linewidth=2,label=["Ori" "Dir"])



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

ha = ccond[!,:HueAngle]
hac = map(i->cond[findfirst(cond[!,:HueAngle].==i),:Color],ha)

plotunitposition(unitposition[vi,:],color=map(i->RGBA(clamp.(cond[findclosestangle(deg2rad(i.oh),deg2rad.(cond[!,:HueAngle])),:Color],0,1)...),ufs[:HueAngle][vi]),
markersize=map(i->(1-i.hcv)*5+2,ufs[:HueAngle][vi]),layer=layer,title=title)

foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])

ys,ns,ws,is = histrv(map(i->i.oh,ufs[:HueAngle][vi]),0,360,nbins=10)
title = "$(colorspace)_$(hues)_OptimalHue_Distribution"
plot(deg2rad.(mean.(ws)),ns,title=title,projection=:polar)
foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
