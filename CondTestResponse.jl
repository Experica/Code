using NeuroAnalysis,FileIO,Statistics,DataFrames,StatsPlots,StatsBase

# import Base: sum,getproperty
# getproperty(x::Missing,f::Symbol) = missing
# sum(x::Missing) = missing

function collectcondtest(indir;unit=DataFrame(),exid="OriSF",ccode=ccode,datafile="factorresponse.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            ctexid = load(joinpath(root,datafile),"exid")
            ctexid == exid || continue
            fr = load(joinpath(root,datafile))
            siteid = fr["siteid"]
            id = map((u,g)->g ? "$(siteid)_SU$u" : "$(siteid)_MU$u",fr["unitid"],fr["unitgood"])
            c = getccode(fr["exenv"]["color"],ccode)
            ks = ["responsive","modulative","enoughresponse","fa"]
            isempty(fr["f1f0"]) || push!(ks,"f1f0")
            kv1 = ("$(k)!$c"=>fr[k] for k in ks)
            kv2 = ("frf_$(k)!$c"=>fr["frf"][k] for k in keys(fr["frf"]))
            # df3 = ("pzfr_$(k)!$c"=>map((i,p,m)->(m[i...].-mean(p[i...]))/std(p[i...]),fr["optfri"][k],fr["pfms"],fr["fms"]) for k in keys(fr["fa"]))
            # df4 = ("f_$(k)!$c"=>fill(fr["fa"][k],length(id)) for k in keys(fr["fa"]))
            # df = DataFrame(df1...,df2...,df3...,df4...)
            df = DataFrame(kv1...,kv2...)
            df.siteid .= siteid
            df.id = id
            append!(unit,df,cols=:union)
        end
    end
    return unit
end

resultroot = "Z:/"
figfmt = [".png"]

## Collect All Color Tests
# colorunit = collectcondtest(joinpath(resultroot,"AG2","AG2_V1_ODR1"),exid="Color")
colortests = collectcondtest(resultroot,exid="Color")


names(colorunits)

dklcg = cgrad(ColorMaps["lidkl_mcchue_l0"].colors)
hslcg = cgrad(ColorMaps["hsl_mshue_l0.4"].colors)

function colorunit(unit;ccode="DKL_L0")
    select!(dropmissing(unit,"responsive!$ccode"),Not(r"\w*!\w*"),
        # ["responsive!$ccode","modulative!$ccode"] => ByRow((i,j)->i|j) => :sr,
        # ["responsive!$ccode","modulative!$ccode","enoughresponse!$ccode"] => ByRow((i,j,k)->i&j&k) => :rme,
        # ["pzfr_HueAngle!$ccode","f_HueAngle!$ccode"] => ByRow((r,f)->r[f.==0][1]-r[f.==180][1]) => :L_M,
        # ["pzfr_HueAngle!$ccode","f_HueAngle!$ccode"] => ByRow((r,f)->r[f.==90][1]-r[f.==270][1]) => :S_LM,
        "frf_Angle!$ccode" => ByRow(i->(;i.fit.mfit.model,i.fit.mfit.fun,i.fit.mfit.param)) => :fit,
        "frf_Angle!$ccode" => ByRow(i->i.up) => :up,
        "frf_Angle!$ccode" => ByRow(i->i.cm) => :cm,
        "frf_Angle!$ccode" => ByRow(i->i.cv) => :cv,
        "frf_Angle!$ccode" => ByRow(i->i.acm) => :acm,
        "frf_Angle!$ccode" => ByRow(i->i.acv) => :acv,
        "frf_Angle!$ccode" => ByRow(i->i.fit.pa) => :pa,
        "frf_Angle!$ccode" => ByRow(i->i.fit.asi2) => :asi,
        "frf_Angle!$ccode" => ByRow(i->sum(i.fit.ahw)) => :ahw,
        "frf_Angle!$ccode" => ByRow(i->i.fit.mfit.r) => :cor
        "frf_Angle!$ccode" => ByRow(i->1-i.fit.mfit.r2) => :fvu)
end

dklunit = colorunit(colorunits,ccode="DKL_L0")
hslcell = colorcell(cell,ccode="HSL_Ym")

dklcell |> Voyager()

hslcell |>
    [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"sr:n",title="Color"});
     @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"sr:n",title="Color"});
     @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"sr:n",title="Color"})]

dklcell |> @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"rme:n"})
dklcell |> @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"rme:n"})

filter(r->r.rme,dklcell) |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"L_M"},color={"layer",scale={scheme=:category10}},
            column={"site",title=""})
filter(r->r.rme,dklcell) |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"S_LM"},color={"layer",scale={scheme=:category10}},
            column={"site",title=""})

filter(r->r.rme,dklcell) |> @vlplot(:bar,transform=[{calculate="datum.L_M > 0 ? 'L+' : 'L-'","as"="Dominance"},{calculate="datum.L_M > 0 ? 1 : -1","as"="sign"}],
        y={"layer",title="Cortical Layer"},x={"sign",title="Number of rme Cells"},color={"Dominance",scale={range="#" .* hex.(get(dklcg,[0,0.5]),:rgb)}},column={"site",title=""})
filter(r->r.rme,dklcell) |> @vlplot(:bar,transform=[{calculate="datum.S_LM > 0 ? 'S+' : 'S-'","as"="Dominance"},{calculate="datum.S_LM > 0 ? 1 : -1","as"="sign"}],
        y={"layer",title="Cortical Layer"},x={"sign",title="Number of rme Cells"},color={"Dominance",scale={range="#" .* hex.(get(dklcg,[0.25,0.75]),:rgb)}},column={"site",title=""})


filter(r->r.sr,hslcell) |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"hcv"},color={"layer",scale={scheme=:category10}})
filter(r->r.rme,dklcell)
dklcell |> @vlplot(:tick,y={"layer",title="Cortical Layer"},x={"hcv"},color={"layer",scale={scheme=:category10}})

filter(r->r.rme && r.hcv <= 0.8,dklcell)
hslcell |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"oh",axis={values=0:90:360}},
            color={"layer",scale={scheme=:category10}},column={"site",title=""})

dklcell |> @vlplot(:bar,transform=[{bin={step=30},field="oh","as"="boh"}],
            y={"count()",title="Number of Cells"},x={"oh",bin={step=30},title="Optimal Hue",axis={values=0:90:360}},
            color={"boh",scale={range="#" .* hex.(get(dklcg,range(1/24,length=12,step=1/12)),:rgb)},legend=false},column={"site",title=""})

hslcell |> @vlplot(:bar,transform=[{bin={step=30},field="oh","as"="boh"}],
            y={"count()",title="Number of Cells"},x={"oh",bin={step=30},title="Optimal Hue",axis={values=0:90:360}},
            color={"boh",scale={range="#" .* hex.(get(hslcg,range(1/24,length=12,step=1/12)),:rgb)},legend=false},column={"site",title=""})

filter(r->r.rme && r.hcv <= 0.8,dklcell) |> @vlplot(:tick,y={"layer",title="Cortical Layer"},x={"oh",axis={values=0:90:360}},
            color={"layer",scale={scheme=:category10}},column={"site",title=""})







dklcs = "#" .* hex.(dklcg.colors.colors,:rgb)
hslcs = "#" .* hex.(hslcg.colors.colors,:rgb)


dklcell |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth"},x={"hcv",title="Number of Cells"},
                        color={"oh",scale={range=dklcs}},column={"site"})

dklcells |> @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"dkl_oh",scale={range=dklcs}},column={"site"})

hslcell |> @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Cells"},
                        color={"oh",scale={range=hslcs}},column={"site"})

hslcell |> @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"oh",scale={range=hslcs}},column={"site"})


plothuehist=(df;w=700,h=700,hue="oh",cs="dkl",cg=dklcg)->begin
    sites = sort(unique(df.site));n=length(sites);cshue="$(cs)_$hue"
    p=plot(layout=(n,1),legend=false,grid=true,size=(w,n*h))
    for i in 1:n
        ys,ns,ws,is = epochspiketrain(df[df.site.==sites[i],cshue],0:30:360)
        if !isnothing(cg)
            x = 2*π*cg.values
            y = fill(1.1maximum(ns),length(x))
            scatter!(p,subplot=i,x,y,projection=:polar,color=cg.colors.colors,markerstrokewidth=0,markersize=12)
        end
        plot!(p,subplot=i,deg2rad.(mean.([ws;ws[1]])),[ns;ns[1]],title="$(sites[i])--$cshue",projection=:polar,linewidth=5,linecolor=:gray30)
    end
    p
end

plothuehist(dklcells,hue="oh",cs="dkl",cg=dklcg)
foreach(ext->savefig(joinpath(resultroot,"DKLOptHueHist$ext")),figfmt)
plothuehist(hslcells,hue="oh",cs="hsl",cg=hslcg)
foreach(ext->savefig(joinpath(resultroot,"HSLOptHueHist$ext")),figfmt)




plothueposition=(;w=800,h=700)->begin
    p=plot(layout=(1,testn),legend=false,grid=false,size=(testn*w,h))
    for i in 1:testn
        xlims = extrema(upos[i][:,1]).+[-2,1]
        if !isnothing(layer)
            lx = xlims[1]+1
            hline!(p,subplot=i,[layer[k][1] for k in keys(layer)],linestyle=:dash,
            annotations=[(lx,layer[k][1],text(k,7,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray70,legend=false)
        end
        scatter!(p,subplot=i,upos[i][:,1],upos[i][:,2],title=testlogs[i],markersize=(1 .-hcvs[i])*5 .+2, color=cms[i][ohs[i]/360],
        xlims=xlims,markerstrokewidth=0)
    end
    p
end

plothueposition()
foreach(ext->savefig(joinpath(siteresultdir,"OptHue_UnitPosition$ext")),figfmt)


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
