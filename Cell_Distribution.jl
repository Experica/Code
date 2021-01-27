using NeuroAnalysis,DataFrames,VegaLite,DataVoyager,StatsBase,StatsPlots,LightGraphs,MetaGraphs,GraphRecipes,Combinatorics,
Interact,Images,UMAP,Clustering,Distances,Makie

## Load Cell database
resultroot = "../Result"
figfmt = [".png",".svg"]
cell = load(joinpath(resultroot,"cell.jld2"),"cell")

## OriSF
import Base: sum,getproperty
getproperty(x::Missing,f::Symbol) = missing
sum(x::Missing) = missing

function orisfcell(cell;ccode="A")
    select!(dropmissing(cell,"responsive!$ccode"),Not(r"\w*!\w*"),
        ["responsive!$ccode","modulative!$ccode"] => ByRow((i,j)->i|j) => :sr,
        ["responsive!$ccode","modulative!$ccode","enoughresponse!$ccode"] => ByRow((i,j,k)->i&j&k) => :rme,
        "f1f0!$ccode" => ByRow(i->i<1 ? "Complex" : "Simple") => :sctype,
        "frf_Ori_Final!$ccode" => ByRow(i->i.oo) => :oo, "frf_Ori_Final!$ccode" => ByRow(i->i.ocv) => :ocv,
        "frf_Ori_Final!$ccode" => ByRow(i->i.od) => :od, "frf_Ori_Final!$ccode" => ByRow(i->i.dcv) => :dcv,
        "frf_Ori_Final!$ccode" => ByRow(i->i.fit.po) => :po, "frf_Ori_Final!$ccode" => ByRow(i->i.fit.osi2) => :osi,
        "frf_Ori_Final!$ccode" => ByRow(i->sum(i.fit.dhw)) => :dhw, "frf_Ori_Final!$ccode" => ByRow(i->i.fit.gvm.r) => :dfitr,
        "frf_Ori_Final!$ccode" => ByRow(i->i.fit.pd) => :pd, "frf_Ori_Final!$ccode" => ByRow(i->i.fit.dsi2) => :dsi,
        "frf_Ori_Final!$ccode" => ByRow(i->sum(i.fit.ohw)) => :ohw, "frf_Ori_Final!$ccode" => ByRow(i->i.fit.vmn2.r) => :ofitr,
        "frf_SpatialFreq!$ccode" => ByRow(i->i.fit.psf) => :psf, "frf_SpatialFreq!$ccode" => ByRow(i->i.fit.dog.r) => :sffitr,
        "frf_SpatialFreq!$ccode" => ByRow(i->i.osf) => :osf, "frf_SpatialFreq!$ccode" => ByRow(i->i.fit.sftype) => :sftype)
end

aoscell = orisfcell(cell,ccode="A")
loscell = orisfcell(cell,ccode="L")
moscell = orisfcell(cell,ccode="M")
soscell = orisfcell(cell,ccode="S")

aoscell |> Voyager()

soscell |>
    [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"sr:n",title="OriSF"});
     @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"sr:n",title="OriSF"});
     @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"sr:n",title="OriSF"})]

soscell |> [@vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"ocv"},color={"layer",scale={scheme=:category10}});
            @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"dcv"},color={"layer",scale={scheme=:category10}});
            @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"osi"},color={"layer",scale={scheme=:category10}});
            @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"dsi"},color={"layer",scale={scheme=:category10}})]

soscell |> [@vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"ocv"},color={"sctype",scale={scheme=:set1}});
            @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"dcv"},color={"sctype",scale={scheme=:set1}});
            @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"osi"},color={"sctype",scale={scheme=:set1}});
            @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"dsi"},color={"sctype",scale={scheme=:set1}})]

soscell |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"oo",axis={values=0:45:180}},
            color={"layer",scale={scheme=:category10}},column={"site",title=""})

soscell |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"od",axis={values=0:90:360}},
            color={"layer",scale={scheme=:category10}},column={"site",title=""})

soscell |> @vlplot(:bar,y={"count()",title="Number of Cells"},x={"oo",bin={step=30},title="Optimal Ori",axis={values=0:30:180}},
            column={"site",title=""})

soscell |> @vlplot(:bar,y={"count()",title="Number of Cells"},x={"od",bin={step=30},title="Optimal Dir",axis={values=0:30:360}},
            column={"site",title=""})


aoscell |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"osf"},
            color={"layer",scale={scheme=:category10}},column={"site",title=""})

soscell |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth (μm)"},x={"osf"},
            color={"sftype",scale={scheme=:dark2}})









## Color
dklcg = cgrad(ColorMaps["lidkl_mcchue_l0"].colors)
hslcg = cgrad(ColorMaps["hsl_mshue_l0.4"].colors)

function colorcell(cell;ccode="DKL_L0")
    select!(dropmissing(cell,"responsive!$ccode"),Not(r"\w*!\w*"),
        ["responsive!$ccode","modulative!$ccode"] => ByRow((i,j)->i|j) => :sr,
        ["responsive!$ccode","modulative!$ccode","enoughresponse!$ccode"] => ByRow((i,j,k)->i&j&k) => :rme,
        ["pzfr_HueAngle!$ccode","f_HueAngle!$ccode"] => ByRow((r,f)->r[f.==0][1]-r[f.==180][1]) => :L_M,
        ["pzfr_HueAngle!$ccode","f_HueAngle!$ccode"] => ByRow((r,f)->r[f.==90][1]-r[f.==270][1]) => :S_LM,
        "frf_HueAngle!$ccode" => ByRow(i->i.oh) => :oh, "frf_HueAngle!$ccode" => ByRow(i->i.hcv) => :hcv,
        "frf_HueAngle!$ccode" => ByRow(i->i.oha) => :oha, "frf_HueAngle!$ccode" => ByRow(i->i.hacv) => :hacv,
        "frf_HueAngle!$ccode" => ByRow(i->i.fit.ph) => :ph, "frf_HueAngle!$ccode" => ByRow(i->i.fit.hsi2) => :hsi,
        "frf_HueAngle!$ccode" => ByRow(i->sum(i.fit.hhw)) => :hhw, "frf_HueAngle!$ccode" => ByRow(i->i.fit.gvm.r) => :hfitr,
        "frf_HueAngle!$ccode" => ByRow(i->i.fit.pha) => :pha, "frf_HueAngle!$ccode" => ByRow(i->i.fit.hasi2) => :hasi,
        "frf_HueAngle!$ccode" => ByRow(i->sum(i.fit.hahw)) => :hahw, "frf_HueAngle!$ccode" => ByRow(i->i.fit.vmn2.r) => :hafitr)
end

dklcell = colorcell(cell,ccode="DKL_L0")
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


## any stimuli responsive
asrcell = select!(outerjoin(dklcell[:,[:id,:sr]],hslcell[:,[:id,:sr]],aoscell[:,[:id,:sr]],loscell[:,[:id,:sr]],moscell[:,[:id,:sr]],soscell[:,[:id,:sr]],on=:id,makeunique=true),
            :id,r"sr\w*" => ByRow((i...)->any(i)) => :asr)

select!(outerjoin(cell,asrcell,on=:id),Not(r"\w*!\w*")) |>
    [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"asr:n",title="Visual"});
     @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"asr:n",title="Visual"});
     @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"asr:n",title="Visual"})]

## color responsive
ucrcell = select!(outerjoin(dklcell[:,[:id,:sr]],hslcell[:,[:id,:sr]],on=:id,makeunique=true),
            :id,r"sr\w*" => ByRow((i...)->any(i)) => :ucr)

select!(dropmissing(outerjoin(cell,ucrcell,on=:id),:ucr),Not(r"\w*!\w*")) |>
    [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"ucr:n",title="Color"});
     @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"ucr:n",title="Color"});
     @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"ucr:n",title="Color"})]

## OriSF responsive
osrcell = select!(outerjoin(aoscell[:,[:id,:sr]],loscell[:,[:id,:sr]],moscell[:,[:id,:sr]],soscell[:,[:id,:sr]],on=:id,makeunique=true),
            :id,r"sr\w*" => ByRow((i...)->any(i)) => :osr)

select!(dropmissing(outerjoin(cell,osrcell,on=:id),:osr),Not(r"\w*!\w*")) |>
    [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"osr:n",title="OriSF"});
     @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"osr:n",title="OriSF"});
     @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"osr:n",title="OriSF"})]

## STA responsive
select(cell,Not(r"\w*!\w*"),:rc=>ByRow(i->ismissing(i) ? false : true)=>:star) |>
    [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"star:n",title="STA"});
     @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"star:n",title="STA"});
     @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"star:n",title="STA"})]

select(cell,Not(r"\w*!\w*"),:rc=>ByRow(i->ismissing(i) ? false : 'S' in skipmissing(i) ? true : false)=>:star) |>
 [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"star:n",title="STA"});
  @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"star:n",title="STA"});
  @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"star:n",title="STA"})]


## Spectral STA responsive
SRTYPE = sort(join.(sort.(combinations(['A','L','M','S']))))
CTCOLORS = cgrad(:rainbow,length(CTYPES),categorical = true).colors.colors

select!(dropmissing(cell,:rc),Not(r"\w*!\w*"),:rc=>ByRow(i->join(skipmissing(i)))=>:sr) |>
    [@vlplot(:bar,y={"sr",title="Spectral Responsive"},x={"count()",title="Number of Cells"});
    @vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"sr:n",scale={scheme=:category20},title="Spectral"});
    @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"sr:n",scale={scheme=:category20},title="Spectral"});
    @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"sr:n",scale={scheme=:category20},title="Spectral"})]


crdf = DataFrame(ConeResponseType=filter(i->!isempty(i),map(t->match(r"A?([L,M,S]*)",t).captures[1],mrdf.ModelResponseType)))
p = crdf |> @vlplot(:bar,y=:ConeResponseType,x="count()",width=600,height=400)

## STA
oddeven(p) = p < 0.5 ? 4abs(p-0.25) : 4abs(p-0.75)
onoff(p) = p < 0.5 ? 1 - 4abs(p-0.25) : 4abs(p-0.75) - 1
elness(r) = r < 1 ? 1 - 1/r : r - 1
roundness(r) = abs(elness(r))
elsign(r) = sign(elness(r))
function fitpredict(fit;ppd=45,tight=false)
   if tight
       param = deepcopy(fit.param)
       if fit.model == :dog
           param[[2,4]].=0
           radius = 3*(param[6] < param[3] ? param[3] : param[6])
           fit = (;fit.model,fit.fun,param,radius)
       elseif fit.model == :gabor
           param[[2,4]].=0
           radius = 3*maximum(param[[3,5]])
           fit = (;fit.model,fit.fun,param,radius)
       end
   end
   x=y= -fit.radius:1/ppd:fit.radius
   predict(fit,x,y,yflip=true)
end

function stadogcell(cell;model=:dog,imgsize=(48,48))
    df = select!(dropmissing!(flatten(dropmissing(cell,:rc),["rc","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->imresize(fitpredict(i,tight=true),imgsize))=>"$(model)img",
        "sta!$model"=>ByRow(i->i.r)=>:r,"sta!$model"=>ByRow(i->i.bic)=>:bic,"sta!$model"=>ByRow(i->1-i.r2)=>:fvu,
        "sta!$model"=>ByRow(i->i.param[1])=>:ampe,"sta!$model"=>ByRow(i->i.param[5])=>:ampi,
        "sta!$model"=>ByRow(i->i.param[2])=>:cx,"sta!$model"=>ByRow(i->i.param[4])=>:cy,"sta!$model"=>ByRow(i->4*(i.param[6] < i.param[3] ? i.param[3] : i.param[6]))=>:diameter,
        "sta!$model"=>ByRow(i->elness(i.param[6]/i.param[3]))=>:cs,"sta!$model"=>ByRow(i->roundness(i.param[6]/i.param[3]))=>:oppo,
        "sta!$model"=>ByRow(i->elness(i.param[1]/i.param[5]))=>:onoff,"sta!$model"=>ByRow(i->roundness(i.param[1]/i.param[5]))=>:opposth,
        "sta!$model"=>ByRow(i->elsign(i.param[1]/i.param[5]))=>:sign)
end

function stagaborcell(cell;model=:gabor,imgsize=(48,48))
    df = select!(dropmissing!(flatten(dropmissing(cell,:rc),["rc","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->imresize(fitpredict(i,tight=true),imgsize))=>"$(model)img",
        "sta!$model"=>ByRow(i->i.r)=>:r,"sta!$model"=>ByRow(i->i.bic)=>:bic,"sta!$model"=>ByRow(i->1-i.r2)=>:fvu,
        "sta!$model"=>ByRow(i->i.param[1])=>:amp,"sta!$model"=>ByRow(i->i.param[2])=>:cx,"sta!$model"=>ByRow(i->i.param[4])=>:cy,
        "sta!$model"=>ByRow(i->rad2deg(i.param[6]))=>:ori,"sta!$model"=>ByRow(i->i.param[7])=>:sf,"sta!$model"=>ByRow(i->i.param[8])=>:phase,
        "sta!$model"=>ByRow(i->4*(i.param[5] < i.param[3] ? i.param[3] : i.param[5]))=>:diameter,
        "sta!$model"=>ByRow(i->elness(i.param[3]/i.param[5]))=>:el,"sta!$model"=>ByRow(i->roundness(i.param[3]/i.param[5]))=>:rd,
        "sta!$model"=>ByRow(i->4*i.param[5]*i.param[7])=>:wn,"sta!$model"=>ByRow(i->oddeven(i.param[8]))=>:oddeven,
        "sta!$model"=>ByRow(i->onoff(i.param[8]))=>:onoff,"sta!$model"=>ByRow(i->sign(onoff(i.param[8])))=>:sign)
end

dogcell = stadogcell(cell)
gaborcell = stagaborcell(cell)

gaborcell |> Voyager()

gaborcell |> [@vlplot(:bar,x={"r",bin={step=0.05,extent=[0,1]},title="r"},y={"count()",title="Number of Linear Filter"});
                @vlplot(mark={:line,size=2},x=:r,transform=[{sort=[{field=:r}],window=[{op="count",as="cum"}],frame=[nothing,0]}],
                y={"cum",title="Cumulated Number of Linear Filter"})]

doggaborgoodness = DataFrame(dogr=dogcell.r,gaborr=gaborcell.r,dogbic=dogcell.bic,gaborbic=gaborcell.bic,dogfvu=dogcell.fvu,gaborfvu=gaborcell.fvu)
doggaborgoodness |>
    [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1.5e5}]},x=:x,y=:x}, {:point,x={:gaborbic,title="BIC (Gabor)"},y={:dogbic,title="BIC (DoG)"}}]);
    @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1}]},x=:x,y=:x}, {:point,x={:gaborfvu,title="FVU (Gabor)"},y={:dogfvu,title="FVU (DoG)"}}]);
    @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1}]},x=:x,y=:x}, {:point,x={:gaborr,title="R (Gabor)"},y={:dogr,title="R (DoG)"}}])]

srfcells |> [@vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Spatial RF"},color="rc");
       @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Spatial RF"},color="rc")]

srfcells |> [@vlplot(:point,x={"ar",title="Aspect Ratio"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"});
          @vlplot(:point,x={"nc",title="Number of Cycle"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"});
          @vlplot(:point,x={"odd",title="Oddness"},y={"depth",sort="descending",title="Cortical Depth"},color={"layer"})]

gaborcell |> [@vlplot(mark={:point,size=20},x={"ar",title="Aspect Ratio"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"});
       @vlplot(mark={:point,size=20},x={"nc",title="Number of Cycle"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"});
       @vlplot(mark={:point,size=20},x={"oddeven",title="Oddness"},y={"depth",sort="descending",title="Cortical Depth"},color={"rc"})]

gaborcell |> [@vlplot(:point,x={"diameter",title="Diameter (Deg)"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"ar",title="Aspect Ratio"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"el",title="Ellipticalness"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"ori",title="Ori (Deg)"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"sf",title="SpatialFreq (Cycle/Deg)"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"nc",title="Number of Cycle"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"phase",title="Phase"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"oddeven",title="OddEven"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}});
           @vlplot(:point,x={"onoff",title="OnOff"},y={"depth",sort="descending",title="Cortical Depth (μm)"},color={"layer",scale={scheme=:category10}})]

dogcell |> [@vlplot(mark={:bar,binSpacing=0},x={"diameter",bin={step=0.03},title="Diameter (Deg)"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"ampe",bin={maxbin=30},title="Ae"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"ampi",bin={maxbin=30},title="Ai"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
           @vlplot(mark={:bar,binSpacing=0},x={"cs",bin={step=0.03},title="CenterSurround"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
           @vlplot(mark={:bar,binSpacing=0},x={"oppo",bin={step=0.03},title="Opponency"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
           @vlplot(mark={:bar,binSpacing=0},x={"onoff",bin={step=0.2},title="OnOff"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"opposth",bin={step=0.2},title="Opponency Strength"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
           @vlplot(mark={:bar,binSpacing=0},x={"sign",bin={maxbin=3},title="OnOff Sign"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}})]

gaborcell |> [@vlplot(mark={:bar,binSpacing=0},x={"diameter",bin={step=0.03},title="Diameter (Deg)"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"el",bin={step=0.1},title="Ellipticalness"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"rd",bin={step=0.1},title="Roundness"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"ori",bin={step=3.6},title="Ori (Deg)"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"sf",bin={step=0.1},title="SpatialFreq (Cycle/Deg)"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"wn",bin={step=0.1},title="Wave Number"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"phase",bin={step=0.02},title="Phase"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"oddeven",bin={step=0.02},title="OddEven"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"onoff",bin={step=0.04},title="OnOff"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}});
            @vlplot(mark={:bar,binSpacing=0},x={"sign",bin={maxbin=3},title="OnOff Sign"},y={"count()",title="Number of Linear Filter"},color={"layer",scale={scheme=:category10}})]

gaborcell |> [@vlplot(:bar,x={"ar",bin={maxbins=20},title="Aspect Ratio"},y={"count()",title="Number of Spatial RF"},color={"rc"});
       @vlplot(:bar,x={"nc",bin={maxbins=20},title="Number of Cycle"},y={"count()",title="Number of Spatial RF"},color={"rc"});
       @vlplot(:bar,x={"odd",bin={maxbins=15},title="Oddness"},y={"count()",title="Number of Spatial RF"},color={"rc"})]



## UMAP of STA
function d2clu(x;r=1)
    cs = dbscan(x,r)
    cid = zeros(Int,size(x,2))
    foreach(i->cid[cs[i].core_indices] .= i, 1:length(cs))
    cid
end
function srimg(c,i)
    mapreduce((p,q)->float32.(alphamask(q,color=spectralcm[p])),hcat,c,i)
end
dogimg = map(i->float32.(alphamask(i)),dogcell.dogimg)
gaborimg = map(i->float32.(alphamask(i)),gaborcell.gaborimg)
spectralcm = Dict('A'=>ColorMaps["dkl_mcclumiso"].colors,'L'=>ColorMaps["lms_mccliso"].colors,
                        'M'=>ColorMaps["lms_mccmiso"].colors,'S'=>ColorMaps["lms_mccsiso"].colors)
dogimgc = map((i,c)->float32.(alphamask(i,color=spectralcm[c])),dogcell.dogimg,dogcell.rc)
gaborimgc = map((i,c)->float32.(alphamask(i,color=spectralcm[c])),gaborcell.gaborimg,gaborcell.rc)

@manipulate for i in 1:length(dogimg)
    dogimg[i]
end

@manipulate for i in 1:length(gaborimg)
    gaborimg[i]
end

# DoG
features = [:cs,:oppo,:opposth,:onoff,:sign]
Y = mapfoldl(i->dogcell[:,i],hcat,features)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
    Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1000, 1000))
end

Y2 = umap(permutedims(Y), 2, n_neighbors=15, min_dist=0.12, n_epochs=300,metric=Euclidean(),init=:spectral)
Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=d2clu(Y2,r=1.0),frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
savefig(joinpath(resultroot,"sta_dog_shape_umap_clu.png"))

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_dog_shape_umap.png"),p)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimgc,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_dog_shape_umap_c.png"),p)

save("UMAP_sta.svg",plotunitpositionimage(Y2',dogimg))




dogcell.soppo = dogcell.oppo .> 0.060
Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=dogcell.soppo,frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
savefig(joinpath(resultroot,"tt.png"))

dogcell |> [@vlplot(mark={:bar,binSpacing=0},x={"oppo",bin={step=0.03},title="Opponency"},y={"count()",title="Number of Spatial RF"},color={"soppo:n"});
       @vlplot(mark={:bar,binSpacing=0},x={"opposth",bin={step=0.2},title="Opponency Strength"},y={"count()",title="Number of Spatial RF"},color={"soppo:n"})]

function soppo(s)
    us = unique(s)
    length(us) == 1 ? us[1] : -1
end
function celltype(so,co)
    so == -1 && return -1
    so == 1 ? (co ? 1 : 3) : (co ? 2 : 0)
end
dogtypecell = select!(leftjoin(combine(groupby(dogcell,:id),:soppo=>(i->soppo(i))=>:soppo,
                [:rc,:sign]=>((i,j)->join(map((p,q)->q==1 ? p : "$(p)̲",i,j)))=>:sr,:sign=>(i->[1,-1] ⊆ i)=>:coppo,
                [:rc,:dogimg]=>((c,i)->[srimg(c,i)])=>:srimg), cell,on=:id),Not(r"\w*!\w*"),[:soppo,:coppo]=>ByRow((i,j)->celltype(i,j))=>:celltype)


function plotsrsta(stacell;col=5,dir=nothing,file="")
    n = nrow(stacell)
    n==0 && return nothing
    nc = n < col ? n : col
    nr = n <= col ? 1 : ceil(Int,n/col)
    p=Plots.plot(layout=(nr,nc),frame=:none,size=(nc*450,nr*200),leg=false,titlefontsize=16)
    for i in 1:n
        Plots.plot!(p[i],stacell.srimg[i],title="$(stacell.id[i]) - $(stacell.sr[i])")
    end
    isnothing(dir) ? p : savefig(joinpath(dir,"$file.png"))
end


plotsrsta(filter(r->r.celltype == -1,dogtypecell),col=6,dir=resultroot,file="sta_dog_type-1")


















# Gabor
features = [:rd,:wn,:oddeven]
Y = mapfoldl(i->gaborcell[:,i],hcat,features)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
    Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1000, 1000))
end




Y2 = umap(permutedims(Y), 2, n_neighbors=7, min_dist=0.1, n_epochs=300,metric=Euclidean(),init=:spectral)
Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=d2clu(Y2,r=1.2),frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
savefig(joinpath(resultroot,"sta_gabor_shape_umap_clu.png"))

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap.png"),p)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimgc,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap_c.png"),p)




save("UMAP_sta.svg",plotunitpositionimage(Y2',gaborimg))





gaborcell.soppo = gaborcell.wn .> 0.5
Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=gaborcell.soppo,frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))

gaborcell |> [@vlplot(mark={:bar,binSpacing=0},x={"wn",bin={step=0.1},title="Wave Number"},y={"count()",title="Number of Spatial RF"},color={"soppo:n"});
       @vlplot(mark={:bar,binSpacing=0},x={"oddeven",bin={step=0.02},title="OddEven"},y={"count()",title="Number of Spatial RF"},color={"soppo:n"})]

gaborcell |> [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Linear Filter"},color={"soppo:n",title="Spatial"});
        @vlplot(:bar,y={"depth",bin={step=100},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Linear Filter"},color={"soppo:n",title="Spatial"});
        @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Linear Filter"},color={"soppo:n",title="Spatial"})]






gabortypecell = select!(leftjoin(combine(groupby(gaborcell,:id),:soppo=>(i->soppo(i))=>:soppo,
               [:rc,:sign]=>((i,j)->join(map((p,q)->q==1 ? p : "$(p)̲",i,j)))=>:sr,:sign=>(i->[1,-1] ⊆ i)=>:coppo,
               [:rc,:gaborimg]=>((c,i)->[srimg(c,i)])=>:srimg), cell,on=:id),Not(r"\w*!\w*"),[:soppo,:coppo]=>ByRow((i,j)->celltype(i,j))=>:celltype)

plotsrsta(filter(r->r.celltype == 3,gabortypecell),col=6,dir=resultroot,file="sta_gabor_type3")


gabortypecell |> [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"celltype:n"});
        @vlplot(:bar,y={"depth",bin={step=200},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"celltype:n"});
        @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"celltype:n"})]




















## Relation between spectral STAs
function stapair(sta,c='L',c1='M')
   transform!(innerjoin(filter(r->r.rc==c,sta),filter(r->r.rc==c1,sta),on=:id,makeunique=true),
           [:cx,:cx_1]=>ByRow((i,j)->i-j)=>:dx,[:cy,:cy_1]=>ByRow((i,j)->i-j)=>:dy)
end
function plotdogpair(sta,c='L',c1='M')
   stapair(sta,c,c1) |>
       [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1.8}]},x=:x,y=:x}, {:point,x={:oppo,title=c},y={:oppo_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="Opponency");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=9}]},x=:x,y=:x}, {:point,x={:opposth,title=c},y={:opposth_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="Opponency Strength");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1.2}]},x=:x,y=:x}, {:point,x={:diameter,title=c},y={:diameter_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="Diameter");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=-0.5,y=0},{x=0.5,y=0}]},x=:x,y=:y},{mark={:line,size=1,color="black"},data={values=[{x=0,y=-0.5},{x=0,y=0.5}]},x=:x,y=:y},
       {:point,x={:dx,title="ΔX ($c-$c1)"},y={:dy,title="ΔY ($c-$c1)"},color={"layer",scale={scheme=:category10}}}],title="Center Displacement")]
end
function plotgaborpair(sta,c='L',c1='M')
   stapair(sta,c,c1) |>
       [@vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=-6},{x=6}]},x=:x,y=:x}, {:point,x={:el,title=c},y={:el_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="Ellipticalness");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=6}]},x=:x,y=:x}, {:point,x={:rd,title=c},y={:rd_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="Roundness");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=5}]},x=:x,y=:x}, {:point,x={:wn,title=c},y={:wn_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="Wave Number");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=1}]},x=:x,y=:x}, {:point,x={:oddeven,title=c},y={:oddeven_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="OddEven");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=-1,y=1},{x=1,y=-1}]},x=:x,y=:y}, {:point,x={:onoff,title=c},y={:onoff_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="OnOff");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=2}]},x=:x,y=:x}, {:point,x={:diameter,title=c},y={:diameter_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="Diameter");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0},{x=180}]},x=:x,y=:x}, {:point,x={:ori,title=c},y={:ori_1,title=c1},color={"layer",scale={scheme=:category10}}}],title="Orientation");
       @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=-0.5,y=0},{x=0.5,y=0}]},x=:x,y=:y},{mark={:line,size=1,color="black"},data={values=[{x=0,y=-0.5},{x=0,y=0.5}]},x=:x,y=:y},
       {:point,x={:dx,title="ΔX ($c-$c1)"},y={:dy,title="ΔY ($c-$c1)"},color={"layer",scale={scheme=:category10}}}],title="Center Displacement")]
end

plotdogpair(dogcell,'A','L')
plotgaborpair(gaborcell,'M','S')


## Cone Weights
function rcw(c,a,s;cone='L')
    lmsi = 1:length(a)
    ai = findfirst(c.=='A')
    isnothing(ai) || (lmsi = setdiff(lmsi,ai))
    ci = findfirst(c.==cone)
    isnothing(ci) ? 0 : s[ci]*a[ci]/sum(a[lmsi])
end

cwcell = select!(leftjoin(combine(groupby(filter(r->r.rc != 'A',gaborcell),:id),[:rc,:amp,:sign]=>((c,a,s)->rcw(c,a,s,cone='L'))=>:wl,
                [:rc,:amp,:sign]=>((c,a,s)->rcw(c,a,s,cone='M'))=>:wm, [:rc,:amp,:sign]=>((c,a,s)->rcw(c,a,s,cone='S'))=>:ws),
                cell,on=:id),Not(r"\w*!\w*"))

cwcell |> @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0,y=1},{x=1,y=0}]},x=:x,y=:y},
                            {mark={:line,size=1,color="black"},data={values=[{x=-1,y=0},{x=0,y=1}]},x=:x,y=:y},
                            {mark={:line,size=1,color="black"},data={values=[{x=-1,y=0},{x=0,y=-1}]},x=:x,y=:y},
                            {mark={:line,size=1,color="black"},data={values=[{x=0,y=-1},{x=1,y=0}]},x=:x,y=:y},
                            {:point,x={:wl,title="L Cone Weight"},y={:wm,title="M Cone Weight"},color={"layer",scale={scheme=:category10}}}])


## Spectral Weights
function rsw(c,a,s,i;ch='L')
    ci = findfirst(c.==ch)
    isnothing(ci) ? 0 : a[ci]/sum(a)
end

function srimg(c,i)
    mapreduce((p,q)->ismissing(p) ? [] : float32.(alphamask(q,color=spectralcm[p])),hcat,c,i)
end

rswcell = select!(leftjoin(combine(groupby(gaborcell,:id),[:rc,:amp,:sign,:sclu]=>((c,a,s,i)->rsw(c,a,s,i,ch='A'))=>:wa,:rc=>(i->join(i))=>:sr,
    [:rc,:amp,:sign,:sclu]=>((c,a,s,i)->rsw(c,a,s,i,ch='L'))=>:wl,[:rc,:amp,:sign,:sclu]=>((c,a,s,i)->rsw(c,a,s,i,ch='M'))=>:wm,
    [:rc,:amp,:sign,:sclu]=>((c,a,s,i)->rsw(c,a,s,i,ch='S'))=>:ws,[:rc,:gaborimg]=>((c,i)->[srimg(c,i)])=>:srimg), cell,on=:id),Not(r"\w*!\w*"))

features = [:wa,:wl,:wm,:ws]
Y = mapfoldl(i->rswcell[:,i],hcat,features)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
    Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:auto,stroke(0)),group=rswcell.sr,frame=:none,leg=:inline)
end

Y2 = umap(permutedims(Y), 2, n_neighbors=10, min_dist=2, n_epochs=300,metric=Euclidean(),init=:spectral)
Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=rswcell.sr,frame=:none,leg=:inline)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=rswcell.srimg,markersize=size.(rswcell.srimg),scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap_c.png"),p)





save("UMAP_c_sta.svg",plotunitpositionimage(Y2',rswcell.srimg))

save("UMAP_c_sta_layer.svg",plotunitlayerimage(rswcell.layer,rswcell.imgs,width=2000,markersize=30))



function plotsrsta(stacell;layer="4Cb",mc=4,dir=nothing)
    vc = filter(r->r.layer==layer,stacell)
    n = nrow(vc)
    n==0 && return nothing
    nc = n < mc ? n : mc
    nr = n <= mc ? 1 : ceil(Int,n/mc)
    p=Plots.plot(layout=(nr,nc),frame=:none,size=(nc*450,nr*200),leg=false,titlefontsize=14)
    for i in 1:n
        Plots.plot!(p[i],vc.srimg[i],title=vc.id[i])
    end
    isnothing(dir) ? p : savefig(joinpath(dir,"sta_$layer.png"))
end

for l in levels(rswcell.layer)
    plotsrsta(rswcell,layer=l,dir=resultroot)
end



## Cell Graph
cellgraph = load(joinpath(resultroot,"cellgraph.jld2"),"cellgraph")
foreach(i->set_prop!(cellgraph,i,:layer,cell.layer[findfirst(cell.id.==get_prop(cellgraph,i,:id))]),1:nv(cellgraph))


layercolor=Dict("Out"=>RGB(0.1,0.1,0.1),
                "1"=>RGB(0.2,0.9,0.9),
                "23"=>RGB(0.12,0.46,0.7),
                "2"=>[2550,1800],
                "3"=>[2550,1800],
                "4AB"=>RGB(1,0.5,0.05),
                "4A"=>[2575,1800],
                "4B"=>[2430,1800],
                "4C"=>[1800,1800],
                "4Ca"=>RGB(0.17,0.62,0.17),
                "4Cb"=>RGB(0.83,0.15,0.15),
                "56"=>RGB(0.58,0.4,0.74),
                "5"=>[1300,1800],
                "6"=>[1500,1800],
                "WM"=>RGB(0.9,0.9,0.9))

graphplot(cellgraph,nodeshape=:circle,edgecolor=:gray10,linewidth=0.3,curvature_scalar=0.01,arrow=0.13,method=:stress,
        markerstrokecolor=:gray10,markerstrokewidth=0.3,nodecolor=map(i->layercolor[get_prop(cellgraph,i,:layer)],1:nv(cellgraph)))
foreach(ext->savefig(joinpath(resultroot,"graph_layer$ext")),figfmt)

graphsummary = transform(transform(transform(combine(groupby(cell,:site),nrow=>:n),
    :n=>ByRow(i->binomial(i,2))=>:np,:site=>ByRow(i->collect(filter_vertices(cellgraph,:site,i)))=>:vs),
    :vs=>ByRow(length)=>:nv,:vs=>ByRow(i->length(collect(edges(induced_subgraph(cellgraph,i)[1]))))=>:ne),
    [:ne,:np]=>ByRow((i,j)->i/j*100)=>:pp)

graphsummary |> @vlplot(:bar,y={"site",title="Recording Site"},x={"pp",title="Projections / Pairs (%)"})



function layerprojection(cellgraph;src="4Cb",dst="23")
    df = DataFrame(map(e->(src=get_prop(cellgraph,e.src,:id),dst=get_prop(cellgraph,e.dst,:id)),
    filter_edges(cellgraph,(g,e)->get_prop(g,e.src,:layer)==src && get_prop(g,e.dst,:layer)==dst)))
    df[!,:projtype] .= "$(src) ⟶ $(dst)";df
end
function graphlayerprojection(cellgraph)
    vls=unique(map(i->get_prop(cellgraph,i,:layer),vertices(cellgraph)))
    mapreduce(i->layerprojection(cellgraph,src=i[1],dst=i[2]),append!,permutations(vls,2),
    init=mapreduce(i->layerprojection(cellgraph,src=i,dst=i),append!,vls))
end

projsummary = graphlayerprojection(cellgraph)

projsummary |> @vlplot(:bar,y={"projtype",title="Layer Projection"},x={"count()",title="Number of Projections"})



function plotprojsta(ps,stacell;src=nothing,dst=nothing,ptype="4Cb ⟶ 23",dir=nothing)
    mi(x) = isnothing(x) ? missing : x
    isnothing(src) || isnothing(dst) || (ptype = "$(src) ⟶ $(dst)")
    vpj=dropmissing!(transform!(filter(r->r.projtype == ptype,ps),:src=>ByRow(i->mi(findfirst(j->j==i,stacell.id)))=>:srci,
            :dst=>ByRow(i->mi(findfirst(j->j==i,stacell.id)))=>:dsti),[:srci,:dsti])
    n = nrow(vpj)
    n==0 && return nothing
    p=Plots.plot(layout=(n,2),frame=:none,size=(2*450,n*200),leg=false,titlefontsize=10)
    for i in 1:n
        Plots.plot!(p[i,1],stacell.srimg[vpj.srci[i]],title=vpj.src[i])
        Plots.plot!(p[i,2],stacell.srimg[vpj.dsti[i]],title=vpj.dst[i])
    end
    isnothing(dir) ? p : savefig(joinpath(dir,"$ptype.png"))
end

pstadir = joinpath(resultroot,"psta")
isdir(pstadir) || mkpath(pstadir)
for pt in levels(projsummary.projtype)
    plotprojsta(projsummary,rswcell,ptype=pt,dir=pstadir)
end



function vergenttree(g,v;dir=:in)
    vs = dir==:in ? inneighbors(g,v) : outneighbors(g,v)
    length(vs) < 2 && return missing
    (vs=get_prop.([g],vs,:id), v=get_prop(g,v,:id),dir)
end

vergentsummary = [collect(skipmissing(vergenttree(cellgraph,i,dir=:in) for i in 1:nv(cellgraph)));
                collect(skipmissing(vergenttree(cellgraph,i,dir=:out) for i in 1:nv(cellgraph)))]

function plottreesta(tree,stacell;dir=nothing)
    mi(x) = isnothing(x) ? missing : x
    vi=findfirst(i->i==tree.v,stacell.id)
    isnothing(vi) && return nothing
    vsi = collect(skipmissing(mi(findfirst(i->i==v,stacell.id)) for v in tree.vs))
    n = length(vsi)
    n<2 && return nothing
    p=Plots.plot(layout=(n+1,1),frame=:none,size=(400,(n+1)*200),leg=false,titlefontsize=10)
    title = "$(tree.v) - $(stacell.layer[vi]) - $(tree.dir)"
    Plots.plot!(p[1,1],stacell.srimg[vi],title=title,titlefontsize=15)
    for i in 1:n
        Plots.plot!(p[i+1,1],stacell.srimg[vsi[i]],title="$(stacell.id[vsi[i]]) - $(stacell.layer[vsi[i]])")
    end
    isnothing(dir) ? p : savefig(joinpath(dir,"$title.png"))
end

iostadir = joinpath(resultroot,"iosta")
isdir(iostadir) || mkpath(iostadir)
for t in vergentsummary
    plottreesta(t,rswcell,dir=iostadir)
end






save("UMAP_c_sta_layer.svg",plotunitlayerimage(rswcell.layer,rswcell.imgs,width=2000,markersize=30))


vi = filter(i->"4Cb"==get_prop(cellgraph,i,:layer),1:nv(cellgraph))

outneighbors(cellgraph,9)

outneighbors(cellgraph,14)



t=bfs_tree(cellgraph,9)
t1=dfs_tree(cellgraph,9)
t=union(t,t1)

t=diffusion(cellgraph,1,100,initial_infections=[9])

tt = reduce(union,t)

sg,vmap = induced_subgraph(cellgraph,tt)
sg,vmap = induced_subgraph(cellgraph,collect(edges(t)))
graphplot(sg,nodeshape=:circle,nodecolor=:lightgray,markersize=3,method=:spectrual,names=vmap,dim=3,linewidth=3)






g = [0 1 1;
     0 0 1;
     0 1 0]

graphplot(g, names=1:3,edgecolor=:red,linewidth=3, curvature_scalar=0.1,arrow=0.4,dim=2,markersize=0.3)
