using NeuroAnalysis,DataFrames,JLD2,StatsBase,StatsPlots,VegaLite,XLSX,Combinatorics,DataStructures,CircStats,HypothesisTests

function collectorisf(indir;unit=DataFrame(),datafile="orisfdataset.jld2",cch=["A","L","M","S"])
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile),"dataset")
            siteid = dataset["siteid"]
            id = ["$(siteid)_$(g ? "S" : "M")U$u" for (u,g) in dataset["ugood"]]
            cs = intersect(dataset["ccode"],cch)
            ci = indexin(cs,dataset["ccode"])

            df = DataFrame(:siteid=>siteid,:id=>id,
                    :ccode => Tuple(cs),
                    :cresponsive=>[v[ci] for v in values(dataset["responsive"])],
                    :cmodulative=>[v[ci] for v in values(dataset["modulative"])],
                    :cenoughresponse=>[v[ci] for v in values(dataset["enoughresponse"])],
                    :f10 => [v[ci] for v in values(dataset["f10"])],
                    :f1phase => [v[ci] for v in values(dataset["f1phase"])],
                    :f1mag => [v[ci] for v in values(dataset["f1mag"])],
                    (f=>[v[ci] for v in values(dataset["frf"][f])] for f in keys(dataset["frf"]))...)
            append!(unit,df,cols=:union)
        end
    end
    return unit
end

resultroot = "Z:/"
orisfunit = collectorisf(resultroot)
jldsave(joinpath(resultroot,"orisfunit.jld2");orisfunit)


## cortical unit OriSF responsive
orisfunit = load(joinpath(resultroot,"orisfunit.jld2"),"orisfunit")
layer,nlbt = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate","nlbt")
layercg = cgrad(:tab10,nlbt,categorical=true)
allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
allcsu = subset(allcu,:good)
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
select!(penetration,[:siteid,:od,:cofd,:ori_lhi,:pid],:cofd=>ByRow(i->ismissing(i) ? i : first.(split(i,", ")) |> join ∘ sort)=>:cofd_type)

orisfcu = innerjoin(allcu,orisfunit,on=[:siteid,:id])
excludesites = ["AG1_V1_ODR8","AG1_V1_ODL17","AG2_V1_ODL18"]
filter!(r->r.siteid ∉ excludesites,orisfcu)
leftjoin!(orisfcu,penetration,on=:siteid)
transform!(orisfcu,:cofd_type=>ByRow(i->any(occursin.(['L','M','S'],i)))=>:cofd_C,
        :cofd_type=>ByRow(i->any(occursin.(['L','M'],i)))=>:cofd_LM,:cofd_type=>ByRow(i->contains(i,'S'))=>:cofd_S,
        :cofd_type=>ByRow(i->contains(i,'A'))=>:cofd_A,:cofd_type=>ByRow(i->i=="N")=>:cofd_N)
transform!(orisfcu,[:cresponsive,:cenoughresponse,:ccode]=>ByRow((r,e,c)->all(.!(r.&e)) ? missing : join(c[findall(r.&e)]))=>:RTYPE)
orisfcu.responsive = map(i->ismissing(i) ? false : true,orisfcu.RTYPE)
orisfcu.CTYPE = map(i->ismissing(i) ? missing : replace(i,r"A([L,M,S]+)"=>s"\1"),orisfcu.RTYPE)
orisfcsu = subset(orisfcu,:good)
orisfrcsu = dropmissing(orisfcsu,:RTYPE)
orisfnrcsu = subset(orisfcsu,:responsive=>.!)
jldsave(joinpath(resultroot,"orisfcsu.jld2");orisfcsu)
orisfdir = joinpath(resultroot,"OriSF");mkpath(orisfdir)
figfmt = [".svg",".png"]


import CairoMakie as mk
function plotdepthhist(unit;g=:responsive,dir=nothing,figfmt=[".png"],layer=nothing,xlabel="Number of Units",
                siteid="All",ylabel="Normalized Cortical Depth",size=(450,750),palette=:tab20,binedges=0:0.02:1)
    sunit = siteid=="All" ? unit : filter(r->r.siteid==siteid,unit)
    title="$(siteid)_DepthHist_$g"
    # try
    group = sunit[!,g]
    ug = unique(group) |> sort!
    cc = cgrad(palette;categorical=true)
    gh = [fit(Histogram,sunit.aligndepth[group.==g],binedges).weights for g in ug]
    sgh = accumulate(.+,gh)
    pgh = pushfirst!(sgh[1:end-1],zeros(Int,length(binedges)-1))
    bw = binedges[2]-binedges[1]
    y=range(bw/2,step=bw,length=length(binedges)-1)
    xmax = maximum(sgh[end])
    xlims = (-0.08xmax,1.01xmax)
    ltx = xlims[1]+0.01(xlims[2]-xlims[1])

    f = mk.Figure(;size,figure_padding=10)
    ax = mk.Axis(f[1,1];title,titlesize=14,xlabel,ylabel,yreversed=true,yticks=[0,1],limits=(xlims,(-0.02,1.02)),xgridvisible=false,ygridvisible=false)
    if !isnothing(layer)
        mk.hlines!(ax,[l[1] for l in values(layer)];color=:gray25,linewidth=0.7)
        mk.text!(ax,[(ltx,mean(layer[k])) for k in keys(layer)];text=collect(keys(layer)),align=(:left,:center),color=:gray10,fontsize=10)
    end
    for i in eachindex(ug)
        mk.barplot!(ax,y,sgh[i],direction=:x,fillto=pgh[i],color=cc[i],label=string(ug[i]),gap=0)
    end
    mk.axislegend(ax,labelsize=10,patchsize=(10,10))

    isnothing(dir) ? f : foreach(ext->save(joinpath(dir,"$title$ext"),f),figfmt)
    # catch
    # end
end

function plotdepthscatter(unit;x=:Diameter,xfun=first,g=:rc,dir=nothing,siteid="All",figfmt=[".png"],layer=nothing,xlabel="$x",ts="",
            ylabel="Normalized Cortical Depth",size=(450,750),palette=:tab20,ma=-0.7,
            gwfun=x->isempty(x) ? 0 : median(x),gwbfun=(x->isempty(x) ? 0 : percentile(x,25),x->isempty(x) ? 0 : percentile(x,75)))
    sunit = siteid=="All" ? unit : filter(r->r.siteid==siteid,unit)
    ts = isempty(ts) ? ts : "_"*ts
    title="$(siteid)_Depth-$(x)_$g$ts"
    xx = xfun.(sunit[!,x])
    yy = sunit.aligndepth
    group = sunit[!,g]
    ixi = isnan.(xx) .| isinf.(xx)
    if any(ixi)
        xx=xx[.!ixi]
        yy = yy[.!ixi]
        group = group[.!ixi]
    end
    xmin,xmax = extrema(xx)
    xr = xmax-xmin
    xlims = (xmin-0.08xr,xmax+0.01xr)
    ltx = xlims[1]+0.01(xlims[2]-xlims[1])
    f = mk.Figure(;size,figure_padding=10)
    ax = mk.Axis(f[1,1];title,titlesize=14,xlabel,ylabel,yreversed=true,yticks=[0,1],limits=(xlims,(-0.02,1.02)),xgridvisible=false,ygridvisible=false)
    if !isnothing(layer)
        mk.hlines!(ax,[l[1] for l in values(layer)];color=:gray25,linewidth=0.7)
        mk.text!(ax,[(ltx,mean(layer[k])) for k in keys(layer)];text=collect(keys(layer)),align=(:left,:center),color=:gray10,fontsize=10)
    end
    ug = unique(group) |> sort!
    cc = cgrad(palette;categorical=true)
    mca,msca = ma > 0 ? (ma,0) : (0,-ma)
    for (idx,g) in enumerate(ug)
        gi = findall(i->i==g,group)
        mk.scatter!(ax,xx[gi],yy[gi];markersize=5,strokewidth=1,color=(cc[idx],mca),strokecolor=(cc[idx],msca),label=string(g))
        if !isnothing(gwfun)
            n,y = unitdensity(yy[gi];w=xx[gi],wfun=gwfun,spacerange=(0,1),bw=0.02,step=0.01,s=1)
            mk.lines!(ax,n,y,color=cc[idx],linewidth=2,label="50%")
        end
        if !isnothing(gwbfun)
            n,y = unitdensity(yy[gi];w=xx[gi],wfun=gwbfun[1],spacerange=(0,1),bw=0.02,step=0.01,s=1)
            mk.lines!(ax,n,y,color=(cc[idx],0.8),linewidth=1.5,label="25%")
            n,y = unitdensity(yy[gi];w=xx[gi],wfun=gwbfun[2],spacerange=(0,1),bw=0.02,step=0.01,s=1)
            mk.lines!(ax,n,y,color=(cc[idx],0.8),linewidth=1.5,label="75%")
        end
    end

    # n,y = unitdensity(sunit.aligndepth;w=xx,wfun,spacerange=(0,1),bw=0.02,step=0.01)
    # mk.lines!(ax,n,y,color=(:gray40,0.9),linewidth=1.5,label="Average")

    mk.axislegend(ax,labelsize=10)
    
    isnothing(dir) ? f : foreach(ext->save(joinpath(dir,"$title$ext"),f),figfmt)
end


plotdepthhist(allcsu;g=:spiketype,layer)
# plotdepthhist(allcsu;g=:spiketype,layer,dir=resultroot,figfmt)

plotdepthhist(orisfcsu;g=:responsive,layer)
plotdepthhist(orisfrcsu;g=:RTYPE,layer)
plotdepthhist(orisfrcsu;g=:CTYPE,layer)

# plotdepthhist(orisfcsu;g=:responsive,layer,dir=orisfdir,figfmt)
# plotdepthhist(orisfrcsu;g=:RTYPE,layer,dir=orisfdir,figfmt)
# plotdepthhist(orisfrcsu;g=:CTYPE,layer,dir=orisfdir,figfmt)

# foreach(siteid->plotdepthhist(orisfcsu;dir=orisfdir,figfmt,layer,g=:responsive,siteid),levels(orisfcsu.siteid))
# foreach(siteid->plotdepthhist(orisfrcsu;dir=orisfdir,figfmt,layer,g=:RTYPE,siteid),levels(orisfrcsu.siteid))
# foreach(siteid->plotdepthhist(orisfrcsu;dir=orisfdir,figfmt,layer,g=:CTYPE,siteid),levels(orisfrcsu.siteid))

# pl = select(orisfrcsu,[:RTYPE,:siteid,:layer]) |>
#     [@vlplot(:bar,y={"layer"},x={"count()",axis={title="Number of Units",grid=false}},color={"RTYPE:n",scale={scheme=:category20}});
#     @vlplot(:bar,y={"siteid"},x={"count()",axis={title="Number of Units",grid=false}},color={"RTYPE:n",scale={scheme=:category20}})]
# foreach(ext->save(joinpath(orisfdir,"Layer_SiteID_Hist_RTYPE$ext"),pl),figfmt)


pl = select(orisfrcsu,[:siteid,:cofd_type]) |> unique |> @vlplot(:bar,y=:cofd_type,x={"count()",axis={title="Number of Penetration"}})
foreach(ext->save(joinpath(orisfdir,"Penetration_COFD_TYPE$ext"),pl),figfmt)
pl = select(orisfrcsu,[:siteid,:cofd_type]) |> @vlplot(:bar,y=:cofd_type,x={"count()",axis={title="Number of Responsive Unit"}})
foreach(ext->save(joinpath(orisfdir,"rUnit_COFD_TYPE$ext"),pl),figfmt)

pl = select(orisfrcsu,[:siteid,:ori_lhi]) |> unique |> @vlplot(:bar,x={"ori_lhi",bin=true},y={"count()",axis={title="Number of Penetration"}})
foreach(ext->save(joinpath(orisfdir,"Penetration_ORI_LHI$ext"),pl),figfmt)
pl = select(orisfrcsu,[:siteid,:ori_lhi]) |> @vlplot(:bar,x={:ori_lhi,bin=true},y={"count()",axis={title="Number of Responsive Unit"}})
foreach(ext->save(joinpath(orisfdir,"rUnit_ORI_LHI$ext"),pl),figfmt)


t |> @vlplot(mark={:boxplot, extent="min-max"},x=:osi,y=:spiketype)

@df t boxplot(:spiketype,:osi,color=:gray)

## Spectral Mixing of PSTH F1
f1mix = (RTYPE,ccode,f1o0,f1mag,f1phase) -> begin
    n = length(RTYPE)
    n == 1 && return missing
    mp=OrderedDict()
    for (i,j) in combinations(RTYPE,2)
        c = i*j
        ci = [findfirst(c->c[1]==ch,ccode) for ch in (i,j)]
        f1m = f1mag[ci]
        f1p = f1phase[ci]
        f10 = f1o0[ci]
        oi = argmax(reduce(.+,f1m))
        of1phase = map(p->p[oi],f1p)
        of10 = map(p->p[oi],f10)
        mof10 = mean(of10)
        pd = mod(circ_dist(of1phase...),π)
        mp[c] = (;(Symbol(ch,"phase")=>p for (ch,p) in zip((i,j),of1phase))...,:mf10=>isnan(mof10) ? missing : clamp(mof10,0,2),
            (Symbol(ch,"f10")=>p for (ch,p) in zip((i,j),of10))...,Symbol(c,"Δphase")=>pd)
    end
    n == 2 && return mp
    for (i,j,k) in combinations(RTYPE,3)
        c = i*j*k
        ci = [findfirst(c->c[1]==ch,ccode) for ch in (i,j,k)]
        f1m = f1mag[ci]
        f1p = f1phase[ci]
        f10 = f1o0[ci]
        oi = argmax(reduce(.+,f1m))
        of1phase = map(p->p[oi],f1p)
        of10 = map(p->p[oi],f10)
        mof10 = mean(of10)
        pd = mod.(circ_dist(of1phase[1],of1phase[2:3]),π)
        mp[c] = (;(Symbol(ch,"phase")=>p for (ch,p) in zip((i,j,k),of1phase))...,:mf10=>isnan(mof10) ? missing : clamp(mof10,0,2),
            (Symbol(ch,"f10")=>p for (ch,p) in zip((i,j,k),of10))...,(Symbol(ch,"Δphase")=>p for (ch,p) in zip((i*j,i*k),pd))...)
    end
    mp
end
mixparam(unit,key) = transform!(subset!(dropmissing(unit,:f1mix),:f1mix=>ByRow(d->haskey(d,key))),:f1mix=>ByRow(d->d[key])=>AsTable)
transform!(orisfrcsu,[:RTYPE,:ccode,:f10,:f1mag,:f1phase]=>ByRow(f1mix)=>:f1mix)

lm = mixparam(orisfrcsu,"LM")
@df lm corrplot(cols([:Lphase,:Mphase]),bin=50,leg=false,size=(850,650),grid=false)

# foreach(k->begin
#     cs = [Symbol(c,"phase") for c in k]
#     @df mixparam(orisfrcsu,k) corrplot(cols(cs),bin=50,leg=false,size=(850,650),grid=false)
#     foreach(ext->savefig(joinpath(orisfdir,"orisfrcsu_$(k)_mix$ext")),figfmt)
# end,["LM","LS","MS","AL","AM","AS"])

plotdepthscatter(lm;x=:LMΔphase,g=:RTYPE,layer)
plotdepthscatter(lm;x=:mf10,g=:RTYPE,layer)
plotdepthscatter(lm;x=:L,g=:RTYPE,layer)

t = filter(r->r.mf10 > 1,dropmissing(lm,:mf10))


plotdepthscatter(t;x=:LMΔphase,g=:RTYPE,layer)

pyplot()



circ_dist(0,1.5π)





## OriSF tuning properties across layers
acolor = RGB(0.3,0.3,0.3)
lcolor = ColorMaps["lms_mccliso"].colors[end]
mcolor = ColorMaps["lms_mccmiso"].colors[end]
scolor = ColorMaps["lms_mccsiso"].colors[end]
cc = Dict('A'=>acolor,'L'=>lcolor,'M'=>mcolor,'S'=>scolor)
aoricg = cgrad(fill(acolor,2))
loricg = cgrad(fill(lcolor,2))
moricg = cgrad(fill(mcolor,2))
soricg = cgrad(fill(scolor,2))

function getfrf(unit;c='A',cch="ALMS")
    ci = findfirst(c,cch)
    frf = transform(unit,:Ori_Final=>ByRow(i->ismissing(i[ci]) ? missing : i[ci])=>identity,:SpatialFreq=>ByRow(i->ismissing(i[ci]) ? missing : i[ci])=>identity) |> dropmissing!
    select!(frf,Not([:Ori_Final,:SpatialFreq]),:Ori_Final=>ByRow(i->
        (;oup=i.oup,ou=i.oup>=0.05,dup=i.dup,du=i.dup>=0.05,ocv=i.ocv,dcv=i.dcv,
        osi=i.fit.osi2[2],dsi=i.fit.osi2[1],hw=mean(i.fit.ohw),pd=mod(i.fit.po+90,360),po=mod(i.fit.po,180)))=>AsTable,
        :SpatialFreq=>ByRow(i->(;sfup=i.up,sfu=i.up>=0.05,psf=i.fit.psf,
        sfhw=i.fit.sfhw,sfpt=i.fit.sfpt,sfbw=i.fit.sfbw))=>AsTable)
    transform(groupby(frf,:siteid),[:osi,:dsi,:ocv,:dcv,:hw].=>zscore.=>[:zosi,:zdsi,:zocv,:zdcv,:zhw])
end

cfrf = Dict(c=>getfrf(orisfrcsu;c) for c in ['A','L','M','S'])

plotdepthhist(cfrf['A'];g=:ou,layer)
plotdepthscatter(cfrf['A'];layer,x=:ocv,g=:responsive)
plotdepthscatter(cfrf['A'];layer,x=:hw,g=:responsive)
plotdepthscatter(cfrf['A'];layer,x=:po,g=:responsive)

plotdepthhist(cfrf['A'];g=:sfpt,layer)
plotdepthscatter(cfrf['A'];layer,x=:sfbw,g=:responsive)
plotdepthscatter(cfrf['A'];layer,x=:psf,g=:responsive)

for c in keys(cfrf)
    # plotdepthscatter(cfrf[c];layer,x=:ocv,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:hw,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:osi,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:zocv,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:zhw,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:zosi,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])

    plotdepthscatter(cfrf[c];layer,x=:osi,g=:cofd_N,dir=orisfdir,figfmt,ts=c,palette=[cc[c],RGB(0,0,1)])
    plotdepthscatter(cfrf[c];layer,x=:zosi,g=:cofd_N,dir=orisfdir,figfmt,ts=c,palette=[cc[c],RGB(0,0,1)])

    # plotdepthscatter(cfrf[c];layer,x=:dcv,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:dsi,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:zdcv,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:zdsi,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])

    plotdepthscatter(cfrf[c];layer,x=:dsi,g=:cofd_N,dir=orisfdir,figfmt,ts=c,palette=[cc[c],RGB(0,0,1)])
    plotdepthscatter(cfrf[c];layer,x=:zdsi,g=:cofd_N,dir=orisfdir,figfmt,ts=c,palette=[cc[c],RGB(0,0,1)])

    # plotdepthscatter(cfrf[c];layer,x=:po,g=:responsive,dir=orisfdir,figfmt,ts=c,palette=[cc[c]])
end


t=cfrf['A']

@df t scatter(:ori_lhi,:osi,group=:layer)

@df filter(r->0.0<= r.aligndepth <0.17,t) corrplot(cols(indexin([:ori_lhi,:ocv,:osi],propertynames(t))))
@df filter(r->0.0<= r.aligndepth <0.17,t) corrplot([:hw,:ocv,:osi])
@df filter(r->r.layer=="6B",t) scatter(:ori_lhi,:osi)
@df filter(r->r.layer=="2/3A" && r.spiketype in ["E","En"],t) scatter(:ori_lhi,:zocv)





using PairPlots



df = select(filter(r->0.0<= r.aligndepth <0.3,t),[:ori_lhi,:dcv,:dsi,:hw])
df = select(filter(r->r.layer=="6B",t),[:ori_lhi,:ocv,:osi,:hw])
df.ori_lhi=float.(df.ori_lhi)
pairplot(df=>(PairPlots.Scatter(markersize=5),PairPlots.TrendLine(),PairPlots.Correlation()))






pl = select(orisfrcsu,[:siteid,:ori_lhi,:cofd_type,:cofd_N,:cofd_A,:cofd_C]) |> unique |> @vlplot(:bar,x={"ori_lhi",bin=true},y={"count()",axis={title="Number of Penetration"}},color=:cofd_A)
foreach(ext->save(joinpath(orisfdir,"Penetration_ORI_LHI$ext"),pl),figfmt)








pl = t |> @vlplot(mark={:boxplot, extent="min-max"},y=:zdsi,column={"layer"},x={"cofd_N:o"},color={"cofd_N:o"})
l="5"
ht = MannWhitneyUTest(filter(r->r.layer==l && r.cofd_N,t).zdsi,filter(r->r.layer==l && !r.cofd_N,t).zdsi)

@df t boxplot(:layer,:osi,group=:cofd_N)

ai = indexin([:ori_lhi,:ocv,:osi],propertynames(t))


tt=transform(groupby(t,:siteid),:po=>(i->i.-rad2deg(circ_mean(deg2rad.(i))[1]))=>:zpo)
tt=transform(groupby(t,:siteid),:po=>(i->mod.(i.-rad2deg(circ_mean(deg2rad.(2i))[1]/2),180))=>:zpo)


plotdepthscatter(t;layer,x=:zosi,g=:spiketype,palette=[cc['A'],RGB(0,0,1)])

tt=transform!(filter(r->r.spiketype in ["E","En","I"],t),:spiketype=>ByRow(i->i=="I" ? "I" : "E")=>:EI)
plotdepthscatter(tt;layer,x=:zdsi,g=:EI,palette=:tab10)

filter!(r->!r.cofd_N,t)
transform!(t,[:cofd_A,:cofd_C]=>ByRow((a,c)->a&!c)=>:cofd_g)
transform!(t,[:cofd_none,:cofd_A,:cofd_c,:cofd_S,:cofd_LM]=>ByRow((n,a,c,s,lm)->n ? "None" : (a & !c) ? "A" : (s & !lm) ? "S" : "LM")=>:cofd_g)





function orisfunit(unit;ccode='A')
    df = dropmissing(unit,"responsive!$ccode")
    df.c .= ccode
    select!(df,Not(r"\w*!\w*"),
        ["responsive!$ccode","enoughresponse!$ccode"] => ByRow(&) => :er,
        "f1f0!$ccode" => :f1f0,
        "fa_Ori_Final!$ccode" => :ori,
        "fa_SpatialFreq!$ccode" => :sf,
        ["fms!$ccode","maxfri_Ori_Final!$ccode"]=>ByRow((fm,fi)->fm[fi...])=>:orir,
        ["fms!$ccode","maxfri_SpatialFreq!$ccode"]=>ByRow((fm,fi)->fm[fi...])=>:sfr,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : (;i.fit.mfit.model,i.fit.mfit.fun,i.fit.mfit.param)) => :ofit,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : i.oup) => :oup,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : i.ocm) => :ocm,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : i.ocv) => :ocv,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : i.dup) => :dup,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : i.dcm) => :dcm,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : i.dcv) => :dcv,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : i.fit.po) => :po,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? [missing,missing] : i.fit.osi2) => [:dsi,:osi],
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : sum(i.fit.ohw)) => :ohw,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : i.fit.mfit.r) => :ocor,
        "frf_Ori_Final!$ccode" => ByRow(i->ismissing(i) ? missing : 1-i.fit.mfit.r2) => :ofvu,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : (;i.fit.mfit.model,i.fit.mfit.fun,i.fit.mfit.param)) => :sffit,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : i.up) => :sfup,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : i.msf) => :msf,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : i.fit.psf) => :psf,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : i.fit.sftype) => :sftype,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : i.fit.sfbw) => :sfbw,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : i.fit.sfpw) => :sfpw,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : i.fit.mfit.r) => :sfcor,
        "frf_SpatialFreq!$ccode" => ByRow(i->ismissing(i) ? missing : 1-i.fit.mfit.r2) => :sffvu)
end

allorisfunit = mapreduce(c->orisfunit(orisftest,ccode=c),(i,j)->append!(i,j,cols=:union),['A','L','M','S'])
transform!(allorisfunit,:f1f0=>ByRow(i->i<=1 ? "Complex" : "Simple")=>:sctype,
            :oup=>(i->i.<0.05)=>:oselective,
            :dup=>(i->i.<0.05)=>:dselective,
            :sfup=>(i->i.<0.05)=>:sfselective,
            :po=>:pd,
            :po=>(i->mod.(i,180))=>:po)
jldsave(joinpath(resultroot,"allorisfunit.jld2");allorisfunit)
allorisfunit = load(joinpath(resultroot,"allorisfunit.jld2"),"allorisfunit")
acolor = RGB(0.3,0.3,0.3)
lcolor = ColorMaps["lms_mccliso"].colors[end]
mcolor = ColorMaps["lms_mccmiso"].colors[end]
scolor = ColorMaps["lms_mccsiso"].colors[end]
aoricg = cgrad(fill(acolor,2))
loricg = cgrad(fill(lcolor,2))
moricg = cgrad(fill(mcolor,2))
soricg = cgrad(fill(scolor,2))

## orisf responsive
orisfcu = innerjoin(allorisfunit,allcu,on=[:siteid,:id])
leftjoin!(orisfcu,penetration,on=:siteid)
orisfcsu = subset(orisfcu,:good)
aorisfcsu = filter(r->r.c=='A',orisfcsu)
lorisfcsu = filter(r->r.c=='L',orisfcsu)
morisfcsu = filter(r->r.c=='M',orisfcsu)
sorisfcsu = filter(r->r.c=='S',orisfcsu)
aorisfrcsu = subset(aorisfcsu,:er)
lorisfrcsu = subset(lorisfcsu,:er)
morisfrcsu = subset(morisfcsu,:er)
sorisfrcsu = subset(sorisfcsu,:er)

allorisfcsu = combine(groupby(orisfcsu,[:id,:siteid]),
                :aligndepth=>first=>identity,
                :layer=>first=>identity,
                :er=>any=>:eresponsive,
                [:c,:er]=>((c,r)->join(c[r]))=>:RTYPE,
                [:c,:er,:ocm]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last,
                [:c,:er,:dcm]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last,
                [:c,:er,:msf]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last,
                [:c,:er,:po]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last,
                [:c,:er,:pd]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last,
                [:c,:er,:psf]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last,
                [:c,:er,:f1f0]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last,
                [:c,:er,:sctype]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last)
allorisfrcsu = subset(allorisfcsu,:eresponsive)

plotdepthhist(allorisfcsu;g=:eresponsive,layer)
plotdepthhist(allorisfrcsu;g=:RTYPE,layer)


## orisf tuning
plotorituning = (unit;type=:data,cg=aoricg,dir=nothing,figfmt=[".png"],layer=nothing,title="OriTuning") -> begin
    x = 0:0.04:2π # 2.3 deg
    if type == :fit
        tc = hcat(map(f->predict(f,x),unit.ofit)...)
    elseif type == :data
        tc = hcat(map((l,r)->Spline1D([deg2rad.(l);2π],[r;r[1]],k=1)(x),unit.ori,unit.orir)...)
    end
    tc./=maximum(tc,dims=1)
    p=plot(leg=false,proj=:polar,yticks=[],xticks=range(0,3π/2,length=4),xformatter=x->round(Int,rad2deg(x)),size=(500,500))
    for i in 1:size(tc,2)
        ccg = cgrad(map(j->coloralpha(color(get(cg,j/length(x))),tc[j,i]),eachindex(x)))
        @views plot!(p,x,tc[:,i],color=ccg,lz=x)
    end
    if !isnothing(layer)
        ann = [(1,mean(layer[k]),text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
    end
    isnothing(dir) && return p
    mkpath(dir);foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotsftuning = (unit;type=:data,cg=aoricg,dir=nothing,figfmt=[".png"],layer=nothing,title="SFTuning") -> begin
    x = 0.2:0.05:6.4
    if type == :fit
        tc = hcat(map(f->predict(f,x),unit.sffit)...)
    elseif type == :data
        tc = hcat(map((l,r)->Spline1D(l,r,k=1)(x),unit.sf,unit.sfr)...)
    end
    tc./=maximum(tc,dims=1)
    p=plot(leg=false,yticks=[],xticks=0.1 .* 2 .^ (1:6))
    for i in 1:size(tc,2)
        ccg = cgrad(map(j->coloralpha(color(get(cg,j/length(x))),tc[j,i]),eachindex(x)))
        @views plot!(p,x,tc[:,i],color=ccg,lz=x)
    end
    if !isnothing(layer)
        ann = [(1,mean(layer[k]),text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
    end
    isnothing(dir) && return p
    mkpath(dir);foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

foreach(siteid->plotorituning(filter(r->r.siteid==siteid,aorisfrcsu);
    cg=aoricg,dir=joinpath(resultroot,"Tuning"),type=:data,title="$(siteid)_A_OriTuning_data"),
    levels(aorisfrcsu.siteid))
foreach(siteid->plotsftuning(filter(r->r.siteid==siteid,aorisfrcsu);
    cg=aoricg,dir=joinpath(resultroot,"Tuning"),type=:data,title="$(siteid)_A_SFTuning_data"),
    levels(aorisfrcsu.siteid))

plotorituning(filter(r->r.siteid=="AG2_V1_ODL4",sorisfrcsu),cg=soricg,type=:fit)
plotsftuning(filter(r->r.siteid=="AG2_V1_ODL4",sorisfrcsu),cg=soricg,type=:data)


p=select(sorisfrcsu,[:ocm,:ocv,:dcm,:dcv,:ocor,:ofvu,:sfcor,:sffvu,:msf,:layer,:siteid]) |>
[@vlplot(:bar,y={"siteid"},x={"count()"});
@vlplot(:tick,y={"layer"},x={"ocv"});
@vlplot(:tick,y={"layer"},x={"ocm",axis={values=0:90:180}});
@vlplot(:tick,y={"layer"},x={"dcv"});
@vlplot(:tick,y={"layer"},x={"dcm",axis={values=0:90:360}});
@vlplot(:bar,y={"count()"},x={"ocor",bin={step=0.05}});
@vlplot(:bar,y={"count()"},x={"ofvu",bin={step=0.05}});
@vlplot(:bar,y={"count()"},x={"sfcor",bin={step=0.05}});
@vlplot(:bar,y={"count()"},x={"sffvu",bin={step=0.05}});
@vlplot(:tick,y={"layer"},x={"msf"})]


plotdepthscatter(aorisfrcsu;x=:ocv,layer,xlabel="ocv",leg=false,color=acolor)
plotdepthscatter(aorisfrcsu;x=:ocm,layer,xlabel="ocm",leg=false,color=acolor,xticks=0:90:180)
plotdepthscatter(aorisfrcsu;x=:dcv,layer,xlabel="dcv",leg=false,color=acolor)
plotdepthscatter(aorisfrcsu;x=:dcm,layer,xlabel="dcm",leg=false,color=acolor,xticks=0:90:360)
plotdepthscatter(aorisfrcsu;x=:f1f0,layer,xlabel="f1f0",leg=false,color=acolor)

plotdepthhist(aorisfrcsu;g=:sctype,layer)
plotdepthhist(sorisfrcsu;g=:sftype,layer)

p=subset(select(aorisfrcsu,[:oselective,:ocm,:ocv,:po,:ohw,:layer,:siteid,:aligndepth]),:oselective) |>
@vlplot(:point,y={"aligndepth",sort="descending"},x={"po",axis={values=0:90:180}},color={"layer",scale={scheme=:category10}},
    columns=9,wrap={"siteid"})

p=subset(select(aorisfrcsu,[:dselective,:dcm,:dcv,:pd,:ohw,:layer,:siteid,:aligndepth]),:dselective) |>
@vlplot(:point,y={"aligndepth",sort="descending"},x={"pd",axis={values=0:90:360}},color={"layer",scale={scheme=:category10}},
    columns=9,wrap={"siteid"})

p=subset(select(aorisfrcsu,[:sfselective,:msf,:psf,:layer,:siteid,:aligndepth]),:sfselective) |>
@vlplot(:point,y={"aligndepth",sort="descending"},x={"psf"},color={"layer",scale={scheme=:category10}},
    columns=9,wrap={"siteid"})


## Relation of OriSF between Spectral Channels
function plotorisfpair(unit,sp;w=350,h=350,sflim=nothing,f1f0lim=nothing)
    punit = filter(r->contains(r.RTYPE,sp.first) && contains(r.RTYPE,sp.second),unit)

    p = plot(layout=(4,1),leg=false,size=(w,4h))
    pn = :po
    x = getindex.(punit[!,pn],sp.first)
    y = getindex.(punit[!,pn],sp.second)
    xticks=yticks=0:45:180
    scatter!(p[1],x,y;xlabel="$(sp.first) (deg)",ylabel="$(sp.second) (deg)",title="$pn",
    xticks,yticks,ms=3,msw=0,ma=0.7,ratio=1)
    plot!(p[1],[0,180],[0,180],color=:gray30)

    pn = :pd
    x = getindex.(punit[!,pn],sp.first)
    y = getindex.(punit[!,pn],sp.second)
    xticks=yticks=0:45:360
    scatter!(p[2],x,y;xlabel="$(sp.first) (deg)",ylabel="$(sp.second) (deg)",title="$pn",
    xticks,yticks,ms=3,msw=0,ma=0.7,ratio=1)
    plot!(p[2],[0,360],[0,360],color=:gray30)

    pn = :psf
    x = getindex.(punit[!,pn],sp.first)
    y = getindex.(punit[!,pn],sp.second)
    lim = isnothing(sflim) ? max(maximum(x),maximum(y)) : sflim
    scatter!(p[3],x,y;xlabel="$(sp.first) (cycle/deg)",ylabel="$(sp.second) (cycle/deg)",title="$pn",
    ms=3,msw=0,ma=0.7,ratio=1)
    plot!(p[3],[0,lim],[0,lim],color=:gray30)

    pn = :f1f0
    x = log2.(getindex.(punit[!,pn],sp.first))
    y = log2.(getindex.(punit[!,pn],sp.second))
    lim = isnothing(f1f0lim) ? max(maximum(abs.(x)),maximum(abs.(y))) : f1f0lim
    scatter!(p[4],x,y;xlabel="$(sp.first)",ylabel="$(sp.second)",title="Log₂(F1/F0)",
    ms=3,msw=0,ma=0.7,ratio=1)
    plot!(p[4],[-lim,lim],[-lim,lim],color=:gray30)
    vline!(p[4],[0],color=:gray30)
    hline!(p[4],[0],color=:gray30)
end

plotorisfpair(allorisfrcsu,'A'=>'L')
plot((plotorisfpair(allorisfrcsu,i=>j,sflim=6,f1f0lim=9) for (i,j) in combinations(['A','L','M','S'],2))...,
            layout=(1,6),size=(6*350,4*350))









## non-sta but orisf responsive
stanrcsu = load(joinpath(resultroot,"stanrcsu.jld2"),"stanrcsu")
nstaorisfrcsu = innerjoin(select(stanrcsu,[:id,:roi]),allorisfrcsu,on=:id)
anstaorisfrcsu = innerjoin(select(stanrcsu,[:id,:roi]),aorisfrcsu,on=:id)
lnstaorisfrcsu = innerjoin(select(stanrcsu,[:id,:roi]),lorisfrcsu,on=:id)
mnstaorisfrcsu = innerjoin(select(stanrcsu,[:id,:roi]),morisfrcsu,on=:id)
snstaorisfrcsu = innerjoin(select(stanrcsu,[:id,:roi]),sorisfrcsu,on=:id)

plotdepthhist(nstaorisfrcsu;g=:RTYPE,layer)


plotdepthscatter(anstaorisfrcsu;x=:ocv,layer,xlabel="ocv",leg=false,color=acolor)
plotdepthscatter(anstaorisfrcsu;x=:dcv,layer,xlabel="dcv",leg=false,color=acolor)
plotdepthhist(anstaorisfrcsu;g=:sftype,layer)
plotdepthhist(anstaorisfrcsu;g=:sctype,layer,palette=[:tomato,:deepskyblue])





batchunit = filter(r->r.siteid=="AG1_V1_ODL1"&&r.layer=="3",anstaorisfrcsu)


batchunit.roi[3]=(centerdeg=[1.1,1.0],radiideg=(0.8,0.8),radiusdeg=0.8)

batchunit=batchunit[3:3,:]



batchunit=nothing


savefig("test.png")
save("test.png",p)


