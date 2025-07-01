using NeuroAnalysis,DataFrames,JLD2,StatsBase,StatsPlots,VegaLite,XLSX,Combinatorics,DataStructures,CircStats,HypothesisTests,AlgebraOfGraphics

function collectflash2color_od(indir;unit=DataFrame(),datafile="flash2colordataset_od.jld2",eye=["N","D"])
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile),"dataset")
            siteid = dataset["siteid"]
            id = ["$(siteid)_$(g ? "S" : "M")U$u" for (u,g) in dataset["ugood"]]
            ei = indexin(eye,dataset["eye"])

            df = DataFrame(:siteid=>siteid,:id=>id,
                    :ccode => Tuple(eye),
                    :cresponsive=>[v[ei] for v in values(dataset["responsive"])],
                    :cmodulative=>[v[ei] for v in values(dataset["modulative"])],
                    :cenoughresponse=>[v[ei] for v in values(dataset["enoughresponse"])],
                    :cr=>[map(i->i.r,v[ei]) for v in values(dataset["frs"])],
                    :cm=>[map(i->i.m,v[ei]) for v in values(dataset["frs"])],
                    :cse=>[map(i->i.se,v[ei]) for v in values(dataset["frs"])],
                    :fa=>Tuple(dataset["fa"][ei]),
                    :fac=>Tuple(dataset["fac"][ei]),
                    :fac_od=>Tuple(dataset["fac_od"][ei]))
            append!(unit,df,cols=:union)
        end
    end
    return unit
end

resultroot = "Z:/"
f2cunit = collectflash2color_od(resultroot)
jldsave(joinpath(resultroot,"f2cunit_od.jld2");f2cunit)



## cortical unit Flash2Color responsive
cofdcode = Dict("DKL_X"=>"A","DKL_Y"=>"Y","DKL_Z"=>"S","LMS_X"=>"L","LMS_Y"=>"M","LMS_Z"=>"S","DKL_RG"=>"L","DKL_BY"=>"S")

function colorminmax(x,c)
    ismissing(x) && return x
    if c=="A"
        s = gray.(Gray.(values(x)))
    elseif c=="L"
        s = comp1.(LMS.(values(x)))
    elseif c=="M"
        s = comp2.(LMS.(values(x))) 
    elseif c=="S"
        s = comp3.(LMS.(values(x))) 
    end
    s[2] > s[1]
end

function loadcofdi!(cofdi,penetration;s=50)
    d = cofdi["scofdis"][findfirst(i->i==s,cofdi["srange"])]
    d = select(d,:siteid,Not([:siteid,:subject]).=>ByRow(i->ismissing(i) ? i : 2*i[1]-1).=>identity)
    cs = cofdi["color_cofd"]
    cminmax = Dict(c=>colorminmax.(cs[!,c],cofdcode[c]) for c in filter(i->i ∉ ["subject","siteid"],names(cs)))
    foreach(c-> d[!,c] = map((i,j)->ismissing(i) ? i : (j ? i : -i),d[!,c],cminmax[c]),keys(cminmax))
    df = select(d,:siteid,:DKL_X=>:cofdi_A,[:LMS_X,:DKL_RG]=>ByRow((i,j)->ismissing(i) ? j : i)=>:cofdi_L,:LMS_Y=>:cofdi_M,
                [:LMS_Z,:DKL_BY]=>ByRow((i,j)->ismissing(i) ? j : i)=>:cofdi_S)
    df = transform(df,Not(:siteid).=>ByRow(abs).=>(i->"a"*i))
    leftjoin!(penetration,df,on=:siteid)
end
function loadodi!(odi,penetration;s=80)
    d = odi["sodis"][:,findfirst(i->i==s,odi["srange"])]
    df = DataFrame(siteid=odi["siteid"],odi = 2*d .- 1)
    df = transform(df,Not(:siteid).=>ByRow(abs).=>(i->"a"*i))
    leftjoin!(penetration,df,on=:siteid)
end



f2cunit = load(joinpath(resultroot,"f2cunit_od.jld2"),"f2cunit")
layer,nlbt = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate","nlbt")
layercg = cgrad(:tab10,nlbt,categorical=true)
allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
allcsu = subset(allcu,:good)
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
odi = load(joinpath(resultroot,"odi.jld2"))
loadodi!(odi,penetration)
cofdi = load(joinpath(resultroot,"cofdi.jld2"))
loadcofdi!(cofdi,penetration)
select!(penetration,Not([:Subject_ID,:RecordSession,:RecordSite,:coord,:d2b,:layer]),:odi=>ByRow(i->ismissing(i) ? i : i>0 ? "Right" : "Left")=>:od,
        :cofd=>ByRow(i->ismissing(i) ? i : first.(split(i,", ")) |> join ∘ sort)=>:cofd_type)

f2ccu = innerjoin(allcu,f2cunit,on=[:siteid,:id])
excludesites = ["AG1_V1_ODR8","AG1_V1_ODL17","AG2_V1_ODL18"]
filter!(r->r.siteid ∉ excludesites,f2ccu)
filter!(r->r.siteid ∉ excludesites,penetration)
leftjoin!(f2ccu,penetration,on=:siteid)

transform!(f2ccu,[:cresponsive,:ccode]=>ByRow((r,c)->all(.!r) ? missing : join(c[findall(r)]))=>:RTYPE)
transform!(f2ccu,[:cmodulative,:ccode]=>ByRow((r,c)->all(.!r) ? missing : join(c[findall(r)]))=>:MTYPE)
f2ccu.responsive = map(i->ismissing(i) ? false : true,f2ccu.RTYPE)
f2ccu.modulative = map(i->ismissing(i) ? false : true,f2ccu.MTYPE)
f2ccsu = subset(f2ccu,:good)
f2crcsu = dropmissing(f2ccsu,:RTYPE)
f2cmcsu = dropmissing(f2ccsu,:MTYPE)
f2cnrcsu = subset(f2ccsu,:responsive=>.!)
f2cdir = joinpath(resultroot,"Flash2Color_OD");mkpath(f2cdir)
figfmt = [".svg",".png"]


import CairoMakie as mk
function plotdepthhist(unit;g=:responsive,dir=nothing,figfmt=[".png"],layer=nothing,xlabel="Number of Units",ts="",
                siteid="All",ylabel="Normalized Cortical Depth",size=(450,750),palette=:tab20,binedges=0:0.02:1)
    sunit = siteid=="All" ? unit : filter(r->r.siteid==siteid,unit)
    ts = isempty(ts) ? ts : "_"*ts
    title="$(siteid)_DepthHist_$g$ts"
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

function plotdepthscatter(unit;x=:Diameter,xfun=first,upv=nothing,g=:rc,dir=nothing,siteid="All",figfmt=[".png"],layer=nothing,xlabel="$x",ts="",
            ylabel="Normalized Cortical Depth",size=(450,750),palette=:tab20,ma=-0.7,
            gwfun=x->isempty(x) ? 0 : median(x),gwbfun=(x->isempty(x) ? 0 : percentile(x,25),x->isempty(x) ? 0 : percentile(x,75)))
    sunit = siteid=="All" ? unit : filter(r->r.siteid==siteid,unit)
    if isnothing(upv)
        pvstr=""
    else
        sunit = filter(r->r[upv.first]==upv.second,sunit)
        pvstr="_$(upv.first)=$(upv.second)"
    end
    ts = isempty(ts) ? ts : "_"*ts
    title="$(siteid)$(pvstr)_Depth-$(x)_$g$ts"
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


plotdepthhist(f2crcsu;g=:spiketype,layer)
# plotdepthhist(f2crcsu;g=:spiketype,layer,dir=f2cdir,figfmt)

plotdepthhist(f2ccsu;g=:responsive,layer)
plotdepthhist(f2ccsu;g=:modulative,layer)
plotdepthhist(f2crcsu;g=:RTYPE,layer)
plotdepthhist(f2cmcsu;g=:MTYPE,layer)

# plotdepthhist(f2ccsu;g=:responsive,layer,dir=f2cdir,figfmt)
# plotdepthhist(f2ccsu;g=:modulative,layer,dir=f2cdir,figfmt)
# plotdepthhist(f2crcsu;g=:RTYPE,layer,dir=f2cdir,figfmt)
# plotdepthhist(f2cmcsu;g=:MTYPE,layer,dir=f2cdir,figfmt)

# pl = select(f2crcsu,[:RTYPE,:siteid,:layer]) |>
#     [@vlplot(:bar,y={"layer"},x={"count()",axis={title="Number of Units",grid=false}},color={"RTYPE:n",scale={scheme=:category20}});
#     @vlplot(:bar,y={"siteid"},x={"count()",axis={title="Number of Units",grid=false}},color={"RTYPE:n",scale={scheme=:category20}})]
# foreach(ext->save(joinpath(f2cdir,"Layer_SiteID_Hist_RTYPE$ext"),pl),figfmt)



## Flash2Color tuning properties across layers
acolor = RGB(0.3,0.3,0.3)
scc = Dict("A"=>:gray,"T"=>:orange,"N"=>:darkblue,"M"=>:darkgreen,"W"=>:darkred)
cc = Dict("N"=>:darkblue,"D"=>:darkred,"ON"=>:gray60,"OFF"=>:gray40)


function getfrf(unit;csc=["N_A-","N_A+","D_A-","D_A+"],tlim=6)
    df = select(unit,:id,[:cmodulative,:fac_od,:cm,:od,:cr]=>ByRow((mo,cc,cm,od,cr)->begin
    
    fac_od = vcat(cc...)
    m = vcat(cm...)
    r = vcat(cr...)
    ii= indexin(csc,fac_od)
    m=m[ii];r=r[ii]

    rc_N = m[2]/m[1]
    rc_D = m[4]/m[3]
    rd_ON = m[4]/m[2]
    rd_OFF = m[3]/m[1]
    lrc_N = log2(rc_N)
    lrc_D = log2(rc_D)
    lrd_ON = log2(rd_ON)
    lrd_OFF = log2(rd_OFF)
    
    pc_N=(m[2]-m[1])/(m[2]+m[1])
    pc_D=(m[4]-m[3])/(m[4]+m[3])
    pd_ON=(m[4]-m[2])/(m[4]+m[2])
    pd_OFF=(m[3]-m[1])/(m[3]+m[1])

    hc_N=UnequalVarianceTTest(r[2],r[1])
    hc_D=UnequalVarianceTTest(r[4],r[3])
    hd_ON=UnequalVarianceTTest(r[4],r[2])
    hd_OFF=UnequalVarianceTTest(r[3],r[1])
    tc_N = hc_N.t
    tc_D = hc_D.t
    td_ON = hd_ON.t
    td_OFF = hd_OFF.t
    tc_N = clamp(tc_N,-tlim,tlim)
    tc_D = clamp(tc_D,-tlim,tlim)
    td_ON = clamp(td_ON,-tlim,tlim)
    td_OFF = clamp(td_OFF,-tlim,tlim)

    if od=="Right"
        lre_ON = lrd_ON
        lre_OFF = lrd_OFF
        pe_ON = pd_ON
        pe_OFF = pd_OFF
        te_ON = td_ON
        te_OFF = td_OFF
    else
        lre_ON = -lrd_ON
        lre_OFF = -lrd_OFF
        pe_ON = -pd_ON
        pe_OFF = -pd_OFF
        te_ON = -td_ON
        te_OFF = -td_OFF
    end

    (;lrc_N,lrc_D,lrd_ON,lrd_OFF,pc_N,pc_D,pd_ON,pd_OFF,tc_N,tc_D,td_ON,td_OFF,lre_ON,lre_OFF,pe_ON,pe_OFF,te_ON,te_OFF)
    end)=>AsTable)
    df = transform(df,Not(:id).=>ByRow(abs)=>(i->"a"*i))
    df=leftjoin(unit,df,on=:id)
end

cfrf = getfrf(f2crcsu)


plotdepthscatter(cfrf;layer,x=:lrc_D,g=:responsive)
plotdepthscatter(cfrf;layer,x=:lrc_N,g=:responsive)
plotdepthscatter(cfrf;layer,x=:pc_D,g=:responsive)
plotdepthscatter(cfrf;layer,x=:pc_N,g=:responsive)
plotdepthscatter(cfrf;layer,x=:tc_D,g=:responsive)
plotdepthscatter(cfrf;layer,x=:tc_N,g=:responsive)

plotdepthscatter(cfrf;layer,x=:lrd_ON,g=:responsive)
plotdepthscatter(cfrf;layer,x=:lrd_OFF,g=:responsive)
plotdepthscatter(cfrf;layer,x=:pd_ON,g=:responsive)
plotdepthscatter(cfrf;layer,x=:pd_OFF,g=:responsive)
plotdepthscatter(cfrf;layer,x=:td_ON,g=:responsive)
plotdepthscatter(cfrf;layer,x=:td_OFF,g=:responsive)

plotdepthscatter(cfrf;layer,x=:lre_ON,g=:responsive)
plotdepthscatter(cfrf;layer,x=:lre_OFF,g=:responsive)
plotdepthscatter(cfrf;layer,x=:pe_ON,g=:responsive)
plotdepthscatter(cfrf;layer,x=:pe_OFF,g=:responsive)
plotdepthscatter(cfrf;layer,x=:te_ON,g=:responsive)
plotdepthscatter(cfrf;layer,x=:te_OFF,g=:responsive)

plotlayerscatter(cfrf;y=:lrc_D,x=:cofdi_A)
plotlayerscatter(cfrf;y=:lrc_N,x=:cofdi_A)
plotlayerscatter(cfrf;y=:pc_D,x=:cofdi_A)
plotlayerscatter(cfrf;y=:pc_N,x=:cofdi_A)
plotlayerscatter(cfrf;y=:tc_D,x=:cofdi_A)
plotlayerscatter(cfrf;y=:tc_N,x=:cofdi_A)

plotlayerscatter(cfrf;y=:lrd_ON,x=:odi)
plotlayerscatter(cfrf;y=:lrd_OFF,x=:odi)
plotlayerscatter(cfrf;y=:pd_ON,x=:odi)
plotlayerscatter(cfrf;y=:pd_OFF,x=:odi)
plotlayerscatter(cfrf;y=:td_ON,x=:odi)
plotlayerscatter(cfrf;y=:td_OFF,x=:odi)

plotlayerscatter(cfrf;y=:alrd_ON,x=:aodi)
plotlayerscatter(cfrf;y=:alrd_OFF,x=:aodi)
plotlayerscatter(cfrf;y=:apd_ON,x=:aodi)
plotlayerscatter(cfrf;y=:apd_OFF,x=:aodi)

plotlayerscatter(cfrf;y=:lre_ON,x=:odi)
plotlayerscatter(cfrf;y=:lre_OFF,x=:odi)
plotlayerscatter(cfrf;y=:pe_ON,x=:odi)
plotlayerscatter(cfrf;y=:pe_OFF,x=:odi)
plotlayerscatter(cfrf;y=:te_ON,x=:odi)
plotlayerscatter(cfrf;y=:te_OFF,x=:odi)

plotlayerscatter(cfrf;y=:lrd_ON,x=:aodi)
plotlayerscatter(cfrf;y=:lrd_OFF,x=:aodi)
plotlayerscatter(cfrf;y=:pd_ON,x=:aodi)
plotlayerscatter(cfrf;y=:pd_OFF,x=:aodi)

plotlayerscatter(cfrf;y=:psi,x=:cofdi_A,g=:spiketype,msca=0)

function plotlayerscatter(unit;x=:ori_lhi,y=:hw,dir=nothing,figfmt=[".png"],size=(200,250),cc=scc,g=nothing,ts="",
    xlims=nothing,ylims=nothing,nol1=true,noatspike=true,markersize=4,strokewidth=1,msca=0.7,color=acolor)
    nol1 && (unit = filter(r->r.layer≠"1",unit))
    noatspike && :spiketype in propertynames(unit) && (unit = filter(r->r.spiketype ∉ ["A","T"],unit))
    unit = filter(r->isfinite(r[x]) & isfinite(r[y]),unit)
    nl = levels(unit.layer) |> length
    ts = isempty(ts) ? ts : "_$ts"
    gs = isnothing(g) ? "" : "_$g"
    title = "LayerScatter_$x-$(y)_$gs$ts"
    plt=data(unit)*mapping(x,y,col=:layer)
    if isnothing(g)
    ng=1.07
    else
    ng=0.8length(levels(unit[!,g]));plt*=mapping(color=g,row=g)
    end
    plt*=(linear()+visual(mk.Scatter;markersize,strokewidth,color=RGBA(0,0,0,0),strokecolor=(color,msca)))
    f = draw(plt,palettes=(color=collect(pairs(cc)),),axis=(limits=(xlims,ylims),),figure=(size=(size[1]*nl,ng*size[2]),))
    isnothing(dir) ? f : foreach(ext->save(joinpath(dir,"$title$ext"),f),figfmt)
end

function plotlayerdist(unit;x=:osi,dir=nothing,figfmt=[".png"],size=(900,600),cc=cc,plottype=:box,g=:Color,
    nol1=true,noatspike=true,ts="")
    nol1 && (unit = filter(r->r.layer≠"1",unit))
    noatspike && :spiketype in propertynames(unit) && (unit = filter(r->r.spiketype ∉ ["A","T"],unit))
    unit = filter(r->!isnan(r[x]),unit)
    ts = isempty(ts) ? ts : "_$ts"
    gs = isnothing(g) ? "" : "_$g"
    if plottype == :box
        p = data(unit)*mapping(:layer,x,dodge=g,color=g)*visual(mk.BoxPlot,show_outliers=false,gap=0.3)
        # f=draw(p,scales(Color=(;palette=collect(pairs(cc)))),figure=(;size=size))
        f=draw(p,palettes=(color=collect(pairs(cc)),),figure=(size=size,))
        title = "LayerBox-$x$gs$ts"
        elseif plottype == :violin
        p = data(unit)*mapping(:layer,x,dodge=g,color=g)*visual(mk.Violin,datalimits=extrema,gap=0.3)
        f=draw(p,palettes=(color=collect(pairs(cc)),),figure=(size=size,))
        title = "LayerViolin-$x$gs$ts"
        elseif plottype==:density
        p = data(unit)*mapping(x,col=:layer,color=g)*AlgebraOfGraphics.density()
        f=draw(p,palettes=(color=collect(pairs(cc)),),figure=(size=size,))
        title = "LayerDensity-$x$gs$ts"
    end
    isnothing(dir) ? f : foreach(ext->save(joinpath(dir,"$title$ext"),f),figfmt)
end


for c in ["N","D"]
    # plotdepthscatter(cfrf;layer,x="lrc_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf;layer,x="pc_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf;layer,x="tc_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    
    plotlayerscatter(cfrf;y="lrc_$c",x=:cofdi_A,dir=f2cdir,figfmt,ts=c,color=cc[c])
    plotlayerscatter(cfrf;y="pc_$c",x=:cofdi_A,dir=f2cdir,figfmt,ts=c,color=cc[c])
    plotlayerscatter(cfrf;y="tc_$c",x=:cofdi_A,dir=f2cdir,figfmt,ts=c,color=cc[c])
end

for c in ["ON","OFF"]
    # plotdepthscatter(cfrf;layer,x="lrd_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf;layer,x="pd_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf;layer,x="td_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf;layer,x="lre_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf;layer,x="pe_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf;layer,x="te_$c",g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])

    plotlayerscatter(cfrf;y="lre_$c",x=:odi,dir=f2cdir,figfmt,ts=c,color=cc[c])
    plotlayerscatter(cfrf;y="pe_$c",x=:odi,dir=f2cdir,figfmt,ts=c,color=cc[c])
    plotlayerscatter(cfrf;y="te_$c",x=:odi,dir=f2cdir,figfmt,ts=c,color=cc[c])
end




jldsave(joinpath(resultroot,"f2crcsu.jld2");f2crcsu=cfrf)

