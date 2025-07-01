using NeuroAnalysis,DataFrames,JLD2,StatsBase,StatsPlots,VegaLite,XLSX,Combinatorics,DataStructures,CircStats,HypothesisTests,AlgebraOfGraphics

function collectflash2color(indir;unit=DataFrame(),datafile="flash2colordataset.jld2",cch=["A","Y","S"])
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
                    :cr=>[map(i->i.r,v[ci]) for v in values(dataset["frs"])],
                    :cm=>[map(i->i.m,v[ci]) for v in values(dataset["frs"])],
                    :cse=>[map(i->i.se,v[ci]) for v in values(dataset["frs"])],
                    :fa=>Tuple(dataset["fa"][ci]),
                    :fac=>Tuple(dataset["fac"][ci]))
            append!(unit,df,cols=:union)
        end
    end
    return unit
end

resultroot = "Z:/"
f2cunit = collectflash2color(resultroot)
jldsave(joinpath(resultroot,"f2cunit.jld2");f2cunit)



## cortical unit Flash2Color responsive
cofdcode = Dict("DKL_X"=>"A","DKL_Z"=>"S","LMS_X"=>"L","LMS_Y"=>"M","LMS_Z"=>"S","DKL_RG"=>"L","DKL_BY"=>"S")

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


f2cunit = load(joinpath(resultroot,"f2cunit.jld2"),"f2cunit")
layer,nlbt = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate","nlbt")
layercg = cgrad(:tab10,nlbt,categorical=true)
allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
allcsu = subset(allcu,:good)
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
cofdi = load(joinpath(resultroot,"cofdi.jld2"))
loadcofdi!(cofdi,penetration)
select!(penetration,Not([:Subject_ID,:RecordSession,:RecordSite,:coord,:d2b,:layer]),:cofd=>ByRow(i->ismissing(i) ? i : first.(split(i,", ")) |> join ∘ sort)=>:cofd_type)

f2ccu = innerjoin(allcu,f2cunit,on=[:siteid,:id])
excludesites = ["AG1_V1_ODR8","AG1_V1_ODL17","AG2_V1_ODL18"]
filter!(r->r.siteid ∉ excludesites,f2ccu)
filter!(r->r.siteid ∉ excludesites,penetration)
leftjoin!(f2ccu,penetration,on=:siteid)

transform!(f2ccu,:cofd_type=>ByRow(i->any(occursin.(['L','M','S'],i)))=>:cofd_C,
        :cofd_type=>ByRow(i->any(occursin.(['L','M'],i)))=>:cofd_LM,:cofd_type=>ByRow(i->contains(i,'S'))=>:cofd_S,
        :cofd_type=>ByRow(i->contains(i,'A'))=>:cofd_A,:cofd_type=>ByRow(i->i=="N")=>:cofd_N)
transform!(f2ccu,[:cofd_N,:cofd_C]=>ByRow((n,c)->n ? 'N' : c ? 'C' : 'A')=>:cofd_NAC)
transform!(f2ccu,[:cofd_N,:cofd_C,:cofd_LM]=>ByRow((n,c,l)->n ? 'N' : (c ? (l ? 'L' : 'S') : 'A'))=>:cofd_NAS)

transform!(f2ccu,[:cresponsive,:ccode]=>ByRow((r,c)->all(.!r) ? missing : join(c[findall(r)]))=>:RTYPE)
transform!(f2ccu,[:cmodulative,:ccode]=>ByRow((r,c)->all(.!r) ? missing : join(c[findall(r)]))=>:MTYPE)
f2ccu.responsive = map(i->ismissing(i) ? false : true,f2ccu.RTYPE)
f2ccu.modulative = map(i->ismissing(i) ? false : true,f2ccu.MTYPE)
f2ccu.CTYPE = map(i->ismissing(i) ? missing : replace(i,r"A([Y,S]+)"=>s"\1"),f2ccu.RTYPE)
f2ccu.CMTYPE = map(i->ismissing(i) ? missing : replace(i,r"A([Y,S]+)"=>s"\1"),f2ccu.MTYPE)
f2ccsu = subset(f2ccu,:good)
f2crcsu = dropmissing(f2ccsu,:RTYPE)
f2cmcsu = dropmissing(f2ccsu,:MTYPE)
f2cnrcsu = subset(f2ccsu,:responsive=>.!)
f2cdir = joinpath(resultroot,"Flash2Color");mkpath(f2cdir)
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
plotdepthhist(f2crcsu;g=:CTYPE,layer)
plotdepthhist(f2cmcsu;g=:MTYPE,layer)
plotdepthhist(f2cmcsu;g=:CMTYPE,layer)

# plotdepthhist(f2ccsu;g=:responsive,layer,dir=f2cdir,figfmt)
# plotdepthhist(f2ccsu;g=:modulative,layer,dir=f2cdir,figfmt)
# plotdepthhist(f2crcsu;g=:RTYPE,layer,dir=f2cdir,figfmt)
# plotdepthhist(f2crcsu;g=:CTYPE,layer,dir=f2cdir,figfmt)
# plotdepthhist(f2cmcsu;g=:MTYPE,layer,dir=f2cdir,figfmt)
# plotdepthhist(f2cmcsu;g=:CMTYPE,layer,dir=f2cdir,figfmt)

# pl = select(f2crcsu,[:RTYPE,:siteid,:layer]) |>
#     [@vlplot(:bar,y={"layer"},x={"count()",axis={title="Number of Units",grid=false}},color={"RTYPE:n",scale={scheme=:category20}});
#     @vlplot(:bar,y={"siteid"},x={"count()",axis={title="Number of Units",grid=false}},color={"RTYPE:n",scale={scheme=:category20}})]
# foreach(ext->save(joinpath(f2cdir,"Layer_SiteID_Hist_RTYPE$ext"),pl),figfmt)


plt = data(f2crcsu)*frequency()*mapping(:RTYPE,:layer,col=:cofd_NAC)
f=draw(plt,figure=(;size=(600,500)))
foreach(ext->save(joinpath(orisfdir,"Penetration_ORI_LHI-LR$ext"),f),figfmt)



## Flash2Color tuning properties across layers
acolor = RGB(0.3,0.3,0.3)
ycolor = ColorMaps["dkl_mcclmiso"].colors[end]
scolor = ColorMaps["dkl_mccslmiso"].colors[end]
cc = Dict('A'=>acolor,'Y'=>ycolor,'S'=>scolor)
scc = Dict("A"=>:gray,"T"=>:orange,"N"=>:darkblue,"M"=>:darkgreen,"W"=>:darkred)
lhcc = Dict("L"=>:lightblue,"H"=>:pink)
tfcc = Dict(true=>:pink,false=>:lightblue)
dcc = Dict('A'=>:gray,'C'=>:pink,'N'=>:lightblue)
scc = Dict('A'=>:gray,'L'=>:pink,'S'=>:mediumpurple,'N'=>:lightblue)
fmean(x) = mean(filter(isfinite,x))
cmean(x) = angle(sum(cis.(x)))
orim(o;n=2)=rad2deg(mod2pi(cmean(n*deg2rad.(o)))/n)
oricv(o;n=2)=circ_var(n*deg2rad.(o)).S
orir(o;n=2) = circ_r(n*deg2rad.(o))
orid(o1,o2;n=2) = rad2deg.(circ_dist(n*deg2rad.(o1),n*deg2rad.(o2))/n)
absorid(o1,o2;n=2) = abs.(orid(o1,o2;n))
sfm(x,w=ones(size(x))) = 2^(sum(w.*log2.(x))/sum(w))
sfbw(h,l) = log2.(h./l)
aoricg = cgrad(fill(acolor,2))
loricg = cgrad(fill(lcolor,2))
moricg = cgrad(fill(mcolor,2))
soricg = cgrad(fill(scolor,2))

function getfrf(unit;c='A',cch="AYS",tlim=6)
    ci = findfirst(c,cch)
    df = transform(unit,[:cmodulative,:fac,:cm,:cr]=>ByRow((mo,cc,cm,cr)->begin
    m = cm[ci][indexin(["$c-","$c+"] ,cc[ci])]
    r = cr[ci][indexin(["$c-","$c+"] ,cc[ci])]

    lr = log2(m[2]/m[1])
    alr = abs(lr)

    psi=(m[2]-m[1])/sum(m)
    apsi = abs(psi)

    h=UnequalVarianceTTest(r[2],r[1])
    t = clamp(h.t,-tlim,tlim)

    maxi = argmax(cm[ci])
    maxm = cm[ci][maxi]
    pc = mo[ci] ? cc[ci][maxi] : "NS"
    (;maxm,lr,alr,psi,apsi,t,pc)
    end)=>AsTable)
    df = transform(groupby(df,:siteid),[:alr,:apsi].=>fmean.=>[:malr,:mapsi])
end

cfrf = Dict(c=>getfrf(f2crcsu;c) for c in ['A','Y','S'])

cfrf['A']




ccpdf = mapreduce(c->insertcols(
    combine(groupby(cfrf[c],:siteid),[:ori_lhi,:ori_lr,:ori_lhi_LH,:ori_lr_LH,:malr,:mapsi,
        :cofd_type,:cofd_N,:cofd_NAC,:cofd_NAS].=>first.=>identity),
        :Color=>c),append!,['A','Y','S'])

c='A'
pdf = combine(groupby(cfrf[c],:siteid),[:malr,:mapsi,
        :cofd_type,:cofd_N,:cofd_NAC,:cofd_NAS].=>first.=>identity)

plt=data(pdf)*mapping(:cofd_NAC,:malr,color=:cofd_NAC)*visual(mk.BoxPlot,show_outliers=false,gap=0.3)
f=draw(plt,palettes=(;color=collect(pairs(dcc))))
foreach(ext->save(joinpath(f2cdir,"Penetration_malr_Dist_$c$ext"),f),figfmt)

plt=data(pdf)*mapping(:cofd_NAS,:mapsi,color=:cofd_NAS)*visual(mk.BoxPlot,show_outliers=false,gap=0.3)
f=draw(plt,palettes=(;color=collect(pairs(scc))))
# f=draw(plt)
foreach(ext->save(joinpath(f2cdir,"Penetration_mapsi_Dist_$c$ext"),f),figfmt)



plt=data(ccpdf)*mapping(:cofd_NAS,:mapsi,color=:cofd_NAS,col=:Color)*visual(mk.BoxPlot,show_outliers=false,gap=0.3)
f=draw(plt,palettes=(;color=collect(pairs(scc))))
foreach(ext->save(joinpath(f2cdir,"Penetration_mapsi_Dist_$c$ext"),f),figfmt)

plt=data(ccpdf)*mapping(:cofd_NAC,:mapsi,color=:cofd_NAC,col=:Color)*visual(mk.BoxPlot,show_outliers=false,gap=0.3)
f=draw(plt,palettes=(;color=collect(pairs(dcc))))











# df = DataFrame(ori_lhi=repeat(pdf.ori_lhi,outer=2),ori_lr=repeat(pdf.ori_lr,outer=2),ori_lhi_spo=repeat(pdf.ori_lhi_spo,outer=2),
# ad=[pdf.dps_os;pdf.dsgi_os],pair=repeat(["ISI-L1~3","L1~3-L4~6"],inner=nrow(pdf)))
# filter!(r->!isnan(r.ad),df)
# plt = data(df)*mapping(:ori_lhi,:ad,color=:pair)*(linear()+visual(mk.Scatter,color=:gray40))
# f=draw(plt)
# foreach(ext->save(joinpath(orisfdir,"Penetration_ORI_LHI-AD$ext"),f),figfmt)

# plt = data(pdf)*mapping(:ori_lhi_spo,:ori_lhi)*(linear()+visual(mk.Scatter,color=:gray40))
# f=draw(plt,figure=(size=(600,500),))
# foreach(ext->save(joinpath(orisfdir,"Penetration_ORI_LHI_SPO-LHI$ext"),f),figfmt)

c='A'

"cofd_$c"

plotdepthhist(cfrf['A'];g=:pc,layer)
plotdepthscatter(cfrf['A'];layer,x=:lr,g=:responsive)
plotdepthscatter(cfrf['A'];layer,x=:psi,g=:responsive)
plotdepthscatter(cfrf['A'];layer,x=:t,g=:responsive)
plotdepthscatter(cfrf['A'];layer,x=:apsi,g=:responsive)

plotlayerscatter(cfrf['A'];y=:psi,x=:cofdi_A,ts='A')
plotlayerscatter(cfrf['Y'];y=:psi,x=:cofdi_L,ts='A')
plotlayerscatter(cfrf['A'];y=:psi,x=:cofdi_A,ts='A',g=:spiketype,msca=0)


for c in keys(cfrf)
    # plotdepthhist(cfrf[c];layer,g=:pc,dir=f2cdir,figfmt,ts=c)
    # plotdepthscatter(cfrf[c];layer,x=:lr,g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:psi,g=:responsive,dir=f2cdir,figfmt,ts=c,palette=[cc[c]])
    # plotdepthscatter(cfrf[c];layer,x=:alr,g=:cofd_N,dir=f2cdir,figfmt,ts=c,palette=[cc[c],RGB(0,0,1)])
    # plotdepthscatter(cfrf[c];layer,x=:alr,g=:cofd_NAC,dir=f2cdir,figfmt,ts=c,palette=[cc[c],RGB(1,0,1),RGB(0,0,1)])
    # plotdepthscatter(cfrf[c];layer,x=:apsi,g=:cofd_N,dir=f2cdir,figfmt,ts=c,palette=[cc[c],RGB(0,0,1)])
    # plotdepthscatter(cfrf[c];layer,x=:apsi,g=:cofd_NAC,dir=f2cdir,figfmt,ts=c,palette=[cc[c],RGB(1,0,1),RGB(0,0,1)])

    plotlayerscatter(cfrf[c];y=:lr,x="cofdi_$(c=='Y' ? 'L' : c)",dir=f2cdir,figfmt,ts=c,color=cc[c])
    plotlayerscatter(cfrf[c];y=:psi,x="cofdi_$(c=='Y' ? 'L' : c)",dir=f2cdir,figfmt,ts=c,color=cc[c])
    plotlayerscatter(cfrf[c];y=:psi,x="cofdi_$(c=='Y' ? 'L' : c)",dir=f2cdir,figfmt,ts=c,g=:spiketype,msca=0,color=cc[c])
    plotlayerscatter(cfrf[c];y=:alr,x="acofdi_$(c=='Y' ? 'L' : c)",dir=f2cdir,figfmt,ts=c,color=cc[c])
    plotlayerscatter(cfrf[c];y=:apsi,x="acofdi_$(c=='Y' ? 'L' : c)",dir=f2cdir,figfmt,ts=c,color=cc[c])
    plotlayerscatter(cfrf[c];y=:apsi,x="acofdi_$(c=='Y' ? 'L' : c)",dir=f2cdir,figfmt,ts=c,g=:spiketype,msca=0,color=cc[c])
end

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



plt = data(cfrf['A'])*mapping(:dmpo_os,color=:ori_lhi_LH,dodge=:ori_lhi_LH)*AlgebraOfGraphics.histogram(bins=-90:15:90)
f=draw(plt)
foreach(ext->save(joinpath(orisfdir,"Dist_dmpo_os-ORI_LHI_LH$ext"),f),figfmt)

plt = data(cfrf['A'])*mapping(:dmpo_os,color=:spiketype,dodge=:spiketype)*AlgebraOfGraphics.histogram(bins=-90:15:90)
f=draw(plt)
foreach(ext->save(joinpath(orisfdir,"Dist_dmpo_os-spiketype$ext"),f),figfmt)

plt = data(filter(r->r.ori_lhi_LH=="L",cfrf['A']))*mapping(:dmpo_os,color=:spiketype,dodge=:spiketype)*AlgebraOfGraphics.histogram(bins=-90:15:90)
f=draw(plt)
foreach(ext->save(joinpath(orisfdir,"Dist_dmpo_os-spiketype_ORI_LHI_L$ext"),f),figfmt)


plt = data(cfrf['A'])*mapping(:dmpd_ds,color=:ori_lhi_LH,dodge=:ori_lhi_LH)*AlgebraOfGraphics.histogram(bins=-180:15:180)
f=draw(plt)
foreach(ext->save(joinpath(orisfdir,"Dist_dmpd_ds-ORI_LHI_LH$ext"),f),figfmt)

plt = data(cfrf['A'])*mapping(:dmpd_ds,color=:spiketype,dodge=:spiketype)*AlgebraOfGraphics.histogram(bins=-180:15:180)
f=draw(plt)
foreach(ext->save(joinpath(orisfdir,"Dist_dmpd_ds-spiketype$ext"),f),figfmt)

plt = data(filter(r->r.ori_lhi_LH=="H",cfrf['A']))*mapping(:dmpd_ds,color=:spiketype,dodge=:spiketype)*AlgebraOfGraphics.histogram(bins=-180:15:180)
f=draw(plt)
foreach(ext->save(joinpath(orisfdir,"Dist_dmpd_ds-spiketype_ORI_LHI_H$ext"),f),figfmt)


pldf=combine(groupby(cfrf['A'],[:siteid,:layer]),[:po,:ou]=>((o,u)->begin
    osi = .!u
    lmpo = orim(o)
    lmpo_os = count(osi)<3 ? NaN : orim(o[osi])
    lcvpo = oricv(o)
    lcvpo_os = count(osi)<3 ? NaN : oricv(o[osi])
    (;lmpo,lmpo_os,lcvpo,lcvpo_os)
    end)=>AsTable,[:psf,:sfu]=>((f,u)->begin
    fsi = .!u
    lmpsf = mean(f)
    lmpsf_fs = count(fsi)<3 ? NaN : mean(f[fsi])
    lsdpsf = std(f)
    lsdpsf_fs = count(fsi)<3 ? NaN : std(f[fsi])
    (;lmpsf,lmpsf_fs,lsdpsf,lsdpsf_fs)
    end)=>AsTable,[:ocv,:hw,:osi,:dsi,:dcv].=>mean.=>x->"l$x",[:ocv_os,:hw_os,:osi_os,:dcv_ds,:dsi_ds].=>(x->begin
        xx = filter(!isnan,x)
        length(xx)<3 ? NaN : mean(xx)
    end).=>x->"l$x")
leftjoin!(pldf,pdf,on=:siteid)
transform!(pldf,[:lmpo,:mpo]=>orid=>:dlmpo,[:lmpo_os,:mpo_os]=>orid=>:dlmpo_os)
transform!(pldf,[:dlmpo,:dlmpo_os].=>ByRow(abs).=>x->"a$x")

plsdf=combine(groupby(cfrf['A'],[:siteid,:layer,:spiketype]),[:po,:ou]=>((o,u)->begin
    osi = .!u
    lmpo = orim(o)
    lmpo_os = count(osi)<3 ? NaN : orim(o[osi])
    lcvpo = oricv(o)
    lcvpo_os = count(osi)<3 ? NaN : oricv(o[osi])
    (;lmpo,lmpo_os,lcvpo,lcvpo_os)
    end)=>AsTable,[:ocv,:hw,:osi,:dsi,:dcv].=>mean.=>x->"l$x",[:ocv_os,:hw_os,:osi_os].=>(x->begin
        xx = filter(!isnan,x)
        length(xx)<3 ? NaN : mean(xx)
    end).=>x->"l$x",
    :psf=>sfm=>:lpsf)
leftjoin!(plsdf,pdf,on=:siteid)
transform!(plsdf,[:lmpo,:mpo]=>orid=>:dlmpo,[:lmpo_os,:mpo_os]=>orid=>:dlmpo_os)
transform!(plsdf,[:dlmpo,:dlmpo_os].=>ByRow(abs).=>x->"a$x")
plsdf = transform(groupby(plsdf,[:layer,:spiketype]),:lmpo_os=>(x->length(filter(!isnan,x))<4 ? missing : 1)=>:vn)
dropmissing!(plsdf,:vn)







plotlayerdist(cfrf['A'],x=:lr,g=:cofd_NAC,cc=dcc)
plotlayerdist(cfrf['A'],x=:alr,g=:cofd_NAC,cc=dcc)
# plotlayerdist(cfrf['A'],x=:dmpo_os,g=:ori_lhi_LH,cc=lhcc,dir=orisfdir,ts='A')
# plotlayerdist(cfrf['A'],x=:admpo_os,g=:ori_lhi_LH,cc=lhcc,dir=orisfdir,ts='A')

plotlayerdist(cfrf['A'],x=:dmpd_ds,g=:ori_lhi_LH,cc=lhcc,plottype=:violin)
plotlayerdist(cfrf['A'],x=:admpd_ds,g=:ori_lhi_LH,cc=lhcc,plottype=:violin)
# plotlayerdist(cfrf['A'],x=:dmpd_ds,g=:ori_lhi_LH,cc=lhcc,dir=orisfdir,ts='A',plottype=:violin)
# plotlayerdist(cfrf['A'],x=:admpd_ds,g=:ori_lhi_LH,cc=lhcc,dir=orisfdir,ts='A',plottype=:violin)

plotlayerdist(pldf,x=:dlmpo_os,g=:ori_lhi_LH,cc=lhcc,plottype=:box)
plotlayerdist(pldf,x=:adlmpo_os,g=:ori_lhi_LH,cc=lhcc,plottype=:box)
plotlayerdist(pldf,x=:lcvpo_os,g=:ori_lhi_LH,cc=lhcc,plottype=:box)
# plotlayerdist(pldf,x=:dlmpo_os,g=:ori_lhi_LH,cc=lhcc,dir=orisfdir,ts='A')
# plotlayerdist(pldf,x=:adlmpo_os,g=:ori_lhi_LH,cc=lhcc,dir=orisfdir,ts='A')
# plotlayerdist(pldf,x=:lcvpo_os,g=:ori_lhi_LH,cc=lhcc,dir=orisfdir,ts='A')



t = filter(r->!startswith(r.siteid,"AG2"),cfrf['Y'])
plotlayerscatter(t,y=:psi,x=:cofdi_L,ts='A')




plotlayerscatter(pldf,y=:dlmpo_os,x=:ori_lhi,ts='A')
plotlayerscatter(pldf,y=:adlmpo_os,x=:ori_lhi,ts='A')
plotlayerscatter(pldf,y=:lcvpo_os,x=:ori_lhi,ts='A')
# plotlayerscatter(pldf,y=:dlmpo_os,x=:ori_lhi,dir=orisfdir,ts='A')
# plotlayerscatter(pldf,y=:adlmpo_os,x=:ori_lhi,dir=orisfdir,ts='A')
# plotlayerscatter(pldf,y=:lcvpo_os,x=:ori_lhi,dir=orisfdir,ts='A')


plotlayerscatter(cfrf['A'],y=:dcv_ds,x=:ori_lhi,ts='A')
plotlayerscatter(cfrf['A'],y=:ocv_os,x=:ori_lhi,ts='A',g=:spiketype,msca=0)
# plotlayerscatter(cfrf['A'],y=:dcv_ds,x=:ori_lhi,dir=orisfdir,ts='A')
# plotlayerscatter(cfrf['A'],y=:dcv_ds,x=:ori_lr,dir=orisfdir,ts='A',g=:spiketype,msca=0)

plotlayerscatter(pldf,y=:losi_os,x=:ori_lhi,ts='A')
plotlayerscatter(pldf,y=:locv,x=:ori_lhi,ts='A',g=:spiketype,msca=0)
plotlayerscatter(plsdf,y=:lhw_os,x=:ori_lhi,ts='A',g=:spiketype,msca=0)
plotlayerscatter(pldf,y=:ldcv_ds,x=:ori_lr,dir=orisfdir,ts='A')
plotlayerscatter(pldf,y=:ldsi_ds,x=:ori_lhi,dir=orisfdir,ts='A',g=:spiketype,msca=0)
plotlayerscatter(plsdf,y=:losi_os,x=:ori_lr,dir=orisfdir,ts='A',g=:spiketype,msca=0)


plotlayerdist(cfrf['A'],x=:dsi,g=:spiketype,cc=scc)
# plotlayerdist(cfrf['A'],x=:dsi_ds,g=:spiketype,cc=scc,dir=orisfdir,ts='A')
plotlayerdist(plsdf,x=:losi_os,g=:spiketype,cc=scc)
# plotlayerdist(plsdf,x=:lhw_os,g=:spiketype,cc=scc,dir=orisfdir,ts='A')

plotlayerdist(cfrf['A'],x=:ocv,g=:cofd_NAC,cc=dcc)
plotlayerdist(cfrf['A'],x=:dcv_ds,g=:cofd_NAC,cc=dcc,dir=orisfdir,ts='A')
plotlayerdist(pldf,x=:losi_os,g=:cofd_NAC,cc=dcc,plottype=:density)
plotlayerdist(pldf,x=:ldsi_ds,g=:cofd_NAC,cc=dcc,dir=orisfdir,ts='A')


plotlayerscatter(cfrf['A'],y=:psf_fs,x=:ori_lhi,ts='A')
# plotlayerscatter(cfrf['A'],y=:psf_fs,x=:ori_lr,ts='A',dir=orisfdir)
plotlayerscatter(pldf,y=:lmpsf_fs,x=:ori_lhi,ts='A')
# plotlayerscatter(pldf,y=:lmpsf_fs,x=:ori_lhi,ts='A',dir=orisfdir)


plotlayerdist(cfrf['A'],x=:psf_fs,g=:cofd_NAC,cc=dcc)
plotlayerdist(cfrf['A'],x=:psf,g=:cofd_NAC,cc=dcc,ts='A',dir=orisfdir)
plotlayerdist(pldf,x=:lmpsf_fs,g=:cofd_NAC,cc=dcc)
plotlayerdist(pldf,x=:lmpsf,g=:cofd_NAC,cc=dcc,ts='A',dir=orisfdir)
plotlayerdist(pldf,x=:lsdpsf_fs,g=:cofd_NAC,cc=dcc,ts='A',dir=orisfdir)

