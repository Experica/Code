using NeuroAnalysis,FileIO,JLD2,DataFrames,StatsPlots,StatsBase,ProgressMeter,VegaLite,Dierckx,Combinatorics,XLSX

ccode = Dict("DKL_X"=>"A","DKL_Y"=>"Y","DKL_Z"=>"S","LMS_Xmcc"=>"L","LMS_Ymcc"=>"M","LMS_Zmcc"=>"S",
             "LMS_X"=>"L","LMS_Y"=>"M","LMS_Z"=>"S","DKL_Hue_L0"=>"DKL_L0","HSL_Hue_Ym"=>"HSL_Ym",
             "DKL_HueL0"=>"DKL_L0","HSL_HueYm"=>"HSL_Ym","DKL_L0"=>"DKL","HSL_Ym"=>"HSL")

getccode(k,ccode) = haskey(ccode,k) ? ccode[k] : k

function collectcondtest(indir;id="OriSF",datafile="factorresponse.jld2")
    rs=[]
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            exenv = load(joinpath(root,datafile),"exenv")
            haskey(exenv,"ID") || continue
            exenv["ID"] == id || continue
            push!(rs,load(joinpath(root,datafile)))
        end
    end
    return rs
end

function joinorisf(rs;ccode=ccode)
    siteid = rs[1]["siteid"]
    fa = rs[1]["fa"]

    dataset = Dict("siteid"=>siteid,"fa"=>fa)
    colors = [r["exenv"]["color"] for r in rs]
    colorcodes = map(i->getccode(i,ccode),colors)
    corder = sortperm(colorcodes)
    dataset["color"] = colors[corder]
    dataset["ccode"] = colorcodes[corder]
    rs = rs[corder]
    dataset["minmaxcolor"] = map(r->r["exenv"]["minmaxcolor"],rs)
    tfs = map(r->r["exenv"]["TemporalFreq"],rs)
    if length(unique(tfs))==1
        tfs = tfs[1]
    else
        @warn "Different TemporalFreq: $tfs for Blocks of OriSF."
    end
    dataset["tf"]=tfs
    gts = map(r->r["exenv"]["GratingType"],rs)
    if length(unique(gts))==1
        gts = gts[1]
    else
        @warn "Different GratingType: $gts for Blocks of OriSF."
    end
    dataset["GratingType"]=gts
    dataset["conddur"]=1500

    uns = map(r->length(r["unitid"]),rs)
    @info "Number of units for each OriSF test: $uns"
    uids = mapreduce(r->r["unitid"],intersect,rs) # only include units that spikes in all tests
    uis = map(r->indexin(uids,r["unitid"]),rs)
    dataset["ugood"] = Dict(uids.=>rs[1]["unitgood"][uis[1]])

    getru = (n;k=nothing) -> isnothing(k) ? map((r,i)->r[n][i],rs,uis) : map((r,i)->r[n][k][i],rs,uis)
    dataset["frs"] = Dict(uids.=>zip(getru("frs")...))
    dataset["responsive"] = Dict(uids.=>zip(getru("responsive")...))
    dataset["modulative"] = Dict(uids.=>zip(getru("modulative")...))
    dataset["enoughresponse"] = Dict(uids.=>zip(getru("enoughresponse")...))
    dataset["maxi"] = Dict(uids.=>zip(getru("maxi")...))
    dataset["fpsth"] = (psth=Dict(uids.=>zip(getru("fpsth";k=:psth)...)), x = rs[1]["fpsth"].x)
    dataset["frf"] = NamedTuple(k=>Dict(uids.=>zip(getru("frf";k)...)) for k in keys(fa))
    
    return dataset
end


function f1f0orisf!(dataset;btw=(-50,0))
    conddur = dataset["conddur"]
    tf = dataset["tf"]
    fpsth,x = dataset["fpsth"]
    bw = x[2]-x[1]
    fs = 1/(bw/1000)
    bi = findall(i->btw[begin]<=i<btw[end],x)
    ci = findall(i->0<=i<conddur,x)
    getF1F0 = (psth)->begin
        y = psth[ci].-mean(psth[bi])
        F1F0 = dft(y,fs,tf,0)
        m0 = abs(F1F0[2])
        m1 = abs(F1F0[1])
        a1 = mod2pi(angle(F1F0[1]))
        (;m0,m1,a1,f10=m1/m0)
    end
    f1f0 = Dict(u=>map(c->getF1F0.(c.m),v) for (u,v) in fpsth)
    f10 = Dict(u=>map(c->getproperty.(c,:f10),v) for (u,v) in f1f0)
    f1phase = Dict(u=>map(c->getproperty.(c,:a1),v) for (u,v) in f1f0)
    f1maxi = Dict(u=>map(c->[Tuple(argmax(getproperty.(c,:m1)))...],v) for (u,v) in f1f0)

    dataset["f1f0"] = f1f0
    dataset["f10"] = f10
    dataset["f1phase"] = f1phase
    dataset["f1maxi"] = f1maxi
    return dataset
end

plotpsth = (dataset,u;isse=true,showfactor=true,dir=nothing,figfmt=[".png"],cs=["A","L","M","S"])->begin
    ccode = dataset["ccode"]
    cs = intersect(cs,ccode)
    ci = indexin(cs,ccode)
    cn = length(ci)

    conddur = dataset["conddur"]
    fa = dataset["fa"]
    fan1,fan2 = length.(values(fa))
    fpsth,x = dataset["fpsth"]
    psth = fpsth[u][ci]
    maxi = dataset["maxi"][u][ci]
    f1maxi = dataset["f1maxi"][u][ci]
    f1phase = dataset["f1phase"][u][ci]
    ug = dataset["ugood"][u]

    colors = map(c->c.maxcolor,dataset["minmaxcolor"][ci])
    "A" in cs && (colors[findfirst("A".==cs)]=RGBA(0.3,0.3,0.3,1))

    p = plot(;layout=(fan1,fan2),grid=false,legend=false,size=(200fan2,100fan1),
            link=:all,tickdir=:out,xlabel="Time (ms)",ylabel="Response (spike/s)")
    for j in 1:fan1, i in 1:fan2
        leftmargin = i == 1 ? 3Plots.mm : -18Plots.mm
        bottommargin = j==fan1 ? 3Plots.mm : -16Plots.mm
        topmargin = j==1 ? 6Plots.mm : :match
        rightmargin = i==fan2 ? 8Plots.mm : :match
        frame = (i ==1 && j==fan1) ? :auto : :none
        
        vspan!(p[j,i],[0,conddur];color=:gray92,leftmargin,bottommargin,topmargin,rightmargin,frame)
        for c in 1:cn
            ribbon = isse ? psth[c].se[j,i] : nothing
            y = psth[c].m[j,i]
            ann = ((0.1,0.9),("-",20,colors[c],:left,:vcenter,rad2deg(f1phase[c][j,i])))
            plot!(p[j,i],x,y;ribbon,color=colors[c],ann)
        end
    end
    
    ax = 0.06.*(0:cn-1)
    ax = ax .- mean(ax) .+ 0.5
    foreach((i,c,x)->annotate!(p[i...],((x,0.96),("∘",16,c,:center))),maxi,colors,ax)
    foreach((i,c,x)->annotate!(p[i...],((x,0.99),("▵",16,c,:center))),f1maxi,colors,ax)
    if showfactor
        foreach(i->annotate!(p[i,fan2],((1.1,0.5),("-ꜛ-",20,:center,fa[1][i]))),eachindex(fa[1]))
        foreach(i->annotate!(p[1,i],((0.5,1.2),("SF=$(round(fa[2][i],digits=1))",12,:center))),eachindex(fa[2]))
    end
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$(ug ? "S" : "M")U$(u)_PSTH$ext")),figfmt)
end

orisfinfo = (siteresultdir;figfmt=[".png"]) -> begin
    rs = collectcondtest(siteresultdir,id="OriSF")
    dataset = joinorisf(rs)
    dataset = f1f0orisf!(dataset)
    jldsave(joinpath(siteresultdir,"orisfdataset.jld2");dataset)
end

resultroot = "Z:/"


plot(rand(10),ann=((0,0),text("-",50,:left,:vcenter,-0.0)))



## process all orisfs of a RecordSite
subject = "AG1";recordsession = "V1";recordsite = "ODL1"
siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
siteresultdir = joinpath(resultroot,subject,siteid)
# dataset = load(joinpath(siteresultdir,"orisfdataset.jld2"),"dataset")


using Interact
@manipulate for u in sort(collect(keys(dataset["ugood"])))
    plotlstas(dataset,u)
end

plotpsth(dataset,748,isse=false,cs=["A","L","M","S"])
plotpsth(dataset,748,isse=false,cs=["L","M"])

## Batch Penetration Sites
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
@showprogress "Batch All OriSFs ... " for r in eachrow(penetration)
    orisfinfo(joinpath(resultroot,r.Subject_ID,r.siteid);figfmt=[".png",".svg"])
end

pythonplot()
pyplot()






























allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
layer = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate")
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1")...)
penetration = select(penetration,[:siteid,:od,:cofd],:cofd=>ByRow(i->i∈["L/M","S/LM","L/M, S/LM"])=>:incofd,
                :cofd=>ByRow(i->i∈["B","W","B, W"])=>:inbw,:cofd=>ByRow(i->i=="None")=>:none)




## Collect All Color Tests
colortest = collectcondtest(resultroot,exid="Color")

function colorunit(unit;color="DKL_L0",ccode=ccode)
    df = dropmissing(unit,"responsive!$color")
    df.c .= getccode(color,ccode)
    select!(df,Not(r"\w*!\w*"),
        ["responsive!$color","enoughresponse!$color"] => ByRow(&) => :er,
        "fa_Angle!$color" => :fl,
        ["fms!$color","maxfri_Angle!$color"]=>ByRow((fm,fi)->fm[fi...])=>:fr,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : (;i.fit.mfit.model,i.fit.mfit.fun,i.fit.mfit.param)) => :fit,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.up) => :up,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.cm) => :cm,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.cv) => :cv,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.aup) => :aup,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.acm) => :acm,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.acv) => :acv,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.fit.pa) => :pa,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.fit.asi2) => :asi,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : sum(i.fit.ahw)) => :ahw,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : i.fit.mfit.r) => :cor,
        "frf_Angle!$color" => ByRow(i->ismissing(i) ? missing : 1-i.fit.mfit.r2) => :fvu)
end

# allcolorunit = mapreduce(c->colorunit(colortest,color=c),(i,j)->append!(i,j,cols=:union),["DKL_L0","HSL_Ym"])
# transform!(allcolorunit,:up=>(i->i.<0.05)=>:s,:aup=>(i->i.<0.05)=>:as)
# jldsave(joinpath(resultroot,"allcolorunit.jld2");allcolorunit)
allcolorunit = load(joinpath(resultroot,"allcolorunit.jld2"),"allcolorunit")
dklcg = cgrad(ColorMaps["lidkl_mcchue_l0"].colors)
hslcg = cgrad(ColorMaps["hsl_mshue_l0.4"].colors)
# dklcm = ColorMaps["lidkl_mcchue_l0"].colors
# hslcm = ColorMaps["hsl_mshue_l0.4"].colors

## color responsive
colorcu = innerjoin(allcu,allcolorunit,on=[:siteid,:id])
leftjoin!(colorcu,penetration,on=:siteid)
colorcsu = subset(colorcu,:good)
dklcsu = filter(r->r.c=="DKL",colorcsu)
hslcsu = filter(r->r.c=="HSL",colorcsu)
dklrcsu = subset(dklcsu,:er)
hslrcsu = subset(hslcsu,:er)

allcolorcsu = combine(groupby(colorcsu,[:id,:siteid]),
            :aligndepth=>first=>identity,
            :layer=>first=>identity,
            :er=>any=>:eresponsive,
            [:c,:er]=>((c,r)->join(c[r]))=>:RTYPE,
            [:c,:er,:cm]=>((c,r,v)->OrderedDict((c[r].=>v[r])...))=>last)
allcolorrcsu = subset(allcolorcsu,:eresponsive)


plotdepthhist = (unit;g=:responsive,dir=nothing,figfmt=[".png"],layer=nothing,title="CorticalDepth_$g",ylabel="Number of Units",palette=:tab20,leg=:best) -> begin
p = groupedhist(unit.aligndepth;group=unit[!,g],barpositions=:stack,bin=0:0.02:1,permute=(:y,:x),xflip=true,grid=false,legendfontsize=6,
    xlims=(0,1),lw=0,size=(350,500),tickor=:out,xlabel="Cortical Depth",ylabel,palette,leg)
if !isnothing(layer)
    ann = [(1,mean(layer[k]),Plots.text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
    hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
end
isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end


plotdepthhist(allcolorcsu;g=:eresponsive,layer)
plotdepthhist(allcolorrcsu;g=:RTYPE,layer)


## color tuning
plotcolortuning = (unit;type=:data,cg=dklcg,dir=nothing,figfmt=[".png"],layer=nothing,title="ColorTuning") -> begin
x = 0:0.04:2π # 2.3 deg
if type == :fit
    tc = hcat(map(f->predict(f,x),unit.fit)...)
elseif type == :data
    tc = hcat(map((l,r)->Spline1D([deg2rad.(l);2π],[r;r[1]],k=1)(x),unit.fl,unit.fr)...)
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

# foreach(siteid->plotcolortuning(filter(r->r.siteid==siteid,dklrcsu);
#     cg=dklcg,dir=joinpath(resultroot,"Tuning"),type=:fit,title="$(siteid)_DKL_ColorTuning_fit"),
#     levels(dklrcsu.siteid))
# foreach(siteid->plotcolortuning(filter(r->r.siteid==siteid,hslrcsu);
#     cg=hslcg,dir=joinpath(resultroot,"Tuning"),type=:fit,title="$(siteid)_HSL_ColorTuning_fit"),
#     levels(hslrcsu.siteid))

plotcolortuning(filter(r->r.siteid=="AG2_V1_ODL4",dklrcsu),cg=dklcg,type=:fit)
plotcolortuning(filter(r->r.siteid=="AG2_V1_ODL4",hslrcsu),cg=hslcg,type=:fit)




# tall=unitdensity(dklcu.aligndepth,spacerange=(0,1),bw=0.02,step=0.02)
# ter = unitdensity(dklcu.aligndepth[dklcu.er],spacerange=(0,1),bw=0.02,step=0.02)
# plot(tall.n,tall.y,yflip=true,size=(350,500))
# plot!(ter.n,ter.y,yflip=true)

# plot(ter.n./tall.n*100,tall.y,yflip=true,size=(350,500))

# bar(tall.y,ter.n./tall.n*100,permute=(:y,:x),xflip=true,size=(350,500),lw=0,barwidths=0.02)
# ann = [(0,mean(layer[k]),text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
# hline!([l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)







# cunit = filter(r->r.er && r.siteid=="AG2_V1_ODR6",dklcu)
# cunit = filter(r->r.er,dklcu)
# cunit = filter(r->r.er,hslcu)
# x = 0:0.04:2π # 2.3 deg
# tc = hcat(map(f->predict(f,x),cunit.fit)...)
# tc./=maximum(tc,dims=1)
# d =repeat(cunit.aligndepth',inner=(length(x),1))

# plot(x,tc,d,leg=false,color=dklcg,lz=x)

# xx = cos.(x).*tc
# yy = sin.(x).*tc

# plot(xx,yy,d,leg=false,frame=:origin,xlims=(-1,1),ylims=(-1,1),zlims=(0,1),color=dklcg,lz=x,zflip=true)



# y = 0:0.005:1
# is = [findall(y[i] .<= cunit.aligndepth .< y[i+1]) for i in 1:length(y)-1]

# tt = [mean(tc[i,j]) for j in is, i in eachindex(x)]
# replace!(tt,NaN=>0)
# tt=clampscale(tt)

# heatmap(tt,yflip=true,size=(350,500))

# ccg = dklcg[range(0,1,length=158)]
# ccg = hslcg[range(0,1,length=158)]
# dmc = coloralpha.(color.(ccg'),tt);

# plot(dmc)





plotdepthscatter = (unit;x=:cv,dir=nothing,siteid="all",figfmt=[".png"],layer=nothing,xlabel="",xticks=:auto,xlims=:auto,
    title="$(siteid)_CorticalDepth",wfun=x->isempty(x) ? 0 : median(x),color=:auto,palette=:tab20,leg=:best) -> begin
cunit = siteid=="all" ? unit : filter(r->r.siteid==siteid,unit)
x = cunit[!,x]
xmin,xmax = xlims==:auto ? extrema(x) : xlims
p = scatter(x,cunit.aligndepth;yflip=true,legendfontsize=6,xticks,xlims,color,
    ms=3,msw=0,ma=0.7,ylims=(0,1),size=(350,500),tickor=:out,ylabel="Cortical Depth",xlabel,palette,leg)
if !isnothing(wfun)
    n,y = unitdensity(cunit.aligndepth;w=x,wfun,spacerange=(0,1),bw=0.05,step=0.05)
    plot!(p,n,y,color=:gray40,label="Average")
end
if !isnothing(layer)
    ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
    hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
end
isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end


plotdepthscatter(hslrcsu;x=:cv,layer,xlabel="cv",leg=false)
plotdepthscatter(hslrcsu;x=:cm,layer,xlabel="cm",leg=false,xticks=0:90:360)

p=select(hslrcsu,[:cm,:cv,:cor,:fvu,:pa,:ahw,:layer,:siteid]) |>
[@vlplot(:bar,y={"siteid"},x={"count()"});
@vlplot(:tick,y={"layer"},x={"cv"});
@vlplot(:tick,y={"layer"},x={"cm",axis={values=0:90:360}});
@vlplot(:bar,y={"count()"},x={"cor",bin={step=0.05}});
@vlplot(:bar,y={"count()"},x={"fvu",bin={step=0.05}});
@vlplot(:tick,y={"layer"},x={"ahw"});
@vlplot(:tick,y={"layer"},x={"pa",axis={values=0:90:360}})]

plotdepthscatter(dklrcsu;siteid="AG2_V1_ODL4",x=:cm,layer,wfun=nothing,xlabel="cm",leg=false,xticks=0:90:360,xlims=(-5,365))

p=subset(select(hslrcsu,[:s,:cm,:cv,:pa,:layer,:siteid,:aligndepth]),:s) |>
@vlplot(:point,y={"aligndepth",sort="descending"},x={"cm",axis={values=0:90:360}},color={"layer",scale={scheme=:category10}},
    columns=9,wrap={"siteid"})





# filter(r->r.rme,dklcell) |> @vlplot(:bar,transform=[{calculate="datum.L_M > 0 ? 'L+' : 'L-'","as"="Dominance"},{calculate="datum.L_M > 0 ? 1 : -1","as"="sign"}],
#         y={"layer",title="Cortical Layer"},x={"sign",title="Number of rme Cells"},color={"Dominance",scale={range="#" .* hex.(get(dklcg,[0,0.5]),:rgb)}},column={"site",title=""})
# filter(r->r.rme,dklcell) |> @vlplot(:bar,transform=[{calculate="datum.S_LM > 0 ? 'S+' : 'S-'","as"="Dominance"},{calculate="datum.S_LM > 0 ? 1 : -1","as"="sign"}],
#         y={"layer",title="Cortical Layer"},x={"sign",title="Number of rme Cells"},color={"Dominance",scale={range="#" .* hex.(get(dklcg,[0.25,0.75]),:rgb)}},column={"site",title=""})



# dklcell |> @vlplot(:bar,transform=[{bin={step=30},field="oh","as"="boh"}],
#             y={"count()",title="Number of Cells"},x={"oh",bin={step=30},title="Optimal Hue",axis={values=0:90:360}},
#             color={"boh",scale={range="#" .* hex.(get(dklcg,range(1/24,length=12,step=1/12)),:rgb)},legend=false},column={"site",title=""})

# hslcell |> @vlplot(:bar,transform=[{bin={step=30},field="oh","as"="boh"}],
#             y={"count()",title="Number of Cells"},x={"oh",bin={step=30},title="Optimal Hue",axis={values=0:90:360}},
#             color={"boh",scale={range="#" .* hex.(get(hslcg,range(1/24,length=12,step=1/12)),:rgb)},legend=false},column={"site",title=""})



# dklcs = "#" .* hex.(dklcg.colors.colors,:rgb)
# hslcs = "#" .* hex.(hslcg.colors.colors,:rgb)


# dklcell |> @vlplot(:point,y={"depth",sort="descending",title="Cortical Depth"},x={"hcv",title="Number of Cells"},
#                         color={"oh",scale={range=dklcs}},column={"site"})

# dklcells |> @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"dkl_oh",scale={range=dklcs}},column={"site"})

# hslcell |> @vlplot(:bar,y={"depth",bin={maxbins=20},sort="descending",title="Cortical Depth"},x={"count()",title="Number of Cells"},
#                         color={"oh",scale={range=hslcs}},column={"site"})

# hslcell |> @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"oh",scale={range=hslcs}},column={"site"})






## Collect All OriSF Tests
orisftest = collectcondtest(resultroot,exid="OriSF")

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


