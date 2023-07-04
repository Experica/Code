using NeuroAnalysis,DataFrames,JLD2,StatsBase,StatsPlots,LinearAlgebra,Images,VegaLite,XLSX,DSP,
    UMAP,Clustering,Distances,Combinatorics,DataStructures
import GLMakie,Contour


function collectsta(indir;unit=DataFrame(),datafile="stadataset.jld2")
    for (root,dirs,files) in walkdir(indir)
        if datafile in files
            dataset = load(joinpath(root,datafile),"dataset")
            siteid = dataset["siteid"]
            delays = dataset["delays"]
            sizedeg = dataset["sizedeg"]
            unitid = collect(keys(dataset["usta"]))
            unitgood = map(u->dataset["ugood"][u],unitid)
            unitresponsive = map(u->dataset["uresponsive"][u],unitid)
            id = map((u,g)->"$(siteid)_$(g ? "S" : "M")U$u",unitid,unitgood)
            roi = map((u,r)->r ? (;dataset["ulroi"][u].centerdeg,dataset["ulroi"][u].radiideg,dataset["ulroi"][u].radiusdeg) : missing,unitid,unitresponsive)
            rc = map((u,r)->r ? map((c,cr)->cr ? c : missing,dataset["ccode"],dataset["ucresponsive"][u]) : missing,unitid,unitresponsive)
            dc = map((u,r)->r ? map((d,cr)->cr ? delays[d.d] : missing,dataset["ulcrpd"][u],dataset["ucresponsive"][u]) : missing,unitid,unitresponsive)
            df = DataFrame(siteid=siteid,id=id,responsive=unitresponsive,roi=roi,rc=rc,dc=dc,sizedeg=sizedeg)
            if haskey(dataset,"ulfit")
                foreach(m->df[!,"sta!$m"]=map((u,r)->r ? map(f->ismissing(f) ? missing : (;(k=>f[k] for k in keys(f))...),dataset["ulfit"][u][m]) : missing,
                        unitid,unitresponsive), keys(first(values(dataset["ulfit"]))))
            end
            append!(unit,df,cols=:union)
        end
    end
    return unit
end

resultroot = "Z:/"
# staunit = collectsta(resultroot)
# jldsave(joinpath(resultroot,"staunit.jld2");staunit)



## cortical unit STA responsive
staunit = load(joinpath(resultroot,"staunit.jld2"),"staunit")
layer,nlbt = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate","nlbt")
layercg = cgrad(:tab10,nlbt,categorical=true)
allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
allcsu = subset(allcu,:good)
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))

stacu = innerjoin(allcu,staunit,on=[:siteid,:id])
leftjoin!(stacu,penetration[:,[:siteid,:od,:cofd,:pid]],on=:siteid)
transform!(stacu,:rc=>ByRow(i->ismissing(i) ? missing : join(skipmissing(i)))=>:RTYPE)
stacu.CTYPE = map(i->ismissing(i) ? missing : replace(i,r"A([L,M,S]+)"=>s"\1"),stacu.RTYPE)
stacsu = subset(stacu,:good)
starcsu = dropmissing(stacsu,:rc)
stanrcsu = subset(stacsu,:responsive=>.!)
jldsave(joinpath(resultroot,"stacsu.jld2");stacsu)
rfdir = joinpath(resultroot,"RF");mkpath(rfdir)
figfmt = [".svg",".png"]

plotlstaroi = (unit;siteid="AG1_V1_ODL1",cg=layercg,dir=nothing,figfmt=[".png"],title="$(siteid)_lsta_roi")->begin
    sunit = sort!(filter(r->r.siteid==siteid,unit),:aligndepth,rev=true) # from bottom to top of cortex
    sizedeg = sunit.sizedeg[1] # stimulus diameters
    # standard 2D coordinates(rightward:x, upward:y)
    stix = [0,1,1,0,0] * sizedeg[2]
    stiy = [1,1,0,0,1] * sizedeg[1]
    # roi is defined on pixel indices of stimulus image
    cx = map(i->i.centerdeg[2],sunit.roi)
    cy = map(i->sizedeg[1]-i.centerdeg[1],sunit.roi)
    rx = map(i->i.radiideg[2],sunit.roi)
    ry = map(i->i.radiideg[1],sunit.roi)
    roix = [-1,1,1,-1,-1] .* rx' .+ cx'
    roiy = [1,1,-1,-1,1] .* ry' .+ cy'
    p = plot(stix,stiy,lw=2,size=(sizedeg[2],sizedeg[1]).*200,color=:black,frame=:grid,tickdir=:out,ratio=1,
            xlabel="X (deg)",ylabel="Y (deg)",leg=false)
    plot!(p,roix,roiy;la=0.9,lw=1,color=cg[sunit.aligndepth]')
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotlstaroi(starcsu;siteid="AG2_V1_ODR13")
foreach(siteid->plotlstaroi(starcsu;siteid,dir=rfdir,figfmt),levels(starcsu.siteid))


plotdepthhist = (unit;g=:responsive,dir=nothing,figfmt=[".png"],layer=nothing,title="DepthHist_$g",ylabel="Number of Units",
                siteid=nothing,xlabel="Normalized Cortical Depth",size=(350,550),palette=:tab20,leg=:best) -> begin
    if isnothing(siteid)
        sunit=unit
    else
        sunit = filter(r->r.siteid==siteid,unit)
        title = siteid * '_' * title
    end
    p = groupedhist(sunit.aligndepth;group=sunit[!,g],barpositions=:stack,bin=0:0.02:1,permute=(:y,:x),xflip=true,grid=false,legendfontsize=6,
        xlims=(-0.02,1.02),xticks=0:0.2:1,lw=0,size,tickor=:out,xlabel,ylabel,palette,leg,title,titlefontsize=10)
    if !isnothing(layer)
        ann = [(0.2,mean(layer[k]),Plots.text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
    end
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotdepthhist(allcsu;g=:spiketype,layer)
plotdepthhist(stacsu;g=:responsive,layer)
plotdepthhist(starcsu;g=:RTYPE,layer)
plotdepthhist(starcsu;g=:CTYPE,layer)

# plotdepthhist(allcsu;g=:spiketype,layer,dir=resultroot,figfmt)
# plotdepthhist(stacsu;g=:responsive,layer,dir=rfdir,figfmt)
# plotdepthhist(starcsu;g=:RTYPE,layer,dir=rfdir,figfmt)
# plotdepthhist(starcsu;g=:CTYPE,layer,dir=rfdir,figfmt)

# foreach(siteid->plotdepthhist(stacsu;dir=rfdir,figfmt,layer,g=:responsive,siteid),levels(stacsu.siteid))
# foreach(siteid->plotdepthhist(starcsu;dir=rfdir,figfmt,layer,g=:RTYPE,siteid),levels(starcsu.siteid))
# foreach(siteid->plotdepthhist(starcsu;dir=rfdir,figfmt,layer,g=:CTYPE,siteid),levels(starcsu.siteid))

# pl = select(starcsu,[:RTYPE,:siteid,:layer]) |>
#     [@vlplot(:bar,y={"layer"},x={"count()",axis={title="Number of Units",grid=false}},color={"RTYPE:n",scale={scheme=:category20}});
#     @vlplot(:bar,y={"siteid"},x={"count()",axis={title="Number of Units",grid=false}},color={"RTYPE:n",scale={scheme=:category20}})]
# foreach(ext->save(joinpath(rfdir,"Layer_SiteID_Hist_RTYPE$ext"),pl),figfmt)











# COFD, Achromatic and None
# penetration = select(penetration,[:siteid,:od,:cofd,:pid],:cofd=>ByRow(i->i∈["L/M","S/LM","L/M, S/LM"])=>:incofd,
#                 :cofd=>ByRow(i->i∈["B","W","B, W"])=>:inbw,:cofd=>ByRow(i->i=="None")=>:none)

# dn = :none
# plotdepthhist(subset(stacsu,dn);g=:responsive,layer)
# plotdepthhist(subset(starcsu,dn);g=:RTYPE,layer)

# p=select(subset(starcsu,dn),[:RTYPE,:layer]) |>
# [@vlplot(:bar,x={"RTYPE"},y={"count()"});
# @vlplot(:bar,y={"layer"},x={"count()"},color={"RTYPE:n",scale={scheme=:category20}})]













## RF Model Spatial Pattern
getstaedog = (unit;model=:edog) -> begin
    select!(dropmissing!(flatten(dropmissing(unit,:rc),["rc","dc","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->(;i.model,i.fun,i.mfun,i.cfun,i.param,i.radii))=>:fit,
        "sta!$model"=>ByRow(i->i.r)=>:PCC,
        "sta!$model"=>ByRow(i->i.bic)=>:BIC,
        "sta!$model"=>ByRow(i->1-i.r2)=>:FVU,
        "sta!$model"=>ByRow(i->i.param[1])=>:aₑ,
        "sta!$model"=>ByRow(i->i.param[7])=>:aᵢ,
        "sta!$model"=>ByRow(i->i.param[2])=>:μx,
        "sta!$model"=>ByRow(i->i.param[4])=>:μy,
        "sta!$model"=>ByRow(i->i.param[3])=>:σₑₓ,
        "sta!$model"=>ByRow(i->i.param[5])=>:σʸₓ,
        "sta!$model"=>ByRow(i->i.param[6])=>:θ,
        "sta!$model"=>ByRow(i->i.param[8])=>:σⁱₑ)
end
getstagabor = (unit;model=:gabor) -> begin
    select!(dropmissing!(flatten(dropmissing(unit,:rc),["rc","dc","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->(;i.model,i.fun,i.mfun,i.cfun,i.param,i.radii))=>:fit,
        "sta!$model"=>ByRow(i->i.r)=>:PCC,
        "sta!$model"=>ByRow(i->i.bic)=>:BIC,
        "sta!$model"=>ByRow(i->1-i.r2)=>:FVU,
        "sta!$model"=>ByRow(i->i.param[1])=>:a,
        "sta!$model"=>ByRow(i->i.param[2])=>:μx,
        "sta!$model"=>ByRow(i->i.param[4])=>:μy,
        "sta!$model"=>ByRow(i->i.param[3])=>:σₓ,
        "sta!$model"=>ByRow(i->i.param[5])=>:σʸₓ,
        "sta!$model"=>ByRow(i->i.param[6])=>:θ,
        "sta!$model"=>ByRow(i->i.param[7])=>:f,
        "sta!$model"=>ByRow(i->i.param[8])=>:p)
end
"generate fit model image"
function fitimage(fit; ppd=45, ppr=65, tight=false, metric=false)
    if tight # image centered on and encompass model, instead of lsta roi used to fit model
        param = deepcopy(fit.param)
        if fit.model == :edog
            param[[2, 4]] .= 0 # move to image center
            σₑₓ = param[3]
            σʸₓ = param[5]
            σⁱₑ = param[8]
            radius = 3 * max(σₑₓ, σʸₓ * σₑₓ, σⁱₑ * σₑₓ, σⁱₑ * σʸₓ * σₑₓ) # 3σ radius
            fit = (; fit.model, fit.fun, fit.mfun, param, radii=(radius, radius))
        elseif fit.model == :gabor
            param[[2, 4]] .= 0 # move to image center
            σₓ = param[3]
            σʸₓ = param[5]
            radius = 3 * max(σₓ, σʸₓ * σₓ) # 3σ radius
            fit = (; fit.model, fit.fun, fit.mfun, param, radii=(radius, radius))
        end
        # ensure ppr pixels for model encompassing radius so that every masked model has the same relative pixel sampling rate,
        # and metric evaluation has the same level of accuracy.
        ppd = ppr/radius
    end
    y = -fit.radii[1]:1/ppd:fit.radii[1]
    x = -fit.radii[2]:1/ppd:fit.radii[2]
    img = predict(fit, x, y, xygrid=true, yflip=true) # standard 2D coordinates(rightward:x, upward:y)
    if metric
        # metric should be calculated from masked model region of image
        mask = predict((; fun=fit.mfun, fit.param), x, y, xygrid=true, yflip=true)
        return (; img, mmmetric(img[mask])...)
    else
        return img
    end
end
"metrics based on masked model region of image"
function mmmetric(mm)
    n = length(mm)
    on = mm .> 0
    aon = count(on) / n # On Area ratio
    V = sum(abs.(mm)) # Total Volumn(pixel area unit)
    von = sum(mm[on]) / V # On Volumn ratio
    return (; aon, von, V)
end
"generate fit model image on new origin"
function fitimage(fit, μx, μy, r; ppr=70, ppd=ppr / r)
    param = deepcopy(fit.param)
    # use new origin(μx,μy)
    param[2] = param[2] - μx
    param[4] = param[4] - μy
    fit = (; fit.model, fit.fun, fit.mfun, param, radii=(r, r))

    y = -fit.radii[1]:1/ppd:fit.radii[1]
    x = -fit.radii[2]:1/ppd:fit.radii[2]
    # standard 2D coordinates(rightward:x, upward:y)
    img = predict(fit, x, y, xygrid=true, yflip=true)
    mask = predict((; fun=fit.mfun, fit.param), x, y, xygrid=true, yflip=true)
    (; img, mask)
end


# Spatial Pattern Parameters
edogcsu = getstaedog(stacsu)
transform!(edogcsu,
    [:σₑₓ,:σʸₓ,:σⁱₑ]=>ByRow((σₑₓ,σʸₓ,σⁱₑ)->5*max(σₑₓ, σʸₓ * σₑₓ, σⁱₑ * σₑₓ, σⁱₑ * σʸₓ * σₑₓ))=>:Diameter,
    :σʸₓ=>ByRow(log2∘inv)=>:Elliptical,
    :σʸₓ=>ByRow(abs∘log2∘inv)=>:Round,
    :σⁱₑ=>ByRow(log2)=>:CenterSurround,
    [:aᵢ,:aₑ]=>ByRow(/)=>:aⁱₑ,
    [:aᵢ,:aₑ]=>ByRow(log2∘inv∘/)=>:OnOff,
    [:aᵢ,:aₑ]=>ByRow(sign∘log2∘inv∘/)=>:OnOffSign,
    :fit=>ByRow(i->fitimage(i,tight=true,metric=true))=>AsTable)

gaborcsu = getstagabor(stacsu)
transform!(gaborcsu,
    [:σₓ,:σʸₓ]=>ByRow((σₓ,σʸₓ)->5*max(σₓ, σʸₓ * σₓ))=>:Diameter,
    :σʸₓ=>ByRow(log2∘inv)=>:Elliptical,
    :σʸₓ=>ByRow(abs∘log2∘inv)=>:Round,
    [:σₓ,:σʸₓ,:f]=>ByRow((σₓ,σʸₓ,f)->5*σʸₓ*σₓ*f)=>:Cycle,
    :p=>ByRow(p->sin(2π*p))=>:OnOff,
    :p=>ByRow(p->(sign∘sin)(2π*p))=>:OnOffSign,
    :p=>ByRow(p->(abs∘sin)(2π*p))=>:OddEven,
    :fit=>ByRow(i->fitimage(i,tight=true,metric=true))=>AsTable)


@df edogcsu scatter(:σⁱₑ,:aⁱₑ,size=(600,450),ylim=(-0.1,4.1),yticks=0:4,xlim=(-0.1,4.1),xticks=0:4,
                    leg=false,ms=3,msw=0,ma=0.7,xlabel="σᵢ/σₑ",ylabel="aᵢ/aₑ")
foreach(ext->savefig(joinpath(rfdir,"csu_edog_σⁱₑ_aⁱₑ$ext")),figfmt)
@df edogcsu scatter(:CenterSurround,:OnOff,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-10,10),
                    xlabel="CenterSurround",ylabel="OnOff")
foreach(ext->savefig(joinpath(rfdir,"csu_edog_cs_onoff$ext")),figfmt)
@df edogcsu corrplot(cols([:Diameter,:CenterSurround,:aon,:von]),bin=50,leg=false,size=(850,650),grid=false)
foreach(ext->savefig(joinpath(rfdir,"csu_edog_SPP$ext")),figfmt)

@df gaborcsu scatter(:p,:Cycle,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-0.1,4.1),yticks=0:4,
                    xlabel="Phase",ylabel="Cycle")
foreach(ext->savefig(joinpath(rfdir,"csu_gabor_p_cyc$ext")),figfmt)
@df gaborcsu scatter(:OnOff,:Cycle,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-0.1,4.1),yticks=0:4,
                    xlabel="OnOff",ylabel="Cycle")
foreach(ext->savefig(joinpath(rfdir,"csu_gabor_onoff_cyc$ext")),figfmt)
@df gaborcsu corrplot(cols([:Diameter,:OnOff,:aon,:von]),bin=50,leg=false,size=(850,650),grid=false)
foreach(ext->savefig(joinpath(rfdir,"csu_gabor_SPP$ext")),figfmt)



## Spatial Pattern Clustering
cmcode = Dict('A'=>ColorMaps["dkl_mcclumiso"],'L'=>ColorMaps["lms_mccliso"],'M'=>ColorMaps["lms_mccmiso"],'S'=>ColorMaps["lms_mccsiso"])
ALMScp=palette([:gray,:lightcoral,:seagreen,:blueviolet])
maskimg(i;size=(51,51),masktype="Gaussian") = float32.(alphamask(imresize(i,size);masktype))
maskimg(i,c;size=(51,51),masktype="Gaussian") = float32.(alphamask(imresize(i,size);color=cmcode[c].colors,masktype))

edogimg = map(maskimg,edogcsu.img)
gaborimg = map(maskimg,gaborcsu.img)
edogimgc = map((i,c)->maskimg(i,c),edogcsu.img,edogcsu.rc)
gaborimgc = map((i,c)->maskimg(i,c),gaborcsu.img,gaborcsu.rc)


# eDoG
features = [:OnOffSign,:aon,:von]
Fz = stack(i->zscore(edogcsu[:,i]),features)
FzD = pairwise(Euclidean(),Fz,dims=1)

Fz2 = umap(Fz', 2, n_neighbors=45, min_dist=2, n_epochs=250,metric=Euclidean())
scatter(Fz2[2,:], Fz2[1,:],ms=2,msw=0,ma=0.8,frame=:none,leg=false,ratio=1,size=(600,600))

fap = GLMakie.scatter(Fz2[2,:], Fz2[1,:], marker=edogimg, markersize=30, figure=(;resolution=(1200,1200)))
fap.axis.aspect = GLMakie.DataAspect()
GLMakie.hidespines!(fap.axis)
GLMakie.hidedecorations!(fap.axis)
save(joinpath(rfdir,"csu_edog_umap_img.png"),fap)

fap = GLMakie.scatter(Fz2[2,:], Fz2[1,:], marker=edogimgc, markersize=30, figure=(;resolution=(1200,1200)))
fap.axis.aspect = GLMakie.DataAspect()
GLMakie.hidespines!(fap.axis)
GLMakie.hidedecorations!(fap.axis)
save(joinpath(rfdir,"csu_edog_umap_imgc.png"),fap)

k=4
kc = kmeans(Fz',k)
clu = assignments(kc)

clu = replace(clu,1=>"eDoG-Off",3=>"eDoG-On",2=>"Off",4=>"On")
scatter(Fz2[2,:], Fz2[1,:],ms=2,msw=0,ma=0.8,group=clu,frame=:none,leg=:inline,ratio=1,size=(600,600),legendfontsize=12)
foreach(ext->savefig(joinpath(rfdir,"csu_edog_clu$ext")),figfmt)

edogcsu.spclu = clu
jldsave(joinpath(resultroot,"edogcsu.jld2");edogcsu)


# Gabor
features = [:OnOff,:aon,:von]
Fz = stack(i->zscore(gaborcsu[:,i]),features)
FzD = pairwise(Euclidean(),Fz,dims=1)

Fz2 = umap(Fz', 2, n_neighbors=45, min_dist=2, n_epochs=250,metric=Euclidean())
scatter(Fz2[2,:], Fz2[1,:],ms=2,msw=0,ma=0.8,frame=:none,leg=false,ratio=1,size=(600,600))

fap = GLMakie.scatter(Fz2[2,:], Fz2[1,:], marker=gaborimg, markersize=30, figure=(;resolution=(1200,1200)))
fap.axis.aspect = GLMakie.DataAspect()
GLMakie.hidespines!(fap.axis)
GLMakie.hidedecorations!(fap.axis)
save(joinpath(rfdir,"csu_gabor_umap_img.png"),fap)

fap = GLMakie.scatter(Fz2[2,:], Fz2[1,:], marker=gaborimgc, markersize=30, figure=(;resolution=(1200,1200)))
fap.axis.aspect = GLMakie.DataAspect()
GLMakie.hidespines!(fap.axis)
GLMakie.hidedecorations!(fap.axis)
save(joinpath(rfdir,"csu_gabor_umap_imgc.png"),fap)

k=5
kc = kmeans(Fz[:,1:2]',k)
clu = assignments(kc)

clu = replace(clu,3=>"Gabor-Even-Off",4=>"Gabor-Even-On",1=>"Gabor-Odd",2=>"On",5=>"Off")
scatter(Fz2[2,:], Fz2[1,:],ms=2,msw=0,ma=0.8,group=clu,frame=:none,leg=:inline,ratio=1,size=(600,600),legendfontsize=12)
foreach(ext->savefig(joinpath(rfdir,"csu_gabor_clu$ext")),figfmt)

gaborcsu.spclu = clu
jldsave(joinpath(resultroot,"gaborcsu.jld2");gaborcsu)



## Model Goodness and Selection
edogcsu = load(joinpath(resultroot,"edogcsu.jld2"),"edogcsu")
gaborcsu = load(joinpath(resultroot,"gaborcsu.jld2"),"gaborcsu")

gn = [:id,:rc,:PCC,:BIC,:FVU]
csumg = innerjoin(edogcsu[:,gn],gaborcsu[:,gn],on=[:id,:rc],renamecols=:edog=>:gabor)
leftjoin!(csumg,edogcsu[:,[:id,:rc,:layer,:aligndepth]],on=[:id,:rc])
dropmissing!(csumg)

@df csumg scatter(:BICgabor,:BICedog,ratio=1,leg=false,ms=2,msw=0,ma=0.7,grid=false,title="BIC")
lim=1.6e5
plot!([0,lim],[0,lim],color=:black,lw=0.5,xlim=(0,lim),ylim=(0,lim),xlabel="Gabor",ylabel="eDoG",xticks=false,yticks=false)
foreach(ext->savefig(joinpath(rfdir,"csu_BIC$ext")),figfmt)

@df csumg scatter(:FVUgabor,:FVUedog,ratio=1,leg=false,ms=2,msw=0,ma=0.7,grid=false,title="FVU")
plot!([0,1],[0,1],color=:black,lw=0.5,xlim=(0,1),ylim=(0,1),xlabel="Gabor",ylabel="eDoG")
foreach(ext->savefig(joinpath(rfdir,"csu_FVU$ext")),figfmt)

@df csumg scatter(:PCCgabor,:PCCedog,ratio=1,leg=false,ms=2,msw=0,ma=0.7,grid=false,title="PCC")
plot!([0,1],[0,1],color=:black,lw=0.5,xlim=(0,1),ylim=(0,1),xlabel="Gabor",ylabel="eDoG")
foreach(ext->savefig(joinpath(rfdir,"csu_PCC$ext")),figfmt)


transform!(csumg,[:BICgabor,:BICedog]=>((i,j)->i.<j)=>:gaborbetter)
plotdepthhist(csumg;g=:gaborbetter,layer,ylabel="Number of Spatial Model")

smcsu = sort!(append!(gaborcsu[csumg.gaborbetter,:],edogcsu[.!csumg.gaborbetter,:],cols=:union),:rc)
# density(smcsu.PCC,leg=false,xticks=0:0.1:1)
transform!(smcsu,:spclu=>ByRow(i->!startswith(i,'O'))=>:spoppo,:PCC=>(i->i.>0.6)=>:goodfit)
sgmcsu = subset(smcsu,:goodfit)

plotdepthhist(smcsu;g=:goodfit,layer,ylabel="Number of Spatial Model")
plotdepthhist(sgmcsu;g=:rc,layer,ylabel="Number of Spatial Model",palette=ALMScp)
plotdepthhist(sgmcsu;g=:spclu,layer,ylabel="Number of Spatial Model")
plotdepthhist(sgmcsu;g=:spoppo,layer,ylabel="Number of Spatial Model")

# plotdepthhist(sgmcsu;g=:rc,layer,ylabel="Number of Spatial Model",palette=ALMScp,dir=rfdir,figfmt)
# plotdepthhist(sgmcsu;g=:spclu,layer,ylabel="Number of Spatial Model",dir=rfdir,figfmt)
# plotdepthhist(sgmcsu;g=:spoppo,layer,ylabel="Number of Spatial Model",dir=rfdir,figfmt)

# foreach(siteid->plotdepthhist(sgmcsu;dir=rfdir,siteid,g=:rc,layer,figfmt,ylabel="Number of Spatial Model"),levels(sgmcsu.siteid))
# foreach(siteid->plotdepthhist(sgmcsu;dir=rfdir,siteid,g=:spclu,layer,figfmt,ylabel="Number of Spatial Model"),levels(sgmcsu.siteid))
# foreach(siteid->plotdepthhist(sgmcsu;dir=rfdir,siteid,g=:spoppo,layer,figfmt,ylabel="Number of Spatial Model"),levels(sgmcsu.siteid))

# pl = select(sgmcsu,[:spclu,:siteid,:layer]) |>
#     [@vlplot(:bar,y={"layer"},x={"count()",axis={title="Number of Spatial Model",grid=false}},color={"spclu:n",scale={scheme=:category20}});
#     @vlplot(:bar,y={"siteid"},x={"count()",axis={title="Number of Spatial Model",grid=false}},color={"spclu:n",scale={scheme=:category20}})]
# foreach(ext->save(joinpath(rfdir,"Layer_SiteID_Hist_spclu$ext"),pl),figfmt)


plotmodelcontour = (unit;siteid="AG1_V1_ODL1",cg=layercg,dir=nothing,figfmt=[".png"],title="$(siteid)_model_contour")->begin
    sunit = sort!(filter(r->r.siteid==siteid,unit),:aligndepth,rev=true) # from bottom to top of cortex
    sizedeg = sunit.sizedeg[1] # stimulus diameters
    # standard 2D coordinates(rightward:x, upward:y)
    stix = [0,1,1,0,0] * sizedeg[2]
    stiy = [1,1,0,0,1] * sizedeg[1]
    α = 0:0.02:2π
    p = plot(stix,stiy;lw=2,color=:black,frame=:grid,ratio=1,tickdir=:out,xlabel="X (deg)",ylabel="Y (deg)",leg=false)
    for r in eachrow(sunit)
        # roi is defined on pixel indices of stimulus image
        c = (r.roi.centerdeg[2], sizedeg[1] - r.roi.centerdeg[1])
        points = map(i->r.fit.cfun(i,r.fit.param) .+ c,α)
        plot!(p,points;la=0.7,lw=0.8,color=cg[r.aligndepth])
    end
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotmodelcontour(sgmcsu,siteid="AG1_V1_ODL1")
foreach(siteid->plotmodelcontour(sgmcsu;siteid,dir=rfdir,figfmt),levels(sgmcsu.siteid))


plotdepthscatter = (unit;x=:Diameter,xfun=first,g=:rc,dir=nothing,siteid="All",figfmt=[".png"],layer=nothing,xlabel="$x",ts="",
            ylabel="Normalized Cortical Depth",wfun=x->isempty(x) ? 0 : median(x),size=(350,550),palette=:tab20,leg=:best) -> begin
    sunit = siteid=="All" ? unit : filter(r->r.siteid==siteid,unit)
    title="$(siteid)_Depth-$(x)_$g$ts"
    xx = xfun.(sunit[!,x])
    xmin,xmax = extrema(xx)
    xr = xmax-xmin
    xlims = (xmin-0.08xr,xmax+0.01xr)
    annx = xlims[1]+0.01(xlims[2]-xlims[1])
    p = scatter(xx,sunit.aligndepth;group=sunit[!,g],yflip=true,grid=false,legendfontsize=6,xlims,
        ylims=(-0.02,1.02),yticks=0:0.2:1,ms=2,msw=0,ma=0.7,size,tickor=:out,ylabel,xlabel,palette,leg,title,titlefontsize=10)

    n,y = unitdensity(sunit.aligndepth;w=xx,wfun,spacerange=(0,1),bw=0.02,step=0.01)
    plot!(p,n,y,color=:gray40,la=0.8,lw=1,label="Average")
    if !isnothing(layer)
        ann = [(annx,mean(layer[k]),Plots.text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
    end
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotdepthscatter(sgmcsu;x=:Diameter,g=:spclu,layer,xlabel="Diameter (deg)")
plotdepthscatter(sgmcsu;x=:Diameter,g=:rc,layer,xlabel="Diameter (deg)",palette=ALMScp)
plotdepthscatter(sgmcsu;x=:Round,g=:spclu,layer)
plotdepthscatter(sgmcsu;x=:Round,g=:rc,layer,palette=ALMScp)
plotdepthscatter(sgmcsu;x=:dc,g=:spclu,layer,xlabel="Delay (ms)")
plotdepthscatter(sgmcsu;x=:dc,g=:rc,layer,xlabel="Delay (ms)",palette=ALMScp)

# plotdepthscatter(sgmcsu;x=:Diameter,g=:spclu,layer,xlabel="Diameter (deg)",dir=rfdir,figfmt)
# plotdepthscatter(sgmcsu;x=:Diameter,g=:rc,layer,xlabel="Diameter (deg)",dir=rfdir,figfmt,palette=ALMScp)
# plotdepthscatter(sgmcsu;x=:Round,g=:spclu,layer,dir=rfdir,figfmt)
# plotdepthscatter(sgmcsu;x=:Round,g=:rc,layer,dir=rfdir,figfmt,palette=ALMScp)
# plotdepthscatter(sgmcsu;x=:dc,g=:spclu,layer,xlabel="Delay (ms)",dir=rfdir,figfmt)
# plotdepthscatter(sgmcsu;x=:dc,g=:rc,layer,xlabel="Delay (ms)",dir=rfdir,figfmt,palette=ALMScp)



## Spectral-Spatial RF Model
ssmcsu = rightjoin(combine(groupby(sgmcsu,:id),:rc=>join=>:MRTYPE,
        [:rc,:Diameter]=>((c,v)->OrderedDict((c.=>v)...))=>:Diameter,
        [:rc,:θ]=>((c,v)->OrderedDict((c.=>v)...))=>:θ,
        [:rc,:μx]=>((c,v)->OrderedDict((c.=>v)...))=>:μx,
        [:rc,:μy]=>((c,v)->OrderedDict((c.=>v)...))=>:μy,
        [:rc,:V]=>((c,v)->OrderedDict((c.=>v)...))=>:V,
        :μx=>mean=>:ssmμx,
        :μy=>mean=>:ssmμy,
        [:μx,:μy,:Diameter]=>((x,y,d)->maximum([abs.(x.-mean(x)).+0.6d;abs.(y.-mean(y)).+0.6d]))=>:ssmradius,
        [:rc,:fit]=>((c,v)->OrderedDict((c.=>v)...))=>:ssm,
        [:rc,:spclu]=>((c,p)->join(c.*'_'.*p,", "))=>:ssp,
        [:rc,:spoppo]=>((c,o)->join(c.*'_'.*replace(o,true=>"O",false=>"N"),", "))=>:sspo),
        stacsu,on=:id)
transform!(ssmcsu,:MRTYPE=>ByRow(i->!ismissing(i))=>:mresponsive)
ssgmcsu = dropmissing(ssmcsu,:MRTYPE)
transform!(ssgmcsu,:MRTYPE=>ByRow(i->replace(i,r"A([L,M,S]+)"=>s"\1"))=>:MCTYPE,
        :ssp=>ByRow(i->replace(i,r"A_\S+, ([L,M,S]_.+)"=>s"\1"))=>:cssp,
        :sspo=>ByRow(i->replace(i,r"A_\S+, ([L,M,S]_.+)"=>s"\1"))=>:csspo,
        :V=>ByRow(d->begin
            ub = reduce(max,values(d))
            OrderedDict(k=>v/ub for (k,v) in d)
        end)=>:nV,
        :V=>ByRow(d->begin
            L = haskey(d,'L') ? d['L'] : 0
            M = haskey(d,'M') ? d['M'] : 0
            S = haskey(d,'S') ? d['S'] : 0
            vl = missing; vm = missing
            if !(L==M==S==0)
                Σ = L+M+S
                vl = L/Σ; vm = M/Σ
            end
            (;vl, vm)
        end)=>AsTable)

plotdepthhist(ssmcsu;g=:mresponsive,layer)
plotdepthhist(ssgmcsu;g=:MRTYPE,layer)
plotdepthhist(ssgmcsu;g=:MCTYPE,layer)
jldsave(joinpath(resultroot,"ssmcsu.jld2");ssmcsu)

# plotdepthhist(ssmcsu;g=:mresponsive,layer,dir=rfdir,figfmt)
# plotdepthhist(ssgmcsu;g=:MRTYPE,layer,dir=rfdir,figfmt)
# plotdepthhist(ssgmcsu;g=:MCTYPE,layer,dir=rfdir,figfmt)

# pl = select(ssgmcsu,[:MRTYPE,:siteid,:layer]) |>
#     [@vlplot(:bar,y={"layer"},x={"count()",axis={title="Number of Units",grid=false}},color={"MRTYPE:n",scale={scheme=:category20}});
#     @vlplot(:bar,y={"siteid"},x={"count()",axis={title="Number of Units",grid=false}},color={"MRTYPE:n",scale={scheme=:category20}})]
# foreach(ext->save(joinpath(rfdir,"Layer_SiteID_Hist_MRTYPE$ext"),pl),figfmt)

# pl = select(ssgmcsu,[:ssp,:sspo,:cssp,:csspo]) |> 
#     [@vlplot(:bar,y={"count()",axis={title="Number of Units"}},x={"ssp"});
#     @vlplot(:bar,y={"count()",axis={title="Number of Units"}},x={"cssp"});
#     @vlplot(:bar,y={"count()",axis={title="Number of Units"}},x={"sspo"});
#     @vlplot(:bar,y={"count()",axis={title="Number of Units"}},x={"csspo"})]
# foreach(ext->save(joinpath(rfdir,"ssgmcsu_ssp$ext"),pl),figfmt)


function plotssmpair(ssm,mp;w=350,h=350,dlim=nothing,clim=nothing,ms=2.5,ma=0.8)
    pssm = filter(r->contains(r.MRTYPE,mp.first) && contains(r.MRTYPE,mp.second),ssm)

    p = plot(layout=(3,1),leg=false,size=(w,3h))
    x = getindex.(pssm[!,:Diameter],mp.first)
    y = getindex.(pssm[!,:Diameter],mp.second)
    lim = isnothing(dlim) ? max(maximum(x),maximum(y)) : dlim
    scatter!(p[1],x,y;xlabel="$(mp.first) (deg)",ylabel="$(mp.second) (deg)",title="Diameter",ms,msw=0,ma,ratio=1)
    plot!(p[1],[0,lim],[0,lim],lw=0.5,color=:black)

    x = rad2deg.(getindex.(pssm[!,:θ],mp.first))
    y = rad2deg.(getindex.(pssm[!,:θ],mp.second))
    xticks=yticks=0:30:180
    scatter!(p[2],x,y;xlabel="$(mp.first) (deg)",ylabel="$(mp.second) (deg)",title="Orientation",xticks,yticks,ms,msw=0,ma,ratio=1)
    plot!(p[2],[0,180],[0,180],lw=0.5,color=:black)

    x = getindex.(pssm[!,:μx],mp.first) .- getindex.(pssm[!,:μx],mp.second)
    y = getindex.(pssm[!,:μy],mp.first) .- getindex.(pssm[!,:μy],mp.second)
    lim = isnothing(clim) ? maximum(abs.([x;y])) : clim
    xlims=ylims=[-lim,lim]
    scatter!(p[3],x,y;xlabel="$(mp.first)x - $(mp.second)x (deg)",ylabel="$(mp.first)y - $(mp.second)y (deg)",title="Center Displacement",
    ms,msw=0,ma,ratio=1,xlims,ylims)
    vline!(p[3],[0],lw=0.5,color=:black)
    hline!(p[3],[0],lw=0.5,color=:black)
end

plotssmpair(ssgmcsu,'A'=>'L')
plot((plotssmpair(ssgmcsu,i=>j,dlim=6,clim=0.5) for (i,j) in combinations(['A','L','M','S'],2))...,layout=(1,6),size=(6*350,3*350))
foreach(ext->savefig(joinpath(rfdir,"ssgmcsu_affine$ext")),figfmt)


# Model Volumn as Cone Weight
csgmcsu = dropmissing(ssgmcsu,[:vl,:vm]) # exclude 'A' only cell
plot([0,1],[1,0],ratio=1,leg=false,xlims=(-0.01,1.01),ylims=(-0.01,1.01),xticks=0:0.2:1,yticks=0:0.2:1,lw=1,color=:black,tickor=:out)
scatter!(csgmcsu.vl,csgmcsu.vm,ms=3,ma=0,msw=1,msa=1,msc=:dodgerblue,xlabel="L",ylabel="M",frame=:origin,grid=false)
foreach(ext->savefig(joinpath(rfdir,"csgmcsu_rcw$ext")),figfmt)


function plotssm(ssm;n=5,rand=false,ni=nothing,type=:image,ncol=5,h=130,w=type==:image ? 450 : 130,path=nothing,bg=:gray,figfmt=[".png"],
                size=(61,61),titlefontsize=type==:image ? 14 : 6,hidetitle=false,cch=['A','L','M','S'],clf=0.1)
    nu = nrow(ssm); n = min(n,nu); n==0 && return
    nc = n < ncol ? n : ncol
    nr = ceil(Int,n/ncol)
    ni = rand ? sample(1:nu,n,replace=false) : isnothing(ni) ? (1:n) : ni

    p=plot(;layout=(nr,nc),size=(nc*w,nr*h),frame=:none,leg=false,titlefontsize)
    for i in eachindex(ni)
        title = hidetitle ? "" : ssm.id[ni[i]]
        ssmim = ssm.ssmim[ni[i]]
        cs = intersect(cch,keys(ssmim))
        if type == :image
            img = hcat(map(c->maskimg(ssmim[c].img,c;size,masktype="None"),cs)...)
            plot!(p[i],img;title,bginside=bg)
        elseif type == :contour
            x = 1:size[2]; y = 1:size[1]; xlims=(x[begin],x[end]);ylims=(y[begin],y[end])
            for c in cs
                img = imresize(ssmim[c].img,size)
                lb,ub = extrema(img)
                lim = max(abs(lb),abs(ub))
                levels=clf*[lb,ub]
                colors = get(cgrad(cmcode[c].colors),[lb,ub],(-lim,lim))
                for j in 1:2
                    for l in Contour.lines(Contour.contour(y,x,img,levels[j]))
                        ys,xs = Contour.coordinates(l)
                        plot!(p[i],xs,ys;color=colors[j],yflip=true,ratio=1,xlims,ylims,lw=1.5,title,bginside=bg)
                    end
                end
            end
        end
    end
    isnothing(path) && return p
    mkpath(dirname(path));foreach(ext->savefig("$path.$type$ext"),figfmt)
end

transform!(ssgmcsu,[:ssm,:ssmμx,:ssmμy,:ssmradius]=>ByRow((m,μx,μy,r)->OrderedDict(k=>fitimage(v,μx,μy,r) for (k,v) in m))=>:ssmim)

plotssm(ssgmcsu;n=13,rand=true)
plotssm(ssgmcsu;n=23,rand=true,type=:contour)

foreach(t->plotssm(filter(r->r.sspo == t,ssgmcsu);n=nrow(ssgmcsu),hidetitle=true,type=:contour,path=joinpath(rfdir,t),figfmt),
        levels(ssgmcsu.sspo))



## Spectral-Spatial RF Model Mixing
ssmmix = (ssmim) -> begin
    cs = collect(keys(ssmim))
    n = length(cs)
    n == 1 && return missing
    mp=OrderedDict()
    for (i,j) in combinations(cs,2)
        c = i*j
        imgi,maski=ssmim[i]
        imgj,maskj=ssmim[j]
        ar = log2(count(maski)/count(maskj)) # ratio of RF area
        o = maski.&maskj
        Ω = maski.|maskj
        oi = imgi[o]; oj = imgj[o]
        ao = count(o)/count(Ω) # ratio of overlapping area to combined RF area
        or = cor(oi,oj) # correlation in overlapping area
        oass = mmmetric(oi.*oj).aon # ratio of same sign area in overlapping area
        ovr = log2(sum(abs.(oi))/sum(abs.(oj))) # ratio of volumn on overlapping area
        mp[c] = (;ar,ao,or,oass,ovr,oi,oj)
    end
    n == 2 && return mp
    for (i,j,k) in combinations(cs,3)
        c = i*j*k
        imgi,maski=ssmim[i]
        imgj,maskj=ssmim[j]
        imgk,maskk=ssmim[k]
        ai = count(maski);aj = count(maskj);ak = count(maskk)
        ra = (;Symbol(i)=>ai/(ai+aj+ak),Symbol(j)=>aj/(ai+aj+ak)) # relative RF areas
        o = maski.&maskj.&maskk
        Ω = maski.|maskj.|maskk
        oi = imgi[o]; oj = imgj[o]; ok = imgk[o]
        ao = count(o)/count(Ω)
        opr = (;Symbol(i*j)=>partialcor(oi,oj,ok),Symbol(i*k)=>partialcor(oi,ok,oj)) # partial correlations in overlapping area
        vi=sum(abs.(oi));vj=sum(abs.(oj));vk=sum(abs.(ok))
        orv = (;Symbol(i)=>vi/(vi+vj+vk),Symbol(j)=>vj/(vi+vj+vk)) # relative volumns on overlapping area
        mp[c] = (;ra,ao,opr,orv,oi,oj,ok)
    end
    mp
end
mixparam(unit,key) = transform!(subset!(dropmissing(unit,:ssmmix),:ssmmix=>ByRow(d->haskey(d,key))),:ssmmix=>ByRow(d->d[key])=>AsTable)
transform!(ssgmcsu,:ssmim=>ByRow(ssmmix)=>:ssmmix)

lm = mixparam(ssgmcsu,"LM")
@df lm corrplot(cols([:ar,:ao,:or,:oass,:ovr]),bin=50,leg=false,size=(850,650),grid=false)

# foreach(k->begin
#     @df mixparam(ssgmcsu,k) corrplot(cols([:ar,:ao,:or,:oass,:ovr]),bin=50,leg=false,size=(850,650),grid=false)
#     foreach(ext->savefig(joinpath(rfdir,"ssgmcsu_$(k)_mix$ext")),figfmt)
# end,["LM","LS","MS","AL","AM","AS"])

plotdepthscatter(lm;x=:ar,g=:MRTYPE,layer)
plotdepthscatter(lm;x=:ao,g=:MRTYPE,layer)
plotdepthscatter(lm;x=:or,g=:MRTYPE,layer)
plotdepthscatter(lm;x=:oass,g=:MRTYPE,layer)
plotdepthscatter(lm;x=:ovr,g=:MRTYPE,layer)

# foreach(k->begin
#     plot((plotdepthscatter(mixparam(ssgmcsu,k);x,g=:MRTYPE,layer) for x in [:ar,:ao,:or,:oass,:ovr])...,layout=(1,5),size=(5*350,550))
#     foreach(ext->savefig(joinpath(rfdir,"ssgmcsu_$(k)_DepthMix$ext")),figfmt)
# end,["LM","LS","MS","AL","AM","AS"])

lms = mixparam(ssgmcsu,"LMS")
histogram(lms.ao,bins=50,xlims=(0,1),lc=:match,leg=false)
plotdepthscatter(lms;x=:ao,g=:MRTYPE,layer)
plotdepthscatter(lms;x=:ra,g=:MRTYPE,layer,xlabel="L Relative Area")
plotdepthscatter(lms;x=:ra,xfun=last,g=:MRTYPE,layer,xlabel="M Relative Area")
plotdepthscatter(lms;x=:orv,g=:MRTYPE,layer,xlabel="L Relative Volumn")
plotdepthscatter(lms;x=:orv,xfun=last,g=:MRTYPE,layer,xlabel="M Relative Volumn")
plotdepthscatter(lms;x=:opr,g=:MRTYPE,layer,xlabel="LM Partial Correlation")
plotdepthscatter(lms;x=:opr,xfun=last,g=:MRTYPE,layer,xlabel="LS Partial Correlation")

# foreach(k->plotdepthscatter(mixparam(ssgmcsu,k);x=:ao,g=:MRTYPE,layer,dir=rfdir,figfmt,ts="_$k"),["ALM","ALS","AMS","LMS"])

function scattermarginalhist(points;bins=50,size=(650,650),lims=(0,1),rline=true,frame=:origin,ms=4,msw=0,ma=1,msa=1,mc=:dodgerblue,msc=:dodgerblue)
    layout = @layout [a       _ 
                      b{0.9w,0.9h} c]
    xlims=ylims=lims
    p=plot(;layout,leg=false,size,tickor=:out)
    rline && plot!(p[2,1],[0,1],[1,0];lw=1,color=:black)
    scatter!(p[2,1],points;ms,msw,ma,msa,mc,msc,frame,ratio=1,xlims,ylims)
    histogram!(p[1,1],first.(points);bins,frame=:none,xlims,lc=:match,color=:gray20)
    histogram!(p[2,2],last.(points);bins,dir=:h,frame=:none,ylims,lc=:match,color=:gray20)
    p
end

scattermarginalhist(lms.ra,mc=layercg[lms.aligndepth])
scattermarginalhist(lms.orv,mc=layercg[lms.aligndepth])
scattermarginalhist(lms.opr,lims=(-1,1),rline=false,frame=:zerolines,mc=layercg[lms.aligndepth])

# foreach(k->begin
#     t = mixparam(ssgmcsu,k)
#     plot(scattermarginalhist(t.ra,mc=layercg[t.aligndepth]),
#     scattermarginalhist(t.orv,mc=layercg[t.aligndepth]),
#     scattermarginalhist(t.opr,lims=(-1,1),rline=false,frame=:zerolines,mc=layercg[t.aligndepth]),
#     layout=(1,3),size=(3*650,650))
#     foreach(ext->savefig(joinpath(rfdir,"ssgmcsu_$(k)_Mix$ext")),figfmt)
# end,["ALM","ALS","AMS","LMS"])

# i=25
# plotssm(lm[i:i,:],cch=['L','M'],type=:contour,hidetitle=true)
# scatter(lm.oi[i],lm.oj[i],ms=2,msw=0,leg=false,ratio=1,frame=:origin,xticks=[0],yticks=[0])



## Spectral-Spatial Model Mixing Clustering
function scatterssm(ssm;posx,posy,resolution=(1200,1200),type=:contour,path=nothing,figfmt=[".png"],pr = [-0.5,0,0.5],aas=0.5,xlabel="",ylabel="",
                    size=(51,51),cch=['A','L','M','S'],bg=:white,clf=0.1,yflip=false,ms=30,frame=:none,masktype="Gaussian")
    f = GLMakie.Figure(;resolution,backgroundcolor=bg)
    ax = GLMakie.Axis(f[1,1];yreversed=yflip,backgroundcolor=bg,xlabel,ylabel)
    ssmimg=[]
    x = 1:size[2]; y = 1:size[1]
    for i in 1:nrow(ssm)
        ssmim = ssm.ssmim[i]
        cs = intersect(cch,keys(ssmim))
        if type == :image
            push!(ssmimg, mapreduce(c->maskimg(ssmim[c].img,c;size,masktype),(i,j)->weighted_color_mean.(0.5,i,j),cs))
        elseif type == :contour
            cimg = fill(RGBA{Float32}(0,0,0,0),size)
            for c in cs
                img = imresize(ssmim[c].img,size)
                lb,ub = extrema(img)
                lim = max(abs(lb),abs(ub))
                levels=clf*[lb,ub]
                colors = float32.(get(cgrad(cmcode[c].colors),[lb,ub],(-lim,lim)))
                for j in eachindex(levels)
                    for l in Contour.lines(Contour.contour(y,x,img,levels[j]))
                        ys,xs = Contour.coordinates(l)
                        foreach((y,x)->cimg[clamp.(round.(Int,y.+pr),1,size[1]),clamp.(round.(Int,x.+pr),1,size[2])].= colors[j],ys,xs)
                    end
                end
            end
            push!(ssmimg, imfilter(cimg,Kernel.gaussian(aas)))
        end
    end
    GLMakie.scatter!(ax,posx,posy,marker=ssmimg,markersize=ms)
    if frame == :none
        GLMakie.hidedecorations!(ax)
        GLMakie.hidespines!(ax)
    end
    isnothing(path) && return f
    mkpath(dirname(path));foreach(ext->save("$path.$type$ext",f),figfmt)
end

t=lm
features = [:ar,:ao,:or,:oass,:ovr]
Fz = stack(f->zscore(t[:,f]),features)

t=lms
features = [:ao,:ra,:opr,:orv]
F = hcat(t[!,features[1]],map(i->stack(t[:,i])',features[2:end])...)
Fz = stack(zscore,eachslice(F,dims=2))

Fz2 = umap(Fz', 2, n_neighbors=25, min_dist=0.5, n_epochs=250,metric=Euclidean())
scatter(Fz2[2,:], Fz2[1,:],ms=2,msw=0,ma=0.8,frame=:none,leg=false,ratio=1,size=(600,600))

k=2
kc = kmeans(Fz',k)
clu = assignments(kc)

scatter(Fz2[2,:], Fz2[1,:],ms=2,msw=0,ma=0.8,group=clu,frame=:none,leg=:inline,ratio=1,size=(600,600),legendfontsize=12)
scatter(Fz2[2,:], Fz2[1,:],ms=2,msw=0,ma=0.8,color=get(cgrad(:coolwarm),t.or,(-1,1)),frame=:none,leg=true,ratio=1,size=(600,600),legendfontsize=12)



scatterssm(lm,posx=Fz2[2,:],posy=Fz2[1,:],cch=['L','M'],ms=40,size=(100,100))

scatterssm(lm,posy=lm.aligndepth,posx=lm.or,cch=['L','M'],resolution=(500,1000),frame=:grid,ylabel="depth")

scatterssm(lm,posx=Fz2[2,:],posy=Fz2[1,:],cch=['L','M'],type=:image,bg=:gray)


function mixclass(csspo, ssmmix)
    c1 = c2 = "N"
    if contains(csspo, ", ")
        c = ['L', 'M', 'S']
        ci = contains.(csspo, c)
        k = join(c[ci])
        mp = ssmmix[k]
        if length(k) == 3
            lm = mp.opr.LM >= 0 ? "/" : "\\"
            ls = mp.opr.LS >= 0 ? "/" : "\\"
            c1 = "3$lm$ls"
            c2 = join(c[ci], lm, ls)
        else
            r = mp.or >= 0 ? "/" : "\\"
            c1 = "2$r"
            c2 = join(c[ci], r)
        end
    else
        if contains(csspo, "A")
            c1 = c2 = "A"
        else
            if contains(csspo, "N")
                c1 = "0"
            else
                c1 = "1"
            end
            c2 = "$(csspo[1])$c1"
        end
    end
    (mixclass1=c1, mixclass2=c2)
end

transform!(ssgmcsu,[:csspo,:ssmmix]=>ByRow(mixclass)=>AsTable)
jldsave(joinpath(resultroot,"ssgmcsu.jld2");ssgmcsu)

plotdepthhist(ssgmcsu;g=:mixclass2,layer)

# plotdepthhist(ssgmcsu;g=:mixclass1,layer,dir=rfdir,figfmt)
# plotdepthhist(ssgmcsu;g=:mixclass2,layer,dir=rfdir,figfmt)

# pl = select(ssgmcsu,[:mixclass1,:siteid,:layer]) |>
#     [@vlplot(:bar,y={"layer"},x={"count()",axis={title="Number of Units",grid=false}},color={"mixclass1:n",scale={scheme=:category20}});
#     @vlplot(:bar,y={"siteid"},x={"count()",axis={title="Number of Units",grid=false}},color={"mixclass1:n",scale={scheme=:category20}})]
# foreach(ext->save(joinpath(rfdir,"Layer_SiteID_Hist_mixclass1$ext"),pl),figfmt)

# pl = select(ssgmcsu,[:mixclass2,:siteid,:layer]) |>
#     [@vlplot(:bar,y={"layer"},x={"count()",axis={title="Number of Units",grid=false}},color={"mixclass2:n",scale={scheme=:category20}});
#     @vlplot(:bar,y={"siteid"},x={"count()",axis={title="Number of Units",grid=false}},color={"mixclass2:n",scale={scheme=:category20}})]
# foreach(ext->save(joinpath(rfdir,"Layer_SiteID_Hist_mixclass2$ext"),pl),figfmt)


