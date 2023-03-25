using NeuroAnalysis,DataFrames,JLD2,StatsBase,StatsPlots,LinearAlgebra,Images,UMAP,
    Clustering,Distances,VegaLite,Combinatorics,XLSX,DataStructures
import GLMakie,Contour

## Collect All STAs
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
            dc = map((u,r)->r ? map((d,cr)->cr ? delays[d.d] : missing,dataset["ulcsdd"][u],dataset["ucresponsive"][u]) : missing,unitid,unitresponsive)
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
figfmt = [".png"]

# staunit = collectsta(resultroot)
# jldsave(joinpath(resultroot,"staunit.jld2");staunit)

staunit = load(joinpath(resultroot,"staunit.jld2"),"staunit")
allcu = load(joinpath(resultroot,"allcu.jld2"),"allcu")
layer = load(joinpath(resultroot,"layertemplate.jld2"),"layertemplate")
penetration = DataFrame(XLSX.readtable(joinpath(resultroot,"penetration.xlsx"),"Sheet1"))
penetration = select(penetration,[:siteid,:od,:cofd],:cofd=>ByRow(i->i∈["L/M","S/LM","L/M, S/LM"])=>:incofd,
                :cofd=>ByRow(i->i∈["B","W","B, W"])=>:inbw,:cofd=>ByRow(i->i=="None")=>:none)

## STA responsive
stacu = innerjoin(allcu,staunit,on=[:siteid,:id])
leftjoin!(stacu,penetration,on=:siteid)
transform!(stacu,:rc=>ByRow(i->ismissing(i) ? missing : join(skipmissing(i)))=>:RTYPE)
stacu.CTYPE = map(i->ismissing(i) ? missing : replace(i,r"A([L,M,S]+)"=>s"\1"),stacu.RTYPE)
stacsu = subset(stacu,:good)
starcsu = dropmissing(stacsu,:rc)
stanrcsu = subset(stacsu,:responsive=>.!)
jldsave(joinpath(resultroot,"stacsu.jld2");stacsu)

plotlstaroi = (unit;siteid="AG1_V1_ODL1",cg=cgrad(:matter),dir=nothing,figfmt=[".png"],title="$(siteid)_lsta_roi")->begin
    cunit = sort!(filter(r->r.siteid==siteid,unit),:aligndepth,rev=true)
    sizedeg = cunit.sizedeg[1]
    stimx = [0,1,1,0,0] * sizedeg[1]
    stimy = [1,1,0,0,1] * sizedeg[2]
    # roi is derived from sta image with pixel origin at topleft corner
    x = map(i->i.centerdeg[2],cunit.roi)
    y = map(i->sizedeg[2]-i.centerdeg[1],cunit.roi)
    rx = map(i->i.radiideg[2],cunit.roi)
    ry = map(i->i.radiideg[1],cunit.roi)
    roix = [-1,1,1,-1,-1] .* rx' .+ x'
    roiy = [1,1,-1,-1,1] .* ry' .+ y'
    p = plot(stimx,stimy,lw=1,color=:black,frame=:grid,tickdir=:out,ratio=1,xlabel="X (deg)",ylabel="Y (deg)",leg=false)
    plot!(p,roix,roiy;la=0.8,color=cg[cunit.aligndepth]')
    isnothing(dir) && return p
    mkpath(dir);foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotlstaroi(starcsu;siteid="AG2_V1_ODR13")
foreach(siteid->plotlstaroi(starcsu;siteid,dir=joinpath(resultroot,"RF")),levels(starcsu.siteid))


plotdepthhist = (unit;g=:responsive,dir=nothing,figfmt=[".png"],layer=nothing,title="CorticalDepth_$g",ylabel="Number of Units",palette=:tab20,leg=:best) -> begin
    p = groupedhist(unit.aligndepth;group=unit[!,g],barpositions=:stack,bin=0:0.02:1,permute=(:y,:x),xflip=true,grid=false,legendfontsize=6,
        xlims=(0,1),lw=0,size=(350,500),tickor=:out,xlabel="Cortical Depth",ylabel,palette,leg)
    if !isnothing(layer)
        ann = [(1,mean(layer[k]),Plots.text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
    end
    isnothing(dir) && return p
    mkpath(dir);foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

# foreach(siteid->plotdepthhist(filter(r->r.siteid==siteid,stacsu);dir=resultroot,g=:responsive,title="$(siteid)_CorticalDepth_$g",layer),levels(stacsu.siteid))
# foreach(siteid->plotdepthhist(filter(r->r.siteid==siteid,stacsu);dir=resultroot,g=:RTYPE,title="$(siteid)_CorticalDepth_$g",layer),levels(stacsu.siteid))

plotdepthhist(stacsu;g=:responsive,layer)
plotdepthhist(starcsu;g=:RTYPE,layer)
plotdepthhist(starcsu;g=:CTYPE,layer)

p=select(starcsu,[:RTYPE,:siteid,:layer]) |>
[@vlplot(:bar,y={"RTYPE"},x={"count()"});
@vlplot(:bar,y={"siteid"},x={"count()"},color={"RTYPE:n",scale={scheme=:category20}});
@vlplot(:bar,y={"layer"},x={"count()"},color={"RTYPE:n",scale={scheme=:category20}})]

# COFD, Achromatic and None
dn = :none
plotdepthhist(subset(stacsu,dn);g=:responsive,layer)
plotdepthhist(subset(starcsu,dn);g=:RTYPE,layer)

p=select(subset(starcsu,dn),[:RTYPE,:layer]) |>
[@vlplot(:bar,x={"RTYPE"},y={"count()"});
@vlplot(:bar,y={"layer"},x={"count()"},color={"RTYPE:n",scale={scheme=:category20}})]


## Spatial Pattern
staedogunit = (unit;model=:edog) -> begin
    select!(dropmissing!(flatten(dropmissing(unit,:rc),["rc","dc","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->(;i.model,i.fun,i.mfun,i.cfun,i.param,i.radii))=>:fit,
        "sta!$model"=>ByRow(i->i.r)=>:cor,
        "sta!$model"=>ByRow(i->i.bic)=>:bic,
        "sta!$model"=>ByRow(i->1-i.r2)=>:fvu,
        "sta!$model"=>ByRow(i->i.param[1])=>:ae,
        "sta!$model"=>ByRow(i->i.param[7])=>:ai,
        "sta!$model"=>ByRow(i->i.param[2])=>:μx,
        "sta!$model"=>ByRow(i->i.param[4])=>:μy,
        "sta!$model"=>ByRow(i->i.param[3])=>:σxe,
        "sta!$model"=>ByRow(i->i.param[5])=>:r,
        "sta!$model"=>ByRow(i->i.param[6])=>:θ,
        "sta!$model"=>ByRow(i->i.param[8])=>:ρ)
end
stagaborunit = (unit;model=:gabor) -> begin
    select!(dropmissing!(flatten(dropmissing(unit,:rc),["rc","dc","sta!$model"]),"sta!$model"),Not(r"\w*!\w*"),
        "sta!$model"=>ByRow(i->(;i.model,i.fun,i.mfun,i.cfun,i.param,i.radii))=>:fit,
        "sta!$model"=>ByRow(i->i.r)=>:cor,
        "sta!$model"=>ByRow(i->i.bic)=>:bic,
        "sta!$model"=>ByRow(i->1-i.r2)=>:fvu,
        "sta!$model"=>ByRow(i->i.param[1])=>:a,
        "sta!$model"=>ByRow(i->i.param[2])=>:μx,
        "sta!$model"=>ByRow(i->i.param[4])=>:μy,
        "sta!$model"=>ByRow(i->i.param[3])=>:σx,
        "sta!$model"=>ByRow(i->i.param[5])=>:r,
        "sta!$model"=>ByRow(i->i.param[6])=>:θ,
        "sta!$model"=>ByRow(i->i.param[7])=>:f,
        "sta!$model"=>ByRow(i->i.param[8])=>:p)
end
function fitpredict(fit;ppd=50,rp=170,tight=false,metric=false)
   if tight
       param = deepcopy(fit.param)
       ppd = max(ppd,rp/min(fit.radii...)) # variable ppd to ensure rp pixels for local roi radius for all colors
       if fit.model == :edog
           param[[2,4]].=0
           σxe = param[3]
           r = param[5]
           ρ = param[8]
           radius = 3*max(σxe,r*σxe,ρ*σxe,ρ*r*σxe) # img radius for predict
           fit = (;fit.model,fit.fun,fit.mfun,param,radii=(radius,radius))
       elseif fit.model == :gabor
           param[[2,4]].=0
           σx = param[3]
           r = param[5]
           radius = 3*max(σx,r*σx)
           fit = (;fit.model,fit.fun,fit.mfun,param,radii=(radius,radius))
       end
   end
   y = -fit.radii[1]:1/ppd:fit.radii[1]
   x = -fit.radii[2]:1/ppd:fit.radii[2]
   img = predict(fit,x,y,xygrid=true,yflip=true)
   if metric
       mask = tight ? predict((;fun=fit.mfun,param),x,y,xygrid=true,yflip=true) : trues(size(img))
       onap,onwp,W=mimetric(img[mask])
       return (;img,onap,onwp,W)
   else
       return  img
   end
end
function mimetric(img)
    on = img.>0
    onap = count(on)/length(img)
    W = sum(abs.(img))
    onwp = sum(img[on])/W
    return (;onap,onwp,W)
end
function fitpredict(fit,μx,μy,r;rp=170,ppd=max(50,rp/r))
    param = deepcopy(fit.param)
    param[2]=param[2]-μx
    param[4]=param[4]-μy
    fit = (;fit.model,fit.fun,fit.mfun,param,radii=(r,r))

    y = -fit.radii[1]:1/ppd:fit.radii[1]
    x = -fit.radii[2]:1/ppd:fit.radii[2]
    img = predict(fit,x,y,xygrid=true,yflip=true)
    mask = predict((;fun=fit.mfun,param),x,y,xygrid=true,yflip=true)
    (;img,mask)
end


edogcsu = load(joinpath(resultroot,"edogcsu.jld2"),"edogcsu")
gaborcsu = load(joinpath(resultroot,"gaborcsu.jld2"),"gaborcsu")

edogcsu = staedogunit(stacsu)
transform!(edogcsu,
    [:σxe,:r,:ρ]=>ByRow((σxe,r,ρ)->5*max(σxe,r*σxe,ρ*σxe,ρ*r*σxe))=>:diameter,
    [:ae,:ai]=>ByRow(/)=>:ρa,
    :r=>ByRow(log2∘inv)=>:el,
    :r=>ByRow(abs∘log2∘inv)=>:rd,
    :ρ=>ByRow(log2)=>:cs,
    [:ae,:ai]=>ByRow(log2∘/)=>:onoff,
    [:ae,:ai]=>ByRow(sign∘log2∘/)=>:onoffsign,
    :fit=>ByRow(i->fitpredict(i,tight=true,metric=true))=>[:img,:onap,:onwp,:W])

gaborcsu = stagaborunit(stacsu)
transform!(gaborcsu,
    [:σx,:r]=>ByRow((σx,r)->5*max(σx,r*σx))=>:diameter,
    [:σx,:r,:f]=>ByRow((σx,r,f)->5*r*σx*f)=>:cyc,
    :r=>ByRow(log2∘inv)=>:el,
    :r=>ByRow(abs∘log2∘inv)=>:rd,
    :p=>ByRow(p->sin(2π*p))=>:onoff,
    :p=>ByRow(p->(sign∘sin)(2π*p))=>:onoffsign,
    :p=>ByRow(p->(abs∘sin)(2π*p))=>:oddeven,
    :fit=>ByRow(i->fitpredict(i,tight=true,metric=true))=>[:img,:onap,:onwp,:W])

gn = [:id,:rc,:cor,:bic,:fvu]
csumg = innerjoin(edogcsu[:,gn],gaborcsu[:,gn],on=[:id,:rc],renamecols=:edog=>:gabor)
leftjoin!(csumg,edogcsu[:,[:id,:rc,:layer,:aligndepth,:incofd,:inbw,:none]],on=[:id,:rc])
dropmissing!(csumg)


@df edogcsu scatter(:ρ,:ρa,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-0.1,4.1),yticks=0:0.5:4,
    xticks=0:0.5:4,xlim=(-0.1,4.1),xlabel="ρ",ylabel="ρa")
@df edogcsu scatter(:cs,:onoff,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-10,10),
    xlabel="CenterSurround",ylabel="OnOff")
@df edogcsu corrplot(cols([:diameter,:cs,:onap,:onwp]),bin=50,leg=false,size=(850,650),grid=false)


@df gaborcsu scatter(:p,:cyc,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-0.1,4.1),yticks=0:0.5:4,
    xlabel="p",ylabel="cyc")
@df gaborcsu scatter(:onoff,:cyc,leg=false,size=(600,450),ms=3,msw=0,ma=0.7,ylim=(-0.1,4.1),yticks=0:0.5:4,
    xlabel="OnOff",ylabel="cyc")
@df gaborcsu corrplot(cols([:diameter,:onoff,:onap,:onwp]),bin=50,leg=false,size=(850,650),grid=false)


@df csumg scatter(:bicgabor,:bicedog,ratio=1,leg=false,ms=3,msw=0,ma=0.7,grid=false)
lim=1.55e5
plot!([0,lim],[0,lim],color=:black,lw=1,xlim=[0,lim],ylim=[0,lim],xlabel="BIC (Gabor)",ylabel="BIC (eDoG)",xticks=false,yticks=false)
foreach(ext->savefig(joinpath(resultroot,"csu_bic$ext")),figfmt)

@df csumg scatter(:fvugabor,:fvuedog,ratio=:equal,leg=false,ms=3,msw=0,ma=0.7,grid=false)
plot!([0,1],[0,1],color=:black,lw=1,xlim=[0,1],xlabel="FVU (Gabor)",ylabel="FVU (eDoG)")
foreach(ext->savefig(joinpath(resultroot,"csu_fvu$ext")),figfmt)

@df csumg scatter(:corgabor,:coredog,ratio=:equal,leg=false,ms=3,msw=0,ma=0.7,grid=false)
plot!([0,1],[0,1],color=:black,lw=1,xlim=[0,1],xlabel="R (Gabor)",ylabel="R (eDoG)")
foreach(ext->savefig(joinpath(resultroot,"csu_cor$ext")),figfmt)



## Spatial Pattern Classification
ccodecm = Dict('A'=>ColorMaps["dkl_mcclumiso"],'L'=>ColorMaps["lms_mccliso"],'M'=>ColorMaps["lms_mccmiso"],'S'=>ColorMaps["lms_mccsiso"])
maskimg(i;imgsize=(51,51),masktype="Gaussian") = float32.(alphamask(imresize(i,imgsize);masktype))
maskimg(i,c;imgsize=(51,51),masktype="Gaussian") = float32.(alphamask(imresize(i,imgsize);color=ccodecm[c].colors,masktype))

edogimg = map(maskimg,edogcsu.img)
gaborimg = map(maskimg,gaborcsu.img)
edogimgc = map((i,c)->maskimg(i,c),edogcsu.img,edogcsu.rc)
gaborimgc = map((i,c)->maskimg(i,c),gaborcsu.img,gaborcsu.rc)


# eDoG
features = [:onoffsign,:onap,:onwp]
Y = hcat(map(i->edogcsu[:,i],features)...)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
    Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Makie.scatter(Y2[2,:], Y2[1,:],markersize=5)
    # Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1000, 1000))
end

Y2 = umap(Y', 2, n_neighbors=50, min_dist=2, n_epochs=200,metric=Euclidean())
scatter(Y2[1,:], Y2[2,:],ms=3,msw=0,ma=0.8,frame=:none,leg=false,ratio=1,size=(600,600))

p = GLMakie.scatter(Y2[1,:], Y2[2,:], marker=edogimg, ms=30)
save(joinpath(resultroot,"sta_dog_shape_umap.png"),p)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=dogimgc,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_dog_shape_umap_c.png"),p)

Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=d2clu(Y2,r=2.5),frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
dogcell.umapclu = d2clu(Y2,r=2.5)


save("UMAP_sta.svg",plotunitpositionimage(Y2',edogimg))


cr = kmeans(Y',4)
clu = assignments(cr)

YD = pairwise(Euclidean(),Y,dims=1)
hc = hclust(YD,linkage=:ward,branchorder=:optimal)
plot(hc,xticks=false)
clu = cutree(hc,k=4)
clu = replace(clu,3=>"eDoG-Off",4=>"eDoG-On",1=>"Off",2=>"On")

scatter(Y2[1,:], Y2[2,:],msw=0,ms=2,ma=0.8,group=clu,frame=:none,leg=:inline,ratio=1,size=(600,600),legendfontsize=12)
foreach(ext->savefig(joinpath(resultroot,"sta_edogcsu_pattern_umap_clu$ext")),figfmt)

edogcsu.umap = map(i->Tuple(Y2[:,i]),1:size(Y2,2))
edogcsu.pattern = clu
jldsave(joinpath(resultroot,"edogcsu.jld2");edogcsu)


# Gabor
features = [:onoff,:onap,:onwp]
Y = hcat(map(i->gaborcsu[:,i],features)...)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

@manipulate for i in 100:500, n in 5:50, d in 0.001:0.001:3
    Y2 = umap(permutedims(Y), 2, n_neighbors=n, min_dist=d, n_epochs=i)
    Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1000, 1000))
end

Y2 = umap(Y', 2, n_neighbors=40, min_dist=2, n_epochs=200,metric=Euclidean())
scatter(Y2[1,:], Y2[2,:],ms=3,msw=0,ma=0.8,frame=:none,leg=false,ratio=1,size=(600,600))

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimg,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap.png"),p)

p = Makie.scatter(Y2[1,:], Y2[2,:], marker=gaborimgc,markersize=30,scale_plot = false,show_axis = false,resolution = (1200, 1200))
save(joinpath(resultroot,"sta_gabor_shape_umap_c.png"),p)

Plots.scatter(Y2[1,:], Y2[2,:], marker=(3,3,:circle,stroke(0)),group=d2clu(Y2,r=1.2),frame=:none,leg=:inline,aspect_ratio=1,size=(600,600))
savefig(joinpath(resultroot,"sta_gabor_shape_umap_clu.png"))
gaborcell.umapclu = d2clu(Y2,r=1.2)

save("UMAP_sta.svg",plotunitpositionimage(Y2',gaborimg))


cr = kmeans(Y[:,1:2]',5)
clu = assignments(cr)

YD = pairwise(Euclidean(),Y[:,1:2],dims=1)
hc = hclust(YD,linkage=:ward,branchorder=:optimal)
plot(hc,xticks=false)
clu = cutree(hc,k=5)
clu = replace(clu,3=>"Gabor-Even-Off",4=>"Gabor-Even-On",5=>"Gabor-Odd",1=>"On",2=>"Off")

scatter(Y2[1,:], Y2[2,:],msw=0,ms=2,ma=0.8,group=clu,frame=:none,leg=:inline,ratio=1,size=(600,600),legendfontsize=12)
foreach(ext->savefig(joinpath(resultroot,"sta_gaborcsu_pattern_umap_clu$ext")),figfmt)

gaborcsu.umap = map(i->Tuple(Y2[:,i]),1:size(Y2,2))
gaborcsu.pattern = clu
jldsave(joinpath(resultroot,"gaborcsu.jld2");gaborcsu)


## Model Selection and Goodness
transform!(csumg,[:bicgabor,:bicedog]=>((i,j)->i.<j)=>:gaborbetter)
plotdepthhist(csumg;g=:gaborbetter,layer,ylabel="Number of Spatial Pattern")

smcsu = sort!(append!(gaborcsu[csumg.gaborbetter,:],edogcsu[.!csumg.gaborbetter,:],cols=:union),:rc)
transform!(smcsu,:pattern=>ByRow(i->!startswith(i,'O'))=>:so,:cor=>(i->i.>0.6)=>:goodfit)
sgmcsu = subset(smcsu,:goodfit)

plotdepthhist(smcsu;g=:goodfit,layer,ylabel="Number of Spatial Pattern")
plotdepthhist(sgmcsu;g=:rc,layer,ylabel="Number of Spatial Pattern")
plotdepthhist(sgmcsu;g=:pattern,layer,ylabel="Number of Spatial Pattern")
plotdepthhist(sgmcsu;g=:so,layer,ylabel="Number of Spatial Pattern")

foreach(siteid->plotdepthhist(sgmcsu;dir=resultroot,siteid,g=:pattern,layer,ylabel="Number of Spatial Pattern"),levels(sgmcsu.siteid))
foreach(siteid->plotdepthhist(sgmcsu;dir=resultroot,siteid,g=:sop,layer,ylabel="Number of Spatial Pattern"),levels(sgmcsu.siteid))

p=select(sgmcsu,[:pattern,:siteid,:layer]) |>
[@vlplot(:bar,y={"pattern"},x={"count()"});
@vlplot(:bar,y={"siteid"},x={"count()"},color={"pattern:n",scale={scheme=:category20}});
@vlplot(:bar,y={"layer"},x={"count()"},color={"pattern:n",scale={scheme=:category20}})]


# vdogcell |> [@vlplot(:bar,x={"r",bin={step=0.05,extent=[0,1]},title="r"},y={"count()",title="Number of Linear Filter"});
#                 @vlplot(mark={:line,size=2},x=:r,transform=[{sort=[{field=:r}],window=[{op="count",as="cum"}],frame=[nothing,0]}],
#                 y={"cum",title="Cumulated Number of Linear Filter"})]



plotmodelcontour = (unit;siteid="AG1_V1_ODL1",cg=cgrad(:matter),dir=nothing,figfmt=[".png"],title="$(siteid)_model_contour")->begin
    cunit = sort!(filter(r->r.siteid==siteid,unit),:aligndepth,rev=true)
    sizedeg = cunit.sizedeg[1]
    stimx = [0,1,1,0,0] * sizedeg[1]
    stimy = [1,1,0,0,1] * sizedeg[2]
    t = 0:0.02:2π
    p = plot(stimx,stimy;lw=1,color=:black,frame=:grid,ratio=1,tickdir=:out,xlabel="X (deg)",ylabel="Y (deg)",leg=false)
    for r in eachrow(cunit)
        # roi is derived from sta image with pixel origin at topleft corner
        xys = map(i->r.fit.cfun(i,r.fit.param).+(r.roi.centerdeg[2], sizedeg[2]-r.roi.centerdeg[1]),t)
        plot!(p,xys,la=0.8,color=cg[r.aligndepth])
    end
    isnothing(dir) && return p
    mkpath(dir);foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotmodelcontour(sgmcsu,siteid="AG1_V1_ODL1")
foreach(siteid->plotmodelcontour(sgmcsu;siteid,dir=joinpath(resultroot,"RF")),levels(srmcsu.siteid))


plotdepthscatter = (unit;x=:diameter,xfun=first,g=:pattern,dir=nothing,siteid="all",figfmt=[".png"],layer=nothing,xlabel="",
    wfun=x->isempty(x) ? 0 : median(x),palette=:tab20,leg=:best) -> begin
    cunit = siteid=="all" ? unit : filter(r->r.siteid==siteid,unit)
    xx = xfun.(cunit[!,x])
    title="$(siteid)_CorticalDepth_$g"
    xmin,xmax = extrema(xx)
    p = Plots.scatter(xx,cunit.aligndepth;group=cunit[!,g],yflip=true,grid=false,legendfontsize=6,
        ms=3,msw=0,ma=0.7,ylims=(0,1),size=(350,500),tickor=:out,ylabel="Cortical Depth",xlabel,palette,leg)

    n,y = unitdensity(cunit.aligndepth;w=xx,wfun,spacerange=(0,1),bw=0.05,step=0.05)
    Plots.plot!(p,n,y,color=:gray40,label="Average")
    if !isnothing(layer)
        ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),Plots.text(k,7,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[1] for l in values(layer)];linecolor=:gray25,label="layer",lw=1,ann)
    end
    isnothing(dir) ? p : foreach(ext->savefig(joinpath(dir,"$title$ext")),figfmt)
end

plotdepthscatter(sgmcsu;x=:diameter,g=:pattern,layer,xlabel="Diameter (deg)")
plotdepthscatter(sgmcsu;x=:diameter,g=:rc,layer,xlabel="Diameter (deg)")
plotdepthscatter(sgmcsu;x=:rd,g=:pattern,layer,xlabel="Roundness")
plotdepthscatter(sgmcsu;x=:rd,g=:rc,layer,xlabel="Roundness")
plotdepthscatter(sgmcsu;x=:dc,g=:pattern,layer,xlabel="Delay (ms)")
plotdepthscatter(sgmcsu;x=:dc,g=:rc,layer,xlabel="Delay (ms)")


## Relation of Spatial Pattern between Spectral Channels
function plotssmpair(ssm,mp;w=350,h=350,dlim=nothing,mulim=nothing)
    pssm = filter(r->contains(r.MRTYPE,mp.first) && contains(r.MRTYPE,mp.second),ssm)

    p = Plots.plot(layout=(3,1),leg=false,size=(w,3h))
    pn = :diameter
    x = getindex.(pssm[!,pn],mp.first)
    y = getindex.(pssm[!,pn],mp.second)
    lim = isnothing(dlim) ? max(maximum(x),maximum(y)) : dlim
    Plots.scatter!(p[1],x,y;xlabel="$(mp.first) (deg)",ylabel="$(mp.second) (deg)",title="Diameter",
    ms=3,msw=0,ma=0.7,ratio=1)
    Plots.plot!(p[1],[0,lim],[0,lim],color=:gray30)

    x = rad2deg.(getindex.(pssm[!,:θ],mp.first))
    y = rad2deg.(getindex.(pssm[!,:θ],mp.second))
    xticks=yticks=0:45:180
    Plots.scatter!(p[2],x,y;xlabel="$(mp.first) (deg)",ylabel="$(mp.second) (deg)",title="Orientation",
    xticks,yticks,ms=3,msw=0,ma=0.7,ratio=1)
    Plots.plot!(p[2],[0,180],[0,180],color=:gray30)

    x = getindex.(pssm[!,:μx],mp.second) .- getindex.(pssm[!,:μx],mp.first)
    y = getindex.(pssm[!,:μy],mp.second) .- getindex.(pssm[!,:μy],mp.first)
    lim = isnothing(mulim) ? maximum(abs.([x;y])) : mulim
    xlims=ylims=[-lim,lim]
    Plots.scatter!(p[3],x,y;xlabel="X$(mp.second) - X$(mp.first) (deg)",ylabel="Y$(mp.second) - Y$(mp.first) (deg)",
    title="Center Displacement",ms=3,msw=0,ma=0.7,ratio=1,xlims,ylims)
    Plots.vline!(p[3],[0],color=:gray30)
    Plots.hline!(p[3],[0],color=:gray30)
end


ssmcsu = rightjoin(combine(groupby(sgmcsu,:id),:rc=>join=>:MRTYPE,
        [:rc,:diameter]=>((c,v)->OrderedDict((c.=>v)...))=>:diameter,
        [:rc,:θ]=>((c,v)->OrderedDict((c.=>v)...))=>:θ,
        [:rc,:μx]=>((c,v)->OrderedDict((c.=>v)...))=>:μx,
        [:rc,:μy]=>((c,v)->OrderedDict((c.=>v)...))=>:μy,
        [:rc,:W]=>((c,v)->OrderedDict((c.=>v)...))=>:W,
        :μx=>mean=>:ssmμx,
        :μy=>mean=>:ssmμy,
        [:μx,:μy,:diameter]=>((x,y,d)->maximum([abs.(x.-mean(x)).+0.6d;abs.(y.-mean(y)).+0.6d]))=>:ssmradius,
        [:rc,:fit]=>((c,v)->OrderedDict((c.=>v)...))=>:ssm,
        [:rc,:pattern]=>((c,p)->join(c.*'_'.*p,", "))=>:ss,
        [:rc,:so]=>((c,o)->join(c.*'_'.*replace(o,true=>"O",false=>"N"),", "))=>:sso),
        stacsu,on=:id)
transform!(ssmcsu,:MRTYPE=>ByRow(i->!ismissing(i))=>:mresponsive)
ssgmcsu = dropmissing(ssmcsu,:MRTYPE)
transform!(ssgmcsu,:MRTYPE=>ByRow(i->replace(i,r"A([L,M,S]+)"=>s"\1"))=>:MCTYPE,
        :ss=>ByRow(i->replace(i,r"A_\S+, ([L,M,S]_.+)"=>s"\1"))=>:css,
        :sso=>ByRow(i->replace(i,r"A_\S+, ([L,M,S]_.+)"=>s"\1"))=>:csso,
        :W=>ByRow(d->begin
            ub = reduce(max,values(d))
            OrderedDict(k=>v/ub for (k,v) in d)
        end)=>:nW,
        :W=>ByRow(d->begin
            L = haskey(d,'L') ? d['L'] : 0
            M = haskey(d,'M') ? d['M'] : 0
            S = haskey(d,'S') ? d['S'] : 0
            L==M==S==0 && return [missing,missing]
            Σ = L+M+S
            [L/Σ,M/Σ]
        end)=>[:wl,:wm])

plotdepthhist(ssmcsu;g=:mresponsive,layer)

plotssmpair(ssgmcsu,'A'=>'L')
plot((plotssmpair(ssgmcsu,i=>j,dlim=6,mulim=0.6) for (i,j) in combinations(['A','L','M','S'],2))...,
            layout=(1,6),size=(6*350,3*350))

# Cone Weight
csgmcsu = dropmissing(ssgmcsu,[:wl,:wm])
plot([0,1],[1,0],ratio=1,leg=false,xlims=(-0.01,1.01),ylims=(-0.01,1.01),color=:black,tickor=:out)
scatter!(csgmcsu.wl,csgmcsu.wm,ms=4,msw=0,color=:dodgerblue,xlabel="L Cone Weight",ylabel="M Cone Weight")

# typecell |>
#     @vlplot(layer=[{mark={:line,size=1,color="black"},data={values=[{x=0,y=1},{x=1,y=0}]},x=:x,y=:y},
#     {mark={:line,size=1,color="black"},data={values=[{x=-1,y=0},{x=0,y=1}]},x=:x,y=:y},
#     {mark={:line,size=1,color="black"},data={values=[{x=-1,y=0},{x=0,y=-1}]},x=:x,y=:y},
#     {mark={:line,size=1,color="black"},data={values=[{x=0,y=-1},{x=1,y=0}]},x=:x,y=:y},
#     {:point,x={:wl,title="L Cone Weight"},y={:wm,title="M Cone Weight"},color={"layer",scale={scheme=:category10}}}])


## Cell Spectral-Spatial Classification
function plotssmimg(ssm;rand=false,type=:img,ni=nothing,maxn=nrow(ssm),col=5,h=130,w=type==:img ? 450 : 130,path=nothing,bg=:gray,
    imgsize=(55,55),titlefontsize=type==:img ? 14 : 6,hidetitle=false,cch=['A','L','M','S'])
    nu = nrow(ssm); n = min(maxn,nu); n==0 && return
    nc = n < col ? n : col
    nr = n <= col ? 1 : ceil(Int,n/col)
    ni = rand ? sample(1:nu,n,replace=false) : isnothing(ni) ? (1:n) : ni
    p=plot(;layout=(nr,nc),size=(nc*w,nr*h),frame=:none,leg=false,titlefontsize)
    for i in eachindex(ni)
        title = hidetitle ? "" : "$(ssm.id[ni[i]])"
        ssmim = ssm.ssmim[ni[i]]
        cs = intersect(cch,keys(ssmim))
        if type == :img
            t = hcat((maskimg(ssmim[c].img,c;imgsize,masktype="None") for c in cs)...)
            plot!(p[i],t;title,bginside=bg,frame=:none)
        elseif type == :contour
            x = 1:imgsize[2]; y = 1:imgsize[1]; xlims=(x[begin],x[end]);ylims=(y[begin],y[end])
            for c in cs
                t = imresize(ssmim[c].img,imgsize)
                lb,ub = extrema(t)
                lim = max(abs(lb),abs(ub))
                levels=[0.1lb,0.1ub]
                colors = get(cgrad(ccodecm[c].colors),[lb,ub],(-lim,lim))
                for j in 1:2
                    for l in Contour.lines(Contour.contour(y,x,t,levels[j]))
                        ys,xs = Contour.coordinates(l)
                        plot!(p[i],xs,ys;color=colors[j],yflip=true,ratio=1,
                            xlims,ylims,lw=2,title,bginside=bg,frame=:none)
                    end
                end
            end
        end
    end
    isnothing(path) && return p
    mkpath(dirname(path));savefig(path)
end


transform!(ssgmcsu,[:ssm,:ssmμx,:ssmμy,:ssmradius]=>ByRow((m,μx,μy,r)->OrderedDict(k=>fitpredict(v,μx,μy,r) for (k,v) in m))=>:ssmim)
p=select(ssgmcsu,[:ss,:sso,:css,:csso]) |> @vlplot(:bar,y={"count()"},x={"ss"})

foreach(t->plotssmimg(filter(r->r.sso == t,ssgmcsu);maxn=nrow(ssgmcsu),col=5,hidetitle=true,type=:img,path=joinpath(resultroot,"RF","$(t).img.png")),
        levels(ssgmcsu.sso))

plotssmimg(filter(r->r.sso == "S_O",ssgmcsu);maxn=1,col=1,hidetitle=true,rand=true)
plotssmimg(filter(r->r.sso == "A_O, L_N",ssgmcsu);maxn=3,col=1,hidetitle=true,type=:contour)
plotssmimg(filter(r->r.sso == "A_O, L_N" || r.sso == "A_N, L_O",ssgmcsu);maxn=20,col=1,hidetitle=true,type=:contour)
plotssmimg(filter(r->r.sso == "L_O, M_O, S_N",ssgmcsu);maxn=5,col=1,hidetitle=true,type=:img)
plotssmimg(filter(r->r.sso == "A_O, L_O, M_N, S_N",ssgmcsu);maxn=1,col=1,hidetitle=true,type=:img)

plotssmimg(filter(r->r.csso == "A_N",ssgmcsu);maxn=5,col=1,hidetitle=true,rand=true,cch=['A'])



function ssmcombine(ssmim)
    cs = collect(keys(ssmim))
    n = length(cs)
    n == 1 && return missing
    mp=OrderedDict()
    for (i,j) in combinations(cs,2)
        k = i*j
        imgi,maski=ssmim[i]
        imgj,maskj=ssmim[j]
        ar = log2(count(maski)/count(maskj))
        o = maski.&maskj
        Ω = maski.|maskj
        op = count(o)/count(Ω)
        cop,_,_ = mimetric(imgi[o].*imgj[o])
        ocor = cor(imgi[o],imgj[o])
        owr = log2(sum(abs.(imgi[o]))/sum(abs.(imgj[o])))
        mp[k] = (;op,cop,ocor,owr,ar,oix=imgi[o],oiy=imgj[o])
    end
    n == 2 && return mp
    for (i,j,p) in combinations(cs,3)
        k = i*j*p
        imgi,maski=ssmim[i]
        imgj,maskj=ssmim[j]
        imgp,maskp=ssmim[p]
        ai = count(maski);aj = count(maskj);ap = count(maskp)
        ra = (;Symbol(i)=>ai/(ai+aj+ap),Symbol(j)=>aj/(ai+aj+ap))
        o = maski.&maskj.&maskp
        Ω = maski.|maskj.|maskp
        op = count(o)/count(Ω)
        # cop,_,_ = mimetric(imgi[o].*imgj[o])
        ocor = (;Symbol(i*j)=>cor(imgi[o],imgj[o]),Symbol(i*p)=>cor(imgi[o],imgp[o]))
        opcor = (;Symbol(i*j)=>partialcor(imgi[o],imgj[o],imgp[o]),Symbol(i*p)=>partialcor(imgi[o],imgp[o],imgj[o]))
        wi=sum(abs.(imgi[o]));wj=sum(abs.(imgj[o]));wp=sum(abs.(imgp[o]))
        orw = (;Symbol(i)=>wi/(wi+wj+wp),Symbol(j)=>wj/(wi+wj+wp))
        mp[k] = (;op,ocor,opcor,orw,ra,oix=imgi[o],oiy=imgj[o],oiz=imgp[o])
    end
    mp
end

transform!(ssgmcsu,:ssmim=>ByRow(ssmcombine)=>:ssmcm)


t=as
plotdepthscatter(t;x=:op,g=:MRTYPE,layer,xlabel="Overlap %")
plotdepthscatter(t;x=:cop,g=:MRTYPE,layer,xlabel="Covariant %")
plotdepthscatter(t;x=:ocor,g=:MRTYPE,layer,xlabel="Correlation")
plotdepthscatter(t;x=:owr,g=:MRTYPE,layer,xlabel="Weight Balance")
plotdepthscatter(t;x=:ar,g=:MRTYPE,layer,xlabel="RF Area Balance")

lm = transform!(subset!(dropmissing(ssgmcsu,:ssmcm),:ssmcm=>ByRow(d->haskey(d,"LM"))),:ssmcm=>ByRow(d->d["LM"])=>[:op,:cop,:ocor,:owr,:ar,:oix,:oiy])
ls = transform!(subset!(dropmissing(ssgmcsu,:ssmcm),:ssmcm=>ByRow(d->haskey(d,"LS"))),:ssmcm=>ByRow(d->d["LS"])=>[:op,:cop,:ocor,:owr,:ar,:oix,:oiy])
ms = transform!(subset!(dropmissing(ssgmcsu,:ssmcm),:ssmcm=>ByRow(d->haskey(d,"MS"))),:ssmcm=>ByRow(d->d["MS"])=>[:op,:cop,:ocor,:owr,:ar,:oix,:oiy])
al = transform!(subset!(dropmissing(ssgmcsu,:ssmcm),:ssmcm=>ByRow(d->haskey(d,"AL"))),:ssmcm=>ByRow(d->d["AL"])=>[:op,:cop,:ocor,:owr,:ar,:oix,:oiy])
am = transform!(subset!(dropmissing(ssgmcsu,:ssmcm),:ssmcm=>ByRow(d->haskey(d,"AM"))),:ssmcm=>ByRow(d->d["AM"])=>[:op,:cop,:ocor,:owr,:ar,:oix,:oiy])
as = transform!(subset!(dropmissing(ssgmcsu,:ssmcm),:ssmcm=>ByRow(d->haskey(d,"AS"))),:ssmcm=>ByRow(d->d["AS"])=>[:op,:cop,:ocor,:owr,:ar,:oix,:oiy])

@df as corrplot(cols([:op,:cop,:ocor,:owr,:ar]),bin=50,leg=false,size=(850,650),grid=false)


lms = innerjoin(lm,ls,ms,on=:id,makeunique=true)
scatter(lms.ar_2,lms.ar_1;ratio=1,leg=false,xlabel="MS",ylabel="LS")
scatter(lms.ar_2,lms.ar;ratio=1,leg=false,xlabel="MS",ylabel="LM")
scatter(lms.ar_1,lms.ar;ratio=1,leg=false,xlabel="LS",ylabel="LM")


lms = transform!(subset!(dropmissing(ssgmcsu,:ssmcm),:ssmcm=>ByRow(d->haskey(d,"LMS"))),:ssmcm=>ByRow(d->d["LMS"])=>[:op,:ocor,:opcor,:orw,:ra,:oix,:oiy,:oiz])
histogram(lms.op,bins=50,xlims=(-0.01,1.01),lc=:match,leg=false,xlabel="op")
scatterhist(lms.ra)
scatterhist(lms.opcor,lims=(-1,1),addline=false,frame=:zerolines)

plotdepthscatter(lms;x=:op,g=:MRTYPE,layer,xlabel="Overlap %")
plotdepthscatter(lms;x=:ra,g=:MRTYPE,layer,xlabel="L Area Balance")
plotdepthscatter(lms;x=:ra,xfun=last,g=:MRTYPE,layer,xlabel="M Area Balance")
plotdepthscatter(lms;x=:orw,g=:MRTYPE,layer,xlabel="L Weight Balance")
plotdepthscatter(lms;x=:orw,xfun=last,g=:MRTYPE,layer,xlabel="M Weight Balance")
plotdepthscatter(lms;x=:opcor,g=:MRTYPE,layer,xlabel="LM Partial Correlation")
plotdepthscatter(lms;x=:opcor,xfun=last,g=:MRTYPE,layer,xlabel="LS Partial Correlation")

function scatterhist(x;bins=50,lims=(0,1),addline=true,frame=:axes)
    layout = @layout [a       _ 
                      b{0.85w,0.85h} c]
    xlims=ylims=lims.+(-0.01,0.01)
    p=plot(;layout,leg=false,link=:both,size=(650,630),tickor=:out)
    addline && plot!(p[2,1],[0,1],[1,0];ratio=1,xlims,ylims,color=:black)
    scatter!(p[2,1],x;ms=4,msw=0,color=:dodgerblue,frame)
    histogram!(p[1,1],map(first,x);bins,xlims,frame=:none,lc=:match)
    histogram!(p[2,2],map(last,x);bins,dir=:h,ylims,frame=:none,lc=:match)
    p
end


i=49#25
scatter(lm.oix[i],lm.oiy[i],ms=1,msw=0,leg=false,ratio=1,frame=:origin,xticks=[0],yticks=[0])
plotssmimg(lm[i:i,:],cch=['L','M'],type=:contour,hidetitle=true)







clms = outerjoin(lm,ls,ms,on=:id,makeunique=true)

mapcols!(i->replace!(i,missing=>0),clms)




p=GLMakie.Figure()
ax=GLMakie.Axis3(p[1,1],#xticks=[-1,0,1],yticks=[-1,0,1],zticks=[-1,0,1],
xlabel="LM",ylabel="LS",zlabel="MS",viewmode=:fit,aspect=:data)
GLMakie.scatter!(ax,lms.ar,lms.ar_1,lms.ar_2)


display(p)

function plotrf(unit;posx=unit.posx,posy=unit.posy,resolution=(1000,1000),type=:contour,path=nothing,
    imgsize=(51,51),cch=['A','L','M','S'],bg=:white,rfscale=1)
    x = 1:imgsize[2]; y = 1:imgsize[1]; imgd=max(imgsize...)
    f = GLMakie.Figure(;resolution,backgroundcolor=bg)
    ax = GLMakie.Axis(f[1,1],yreversed=true,autolimitaspect=1,backgroundcolor=bg)
    GLMakie.hidedecorations!(ax)
    GLMakie.hidespines!(ax)
    for i in 1:nrow(unit)
        ssmim = unit.ssmim[i]
        if type == :img
            cs = intersect(cch,keys(ssmim))
            img = mapreduce(c->maskimg(ssmim[c].img,c;imgsize,masktype="Disk"),(i,j)->weighted_color_mean.(0.5,i,j),cs)
            # img = map(c->maskimg(ssmim[c].img,c;imgsize,masktype="Disk"),cs)
            GLMakie.scatter!(ax,posx[i],posy[i],marker=img,markersize=30)
        elseif type == :contour
            for (c,v) in ssmim
                c in cch || continue
                img = imresize(v.img,imgsize)
                lb,ub = extrema(img)
                lim = max(abs(lb),abs(ub))
                levels=[0.1lb,0.1ub]
                colors = get(cgrad(ccodecm[c].colors),[lb,ub],(-lim,lim))
                for j in 1:2
                    for l in Contour.lines(Contour.contour(y,x,img,levels[j]))
                        ys,xs = Contour.coordinates(l)
                        GLMakie.lines!(ax,xs/imgd*rfscale.+posx[i],ys/imgd*rfscale.+posy[i];color=(colors[j],0.7),linewidth=2)
                    end
                end
            end
        end
    end
    isnothing(path) && return f
    mkpath(dirname(path));save(path,f)
end



t=as
features = [:op,:ocor,:owr,:ar]
Y = hcat(map(i->t[:,i],features)...)
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

t=lms
features = [:op,:opcor,:orw,:ra]
Y = hcat(map(i->[first.(t[:,i]);;last.(t[:,i])],features)...)[:,2:end]
foreach(i->Y[:,i]=zscore(Y[:,i]),1:size(Y,2))

myzscore(x;μ=mean(x),σ=std(x,mean=μ))=(x.-μ)/σ
fm = [0.5,0,0,0]
foreach(i->Y[:,i]=myzscore(Y[:,i],μ=fm[i]),1:size(Y,2))


Y2 = umap(Y', 2, n_neighbors=25, min_dist=4, n_epochs=200,metric=Euclidean())
scatter(Y2[1,:], Y2[2,:],ms=3,msw=0,ma=0.8,frame=:none,leg=false,ratio=1,size=(600,600))

p=plotrf(t,posx=Y2[1,:],posy=Y2[2,:],cch=['A','S'],type=:img,rfscale=4)
display(p)

cr = kmeans(Y',2)
clu = assignments(cr)

YD = pairwise(Euclidean(),Y,dims=1)
hc = hclust(YD,linkage=:ward,branchorder=:optimal)
plot(hc,xticks=false)
clu = cutree(hc,k=5)

scatter(Y2[1,:], Y2[2,:],msw=0,ms=3,ma=0.8,group=clu,frame=:none,yflip=true,leg=:inline,ratio=1,size=(600,600),legendfontsize=12)
scatter(Y2[1,:], Y2[2,:],msw=0,ms=4,ma=0.8,color=get(cgrad(:coolwarm),t.ocor,(-1,1)),frame=:none,yflip=true,leg=false,ratio=1,size=(600,600),legendfontsize=12)


transform!(ssgmcsu,[:csso,:ssmcm]=>ByRow(cellclass)=>[:class1,:class2])

function cellclass(csso,ssmcm)
    c1=c2="N"
    if contains(csso,", ")
        c = ['L','M','S']
        ci = contains.(csso,c)
        k = join(c[ci])
        cm=ssmcm[k]
        if length(k)==3
            lm = cm.opcor.LM >= 0 ? "/" : "\\"
            ls = cm.opcor.LS >= 0 ? "/" : "\\"
            c1=join(["3",lm,ls])
            c2 = join(c[ci],lm,ls)
        else
            t = cm.ocor >= 0 ? "/" : "\\"
            c1 = "2$t"
            c2 = join(c[ci],t)
        end
    else
        if contains(csso,"A")
            c1 =c2= "A"
        else
            if contains(csso,"N")
                c1 = "0"
            else
                c1 = "1"
            end
            c2 = "$(csso[1])$c1"
        end        
    end
    [c1,c2]
end

plotdepthhist(ssgmcsu;g=:class2,layer)

dn = :none
plotdepthhist(subset(ssgmcsu,dn);g=:class1,layer)

p=select(subset(ssgmcsu,dn),[:class1,:layer]) |>
[@vlplot(:bar,x={"class1"},y={"count()"});
@vlplot(:bar,y={"layer"},x={"count()"},color={"class1:n",scale={scheme=:category20}})]




limg = imresize(t.ssmim[22]['L'].img,51,51) 
mimg = imresize(t.ssmim[22]['M'].img,51,51) 
simg = zeros(size(limg)) 

lim = maximum(abs.([limg;mimg]))

limg = limg./lim .+ lmsbg[1]
mimg = mimg./lim .+ lmsbg[2]
simg = simg./lim .+ lmsbg[3]



lmsimg=reshape(permutedims(cat(limg,mimg,simg,dims=3),[3,1,2]),3,:)

rgbimg = LMSToRGB[1:3,1:3]*lmsimg
rgbimg ./= maximum(rgbimg,dims=1) 

rgbpx = LMSToRGB[1:3,1:3]*([0.5;0.5;0.5;;-0.5;-0.5;-0.5;;0.5;-0.5;0;;-0.5;0.5;0].+lmsbg)
rgbpx ./= maximum(rgbpx,dims=1)

colorview(RGB,rgbpx)

cimg = reshape(colorview(RGB,rgbimg),51,:)

heatmap(mimg,yflip=true,ratio=1)


lon = RGB(ccodecm['L'].colors[end])
loff = RGB(ccodecm['L'].colors[1])
mon = RGB(ccodecm['M'].colors[end])
moff = RGB(ccodecm['M'].colors[1])

weighted_color_mean(0.5,lon,mon)
weighted_color_mean(0.5,loff,moff)
weighted_color_mean(0.5,loff,moff)
lon+mon
a=(lon+mon)/2

loff+moff
a=(loff+moff)/2

import  Base:vec
vec(a::RGB)=[a.r,a.g,a.b]

lon+moff
a=(lon+moff)/2
a/=maximum(vec(a))

loff+mon
a=(loff+mon)/2
a/=maximum(vec(a))

RGBToLMS,LMSToRGB=load(joinpath(resultroot,"lmsmatrix.jld2"),"RGBToLMS","RGBToLMS")


lmsbg = (RGBToLMS*[0.5,0.5,0.5,1])[1:3]










savefig("test.png")
save("test.png",p)
pyplot()




typecell |> [@vlplot(:bar,y={"site",title="Recording Site"},x={"count()",title="Number of Cells"},color={"celltype:n"});
        @vlplot(:bar,y={"depth",bin={step=200},sort="descending",title="Cortical Depth (μm)"},x={"count()",title="Number of Cells"},color={"celltype:n"});
        @vlplot(:bar,y={"layer",title="Cortical Layer"},x={"count()",title="Number of Cells"},color={"celltype:n"})]









save("UMAP_c_sta.svg",plotunitpositionimage(Y2',rswcell.srimg))

save("UMAP_c_sta_layer.svg",plotunitlayerimage(rswcell.layer,rswcell.imgs,width=2000,markersize=30))

## non-sta but orisf responsive






