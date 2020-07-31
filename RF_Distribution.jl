using Combinatorics,Query,VegaLite,DataFrames,Clustering,MultivariateStats,UMAP,Makie

## UMAP of responsive sta
lsta = mapreduce(u->mapreduce((c,r,exd)->r ? [dataset["ulsta"][u][:,:,exd.d,c]] : [],
append!,1:testn,dataset["ucresponsive"][u],dataset["ulcexd"][u]),
append!,keys(dataset["ulsta"]))

maxsizepx=mapreduce(i->maximum(size(i)),max,lsta)
map!(i->imresize(i,maxsizepx,maxsizepx),lsta,lsta)
@manipulate for i in eachindex(lsta)
    heatmap(lsta[i],yflip=true,aspect_ratio=:equal,frame=:none,color=:coolwarm)
end

mlsta = map(i->begin
                   m = maximum(abs.(extrema(i)))
                   alphamask(get(cgrad(:coolwarm),i,(-m,m)),radius=0.5,sigma=0.35,masktype="Gaussian")[1]
               end,lsta)
@manipulate for i in eachindex(mlsta)
    mlsta[i]
end

# PCA reduction
X = mapreduce(vec,hcat,lsta)
foreach(i->X[i,:]=zscore(X[i,:]),1:size(X,1))
M = fit(PCA,X,maxoutdim=60,pratio=0.9,mean=0)
Y = MultivariateStats.transform(M,X)

# plot t-SNE sta
@manipulate for i in 1000:2000, p in 5:60
    Y2 = tsne(permutedims(Y), 2, 0, i, p)
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:auto,stroke(0)),frame=:none,leg=false)
end

Y2 = tsne(permutedims(Y), 2, 0, 2000, 10)
p = plotunitpositionimage(Y2,mlsta)
save("tSNE_rlSTA.svg",p)

# plot UMAP sta
@manipulate for i in 200:2000, n in 5:60, d in 0.001:0.001:10
    Y2 = umap(Y, 2, n_neighbors=n, min_dist=d, n_epochs=i)'
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:auto,stroke(0)),frame=:none,leg=false)
end

Y2 = umap(Y, 2, n_neighbors=10, min_dist=0.01, n_epochs=500)'
p = plotunitpositionimage(Y2,mlsta)
save("UMAP_rlSTA.svg",p)

## Collect Responsive RF
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

cells = load(joinpath(resultroot,"cells.jld2"),"cells")

function collectrf(indir;cells=DataFrame(id=[]),model=:gabor,rffile="stadataset.jld2")
    df=DataFrame()
    for (root,dirs,files) in walkdir(indir)
        if rffile in files
            dataset = load(joinpath(root,rffile),"dataset")
            siteid = dataset["siteid"]
            id = ["$(siteid)_SU$u" for u in keys(dataset["urf"])]
            rc = Any[map((c,r)->r ? c : missing,dataset["ccode"],dataset["ucresponsive"][u]) for u in keys(dataset["urf"])]
            rf = Any[map(f->ismissing(f) ? missing : (;(n=>f[n] for n in (:model,:radius,:param,:r))...),dataset["urf"][u][model]) for u in keys(dataset["urf"])]
            vui = map(i->!isempty(skipmissing(i)),rc)
            append!(df,DataFrame(id=id[vui],rc=rc[vui],rf=rf[vui]))
        end
    end
    return outerjoin(cells,unique!(df,:id),on=:id)
end

cells = collectrf(resultroot,cells=cells,model=:gabor)
save(joinpath(resultroot,"cells.jld2"),"cells",cells)


rfcells = dropmissing(cells,:rf)

srfcells = dropmissing!(flatten(rfcells,[:rf,:rc]),:rf)

# rfcells = cells |> 
#         @dropna(:rf) |> 
#         @mutate(rc = map((f,c)->!ismissing(f) && f.r > 0.6 ? c : missing,_.rf,get(_.rc)),rf = map(f->!ismissing(f) && f.r > 0.6 ? f : missing,_.rf)) |> 
#         @filter(!isempty(skipmissing(_.rf))) |> DataFrame

# srfcells = [rfcells |> @mapmany(skipmissing(get(_.rf)),{id=_.id,rf=__}) |> DataFrame rfcells |> @mapmany(skipmissing(get(_.rc)),{rc=__}) |> DataFrame]
# srfcells = leftjoin(srfcells,rfcells[!,Not([:rf,:rc])],on=:id)

srfcells.rfi = rfcolorimage.(rfimage.(srfcells.rf,ppd=40))

@manipulate for i in 1:nrow(srfcells)
    srfcells.rfi[i]
end

S = mapreduce(f->f.param,hcat,srfcells.rf)

Y = deepcopy(S)
foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
@manipulate for i in 200:1000, n in 5:60, d in 0.001:0.001:5
    Y2 = umap(Y, 2, n_neighbors=n, min_dist=d, n_epochs=i)'
    scatter(Y2[:,1], Y2[:,2], marker=(3,3,:auto,stroke(0)),frame=:none,leg=false)
end

Y2 = umap(Y, 2, n_neighbors=20, min_dist=0.015, n_epochs=500)'
p = plotunitpositionimage(Y2,srfcells.rfi)
save("UMAP_rf.svg",p)

## RF Spatial Pattern
Ω(p::Missing)=missing
Ω(p) = p < 0.5 ? 4abs(p-0.25) : 4abs(p-0.75)
SP = mapreduce(i->begin
                    σyoverσx = S[5,i]/S[3,i]
                    f4σy = 4S[5,i]*S[7,i]
                    w = Ω(S[8,i])
                    # w = p<0.25 ? 0.5-p : p>0.75 ? 1.5-p : p
                    # w = 0.5-w
                    [σyoverσx, f4σy, w]
                end,hcat,1:size(S,2))


mmsrf = map(i->RGBAf0.(imresize(i,32,32)),msrf)
Makie.scatter(SP[1,:],SP[2,:],SP[3,:],markersize=0.3,marker=mmsrf)

p=Plots.plot(layout=(3,1),size=(400,3*250),legend=false)
Plots.histogram!(p,subplot=1,SP[1,:],nbins=30,xlabel="σy/σx",ylabel="Number of Cells")
Plots.histogram!(p,subplot=2,SP[2,:],nbins=28,xlabel="4σyf",ylabel="Number of Cells")
Plots.histogram!(p,subplot=3,SP[3,:],nbins=45,xlabel="ω",ylabel="Number of Cells")
p

Y = deepcopy(SP)
foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
Y2 = umap(Y, 2, n_neighbors=25, min_dist=0.02, n_epochs=500)'
p = plotunitpositionimage(Y2,srfcells.rfimage)
save("UMAP_rf.svg",p)

## type of responsiveness
CTYPES = sort(join.(sort.(combinations(['A','L','M','S']))))
CTCOLORS = cgrad(:rainbow,length(CTYPES),categorical = true).colors.colors

rdf = DataFrame(ResponseType=map(ct->join(skipmissing(ct)),crtype))
p = rdf |> @vlplot(:bar,y=:ResponseType,x="count()",width=600,height=400)

mrdf = DataFrame(ModelResponseType=filter!(i->!isempty(i),map(ct->join(skipmissing(ct)),cmrtype)))
p = mrdf |> @vlplot(:bar,y=:ModelResponseType,x="count()",width=600,height=400)

crdf = DataFrame(ConeResponseType=filter(i->!isempty(i),map(t->match(r"A?([L,M,S]*)",t).captures[1],mrdf.ModelResponseType)))
p = crdf |> @vlplot(:bar,y=:ConeResponseType,x="count()",width=600,height=400)
## ALMS space
X = map(fs->vcat(map(f->ismissing(f) || f.r <= 0.6 ? fill(missing,8) : f.param,fs)...),rfs)

Y = replace(hcat(filter(i->!isempty(skipmissing(i)),X)...),missing=>0)
Yg = mrdf.ModelResponseType
foreach(j->begin
    ai = 1:8:size(Y,1)
    maxamp = maximum(Y[ai,j])
    foreach(i->Y[i,j]/=maxamp,ai)
end,1:size(Y,2))

foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
@manipulate for i in 200:1000, n in 5:60, d in 0.001:0.001:5
    Y2 = umap(Y, 2, n_neighbors=n, min_dist=d, n_epochs=i)'
    Plots.scatter(Y2[:,1], Y2[:,2], marker=(3,3,:circle,stroke(0)),group=Yg,frame=:none,leg=:inline)
end

## LMS space
LMSs = hcat(filter(i->!isempty(skipmissing(i)),map(fs->vcat(map(f->ismissing(f) || f.r <= 0.6 ? fill(missing,8) : f.param,fs[2:end])...),rfs))...)

σx = LMSs[3:8:end,:]
p = Makie.scatter(σx[1,:],σx[2,:],σx[3,:],markersize=0.01,axis=(names=(axisnames=("L","M","S"),),),color=:grey20)

plotσx=(σx;title="4σₓ") -> begin
    lmi = [!ismissing(σx[1,i]) && !ismissing(σx[2,i]) for i in 1:size(σx,2)]
    lm = σx[[1,2],lmi]

    lsi = [!ismissing(σx[1,i]) && !ismissing(σx[3,i]) for i in 1:size(σx,2)]
    ls = σx[[1,3],lsi]

    msi = [!ismissing(σx[2,i]) && !ismissing(σx[3,i]) for i in 1:size(σx,2)]
    ms = σx[[2,3],msi]

    lim = maximum(skipmissing(σx))
    p = Plots.plot(layout=(1,3),leg=false,xlims=[0,lim],ylims=[0,lim],title=title)
    Plots.plot!(p,subplot=1,x->x,0,lim)
    Plots.scatter!(p,subplot=1,lm[1,:],lm[2,:],aspect_ratio=1,xlabel="L",ylabel="M")
    Plots.plot!(p,subplot=2,x->x,0,lim)
    Plots.scatter!(p,subplot=2,ls[1,:],ls[2,:],aspect_ratio=1,xlabel="L",ylabel="S")
    Plots.plot!(p,subplot=3,x->x,0,lim)
    Plots.scatter!(p,subplot=3,ms[1,:],ms[2,:],aspect_ratio=1,xlabel="M",ylabel="S")
    p
end

plotσx(4σx)

σyoverσx = mapreduce(i->map(j->LMSs[8j+5,i]/LMSs[8j+3,i],0:2),hcat,1:size(LMSs,2))
p = Makie.scatter(σyoverσx[1,:],σyoverσx[2,:],σyoverσx[3,:],markersize=0.1,axis=(names=(axisnames=("L","M","S"),),),color=:grey20)

plotσx(σyoverσx,title="σy/σx")

f4σy = mapreduce(i->map(j->4LMSs[8j+5,i]*LMSs[8j+7,i],0:2),hcat,1:size(LMSs,2))
p = Makie.scatter(f4σy[1,:],f4σy[2,:],f4σy[3,:],markersize=0.1,axis=(names=(axisnames=("L","M","S"),),),color=:grey20)

plotσx(f4σy,title="f4σy")


ω = mapreduce(i->map(j->Ω(LMSs[8j+8,i]),0:2),hcat,1:size(LMSs,2))


plotσx(ω,title="ω")


θ = LMSs[6:8:end,:]

plotσx(180θ/π,title="θ")


dlm = map(i->sqrt(sum((LMSs[[2,4],i] .- LMSs[[2,4].+8,i]).^2)),1:size(LMSs,2))
Plots.histogram(collect(skipmissing(dlm)),nbins=30,xlabel="LM Distance",leg=false)

dls = map(i->sqrt(sum((LMSs[[2,4],i] .- LMSs[[2,4].+16,i]).^2)),1:size(LMSs,2))
Plots.histogram(collect(skipmissing(dls)),nbins=30,xlabel="LS Distance",leg=false)

dms = map(i->sqrt(sum((LMSs[[2,4].+8,i] .- LMSs[[2,4].+16,i]).^2)),1:size(LMSs,2))
Plots.histogram(collect(skipmissing(dms)),nbins=30,xlabel="MS Distance",leg=false)

savefig("test.svg")








Y = replace(hcat(filter(i->!isempty(skipmissing(i[9:end])),X)...)[9:end,:],missing=>0)
Yg = crdf.ConeResponseType

Z = Y[8:8:end,:]

Plots.scatter(360(Zlm[1,:].-Zlm[2,:]),ones(20),nbins=30,projection=:polar)

## Cone Weights
W = Y[1:8:end,:]
W ./= sum(W,dims=1)

p = Makie.scatter(W[1,:],W[2,:],W[3,:],markersize=0.02,axis=(names=(axisnames=("L","M","S"),),),color=:crimson)
Makie.save("test.png",p)


foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
@manipulate for i in 200:1000, n in 5:60, d in 0.001:0.001:5
    Y2 = umap(Y, 2, n_neighbors=n, min_dist=d, n_epochs=i)'
    Plots.scatter(Y2[:,1], Y2[:,2], marker=(3,3,:circle,stroke(0)),group=Yg,frame=:none,leg=:inline)
end
