using Combinatorics,VegaLite,DataFrames,Clustering,TSne,MultivariateStats,UMAP

## t-SNE and UMAP of responsive sta
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

## UMAP of responsive rf
m=:gabor
rfs = mapreduce(u->mapreduce(f->f.r > 0.6 ? [(;(n=>f[n] for n in (:model,:radius,:param))...)] : [],
append!,skipmissing(dataset["urf"][u][m])),
append!,keys(dataset["urf"]))

mrf = map(f->begin
                i = rfimage(f,f.radius)
                m = maximum(abs.(extrema(i)))
                alphamask(get(cgrad(:coolwarm),i,(-m,m)),radius=0.5,sigma=0.35,masktype="Gaussian")[1]
             end,rfs)
@manipulate for i in eachindex(mrf)
    mrf[i]
end

X = mapreduce(f->f.param,hcat,rfs)
Y = X
foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
@manipulate for i in 200:2000, n in 5:60, d in 0.001:0.001:10
    Y2 = umap(Y, 2, n_neighbors=n, min_dist=d, n_epochs=i)'
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:auto,stroke(0)),frame=:none,leg=false)
end

Y2 = umap(Y, 2, n_neighbors=15, min_dist=0.01, n_epochs=500)'
p = plotunitpositionimage(Y2,mrf)
save("UMAP_rf_$m.svg",p)

## type of responsiveness
CTYPES = sort(join.(sort.(combinations(unique(values(ccode))))))
CTCOLORS = cgrad(:rainbow,length(CTYPES),categorical = true).colors
crtype = Dict(u=>join(sort(dataset["ccode"][dataset["ucresponsive"][u]])) for u in keys(dataset["ulsta"]))

p = DataFrame(ResponseType=collect(values(crtype))) |> @vlplot(:bar,y=:ResponseType,x="count()",width=600,height=400)

## type of model responsiveness
m=:gabor
cmrtype = filter!(p->!isempty(p.second),Dict(u=>join(sort(dataset["ccode"][map(f->!ismissing(f) && f.r > 0.6,dataset["urf"][u][m])])) for u in keys(dataset["urf"])))

p = DataFrame(ModelResponseType=collect(values(cmrtype))) |> @vlplot(:bar,y=:ModelResponseType,x="count()",width=600,height=400)

## rf space of unit
modelnull=Dict(:gabor=>fill(missing,8),:dog=>[])
function rf2unitspace(ps;model=:gabor)
    if model == :gabor
        as = ps[1:8:end]
        maxamp = maximum(skipmissing(as))
        foreach(i->ps[i]/=maxamp,1:8:length(ps))
        ps=replace(ps,missing=>0)
    elseif model == :dog
    end
    ps
end
"map RF models to unit feature space"
function rfunitspace(dataset;model=:gabor)
    corder = sortperm(dataset["ccode"])
    ufs=Dict()
    for u in keys(dataset["urf"])
        ps = vcat(map(f->ismissing(f) || f.r < 0.6 ? modelnull[model] : f.param,dataset["urf"][u][model])[corder]...)
        if !isempty(skipmissing(ps))
            ufs[u] = rf2unitspace(ps,model=model)
        end
    end
    ufs
end


rfspace = rfunitspace(dataset,model=:gabor)
Xg = [cmrtype[u] for u in keys(rfspace)]
X = hcat(values(rfspace)...)
Y = X
foreach(i->Y[i,:]=zscore(Y[i,:]),1:size(Y,1))
@manipulate for i in 1000:6000, p in 2:60
    Y2 = tsne(permutedims(Y), 2, 0, i, p)
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:circle,stroke(0)),group=Xg,frame=:none,leg=:inline)
    # savefig("test.svg")
end

@manipulate for i=200:2000,n in 5:60, d in 0.001:0.001:1000
    Y2 = umap(Y, 2, n_neighbors=n, min_dist=d, n_epochs=i)
    scatter(Y2[:,1], Y2[:,2], marker=(5,5,:circle,stroke(0)),group=Xg,frame=:none,leg=:inline)
end
