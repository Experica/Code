using NeuroAnalysis,Statistics,Plots,Mmap,Images,StatsBase,Combinatorics,LightGraphs,GraphPlot,GraphRecipes,Interact,Dierckx,PyCall


# Prepare Dataset
subject="AE9";recordsession="u004";test="004"
sessionid = join([subject,recordsession],"_")
fileid = join([subject,recordsession,test],"_")

dataroot = "../Data"
datasetroot = "../DataExport"
resultroot = "../Result"
resultdir = joinpath(resultroot,fileid)
isdir(resultdir) || mkdir(resultdir)
lfregex = Regex("^$(uppercase(fileid))[A-Za-z0-9_]*.imec.lf.bin")
lffile = joinpath(dataroot,filter(i->occursin(lfregex,i),readdir(dataroot))[1])
dataset = prepare(joinpath(datasetroot,"$fileid.mat"))

nsavech=dataset["lf"]["meta"]["nSavedChans"]
nsample=dataset["lf"]["meta"]["nFileSamp"]
fs=dataset["lf"]["meta"]["fs"]
nch = dataset["lf"]["meta"]["snsApLfSy"][2]
refmask = refchmaskim(dataset["ap"]["meta"])
mmlfp = Mmap.mmap(lffile,Matrix{Int16},(nsavech,nsample),0)

## Condition Test
ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
spike = dataset["spike"]
eval.([:($(Symbol(k))=spike[$k]) for k in keys(spike)])
condon = ex["CondTest"]["CondOn"]
condoff = ex["CondTest"]["CondOff"]
minconddur=minimum(condoff-condon)
histogram(condoff-condon,bins=20,title="Condition Duration(Set to $conddur ms)")
β

# Laminar CSD
# epochs for each condition test
epochext = 0.0
epochdur = 0.15
epoch = [-epochext epochdur+epochext]
epochs = condon.+epoch
# all epoch LFP segments
ys = subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"])
# all epoch LFP segments for each column of probe
cys=getepochlfpcol(ys,refmask)

# Column 1
mcys=dropdims(mean(cys[1],dims=3),dims=3)
plotanalog(mcys,fs=fs,xext=epochext)
foreach(i->savefig(joinpath(resultdir,"column1_lfp$i")),[".png",".svg"])

mccsd = dropdims(mean(csd(cys[1]),dims=3),dims=3)
plotanalog(imfilter(mccsd,Kernel.gaussian((1,1))),fs=fs,xext=epochext)
foreach(i->savefig(joinpath(resultdir,"column1_csd$i")),[".png",".svg"])

# Column 2
mcys=dropdims(mean(cys[2],dims=3),dims=3)
plotanalog(mcys,fs=fs,xext=epochext)
foreach(i->savefig(joinpath(resultdir,"column2_lfp$i")),[".png",".svg"])

mccsd = dropdims(mean(csd(cys[2]),dims=3),dims=3)
plotanalog(imfilter(mccsd,Kernel.gaussian((1,1))),fs=fs,xext=epochext)
foreach(i->savefig(joinpath(resultdir,"column2_csd$i")),[".png",".svg"])

# Column Mean Normalized CSD
pys = cat(cys...,dims=3)
pcsd = csd(pys)
baseline = pcsd[:,epoch2samplerange([0 0.015],fs),:] # 0-15ms csd as baseline
ncsd=pcsd.-mean(baseline,dims=2) # csd relative to baseline
mncsd = dropdims(mean(ncsd,dims=3),dims=3)
mncsd[1,:].=0;mncsd[end,:].=0
mncsd = imfilter(mncsd,Kernel.gaussian((2,1)))
h=20;depths = h*(1:size(mncsd,1))
plotanalog(mncsd,fs=fs,xext=epochext,color=:RdBu)
foreach(i->savefig(joinpath(resultdir,"column_mean_normalized_csd2$i")),[".png",".svg"])
save(joinpath(resultdir,"csd.jld2"),"csd",mncsd,"depth",depths)









using DSP
ps=mt_pgram(mcys[10,:],fs=fs)
plot(ps.freq[1:20],ps.power[1:20])

# epochs for each condition test
epochext = 0.2
epochdur = 0.0
epoch = [-epochext epochdur+epochext]
epoch = [-epochext 0 ]
epochs = condon[2:2:end].+epoch
# all epoch LFP segments
ys = subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"])
# all epoch LFP segments for each column of probe
cys=getepochlfpcol(ys,refmask)



function powerspectrum(x,fs;maxfreq=100)
    nd=ndims(x)
    if nd==3
        n=size(x,3)
        ep,freqs = powerspectrum(x[:,:,1],fs)
        p = Array{Float64}(undef,size(ep)...,n)
        p[:,:,1]=ep
        for i in 2:n
            p[:,:,i],freqs = powerspectrum(x[:,:,i],fs)
        end
    else
        ps = mt_pgram(x[1,:],fs=fs)
        freqs = filter(i->i<=maxfreq,ps.freq)
        nfreq = length(freqs);nr = size(x,1)
        p = Matrix{Float64}(undef,nr,nfreq)
        p[1,:]=ps.power[1:nfreq]
        for i in 2:nr
            ps = mt_pgram(x[i,:],fs=fs)
            p[i,:]=ps.power[1:nfreq]
        end
    end
    return p,freqs
end

ps,freq = powerspectrum(cys[1],fs)

mps = dropdims(mean(ps,dims=3),dims=3)
plotanalog(ps)


bps,freq = powerspectrum(cys[1],fs)

bmps = dropdims(mean(bps,dims=3),dims=3)
plotanalog(bmps)


zps = ps.-bps
zps = zps./std(bps,dims=3)
zmps = dropdims(mean(zps,dims=3),dims=3)
plotanalog(zmps)






# Unit Spike Trian
epochext = 0.1
@manipulate for u in 1:length(unitspike)
ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
plotspiketrain(ys,timeline=[0,minconddur],title="Unit $u")
end

# Unit Epoch Depth PSTH
epochext = 0.0
epochdur = 0.15
epoch = [-epochext epochdur+epochext]
epochs = condon.+epoch
bw = 0.002
psthbins = -epochext:bw:epochdur+epochext

baseindex = epoch2samplerange([0 0.015],1/bw) # 0-15ms of psth as baseline
unitpsth = map(ust->psth(subrv(ust,epochs,isminzero=true,shift=epochext)[1],psthbins,israte=true,normfun=x->x.-mean(x[baseindex])),unitspike)
depthpsth,x,depths = spacepsth(unitpsth,unitposition,spacebinedges=0:20:depths[end])
depthpsth = imfilter(depthpsth,Kernel.gaussian((1,1)))
plotpsth(depthpsth,x,depths)
foreach(i->savefig(joinpath(resultdir,"depthpsth$i")),[".png",".svg"])
save(joinpath(resultdir,"depthpsth.jld2"),"depthpsth",depthpsth,"depth",depths)










# Unit Position
Layers = load(joinpath(resultroot,"$sessionid.jld2"),"layers")
Layers=nothing
plotunitposition(spike,layers=Layers)
foreach(i->savefig(joinpath(resultdir,"unitposition$i")),[".png",".svg"])
save(joinpath(resultdir,"spike.jld2"),"spike",spike)


# Binary Spike
bepochext = -0.4
bepochdur = 2.2
bepoch = [-bepochext bepochdur+bepochext]
bepochs = ccondon.+bepoch
spikebins = -bepochext:0.001:bepochdur+bepochext

unitbspike = map(ust->float.(histmatrix(subrv(ust,bepochs,isminzero=true,shift=bepochext)[1],spikebins)[1])',unitspike)

# Unit Correlogram
lag=50
ccgs,x,ccgis,projs,eunits,iunits = circuitestimate(unitbspike[unitgood],lag=lag)

@manipulate for i in 1:length(ccgs)
    title = "Correlogram between Unit $(ccgis[i][1]) and $(ccgis[i][2])"
    bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=:x,xtick=[-lag,0,lag])
end

for i in 1:length(ccgs)
    title = "Correlogram between Unit $(ccgis[i][1]) and $(ccgis[i][2])"
    bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=:x,xtick=[-lag,0,lag])
    foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
end
save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits)

# Unit Graph
ug = SimpleDiGraphFromIterator((Edge(i) for i in projs))
nn = nv(ug)
nodecolor = [in(i,eunits) ? RGB(1,0.2,0.2) : in(i,iunits) ? RGB(0.2,0.2,1) : RGB(0.4,0.4,0.4) for i in 1:nn]
p=gplot(ug,unitposition[unitgood,1][1:nn],unitposition[unitgood,2][1:nn],nodelabel=1:nn,edgestrokec="gray30",nodefillc=map(c->coloralpha(c,0.5),nodecolor),
nodelabelc=nodecolor,nodesize=3,arrowangleoffset=10/180*pi,arrowlengthfrac=0.025)

using Compose,Cairo,Fontconfig
draw(PDF(joinpath(resultdir,"circuitgraph.pdf"), 20Compose.cm, 20Compose.
