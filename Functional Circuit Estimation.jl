using NeuroAnalysis,Statistics,Plots,Mmap,Images,StatsBase,LightGraphs,Interact,Dierckx,PyCall


# Prepare Dataset
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"
settimeunit(1)

subject = "AE9";recordsession = "";recordsite = "u003"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
sitedir = joinpath(resultroot,subject,siteid)
testid = "$(siteid)_000"
resultdir = joinpath(sitedir,testid)
isdir(resultdir) || mkpath(resultdir)

dataset = prepare(joinpath(dataexportroot,"$testid.mat"))
nsavech=dataset["lf"]["meta"]["nSavedChans"]
nsample=dataset["lf"]["meta"]["nFileSamp"]
fs=dataset["lf"]["meta"]["fs"]
nch = dataset["lf"]["meta"]["snsApLfSy"][2]
refchmask = refchmaskim(dataset["ap"]["meta"])

lfregex = Regex("^$(uppercase(testid))[A-Za-z0-9_]*.imec.lf.bin")
lffile = joinpath(dataroot,filter(i->occursin(lfregex,i),readdir(dataroot))[1])
mmlfp = Mmap.mmap(lffile,Matrix{Int16},(nsavech,nsample),0)


# Condition Test
ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
spike = dataset["spike"]
eval.([:($(Symbol(k))=spike[$k]) for k in keys(spike)])
condon = ex["CondTest"]["CondOn"]
condoff = ex["CondTest"]["CondOff"]
minconddur=minimum(condoff-condon)
histogram(condoff-condon,nbins=20,title="Condition Duration(Set to $conddur)")


# Laminar CSD
epochdur = 0.15
epoch = [0 epochdur]
epochs = condon.+epoch
ys=reshape2ref(subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),refchmask) # all epoch LFP segments for each channel of probe in shape

col=1
mcys=dropdims(mean(ys[:,col,:,:],dims=3),dims=3)
plotanalog(mcys,fs=fs)
foreach(i->savefig(joinpath(resultdir,"column_$(col)_lfp$i")),[".png",".svg"])

mccsd = dropdims(mean(csd(ys[:,col,:,:]),dims=3),dims=3)
plotanalog(imfilter(mccsd,Kernel.gaussian((1,1))),fs=fs)
foreach(i->savefig(joinpath(resultdir,"column_$(col)_csd$i")),[".png",".svg"])

# Column Mean Normalized CSD
pys = cat((ys[:,i,:,:] for i in 1:2)...,dims=3)
pcsd = csd(pys)
baseline = pcsd[:,epoch2samplerange([0 0.015],fs),:] # 0-15ms csd as baseline
ncsd=pcsd.-mean(baseline,dims=2) # Δcsd relative to baseline
mncsd = dropdims(mean(ncsd,dims=3),dims=3)
mncsd[1,:].=0;mncsd[end,:].=0
mncsd = imfilter(mncsd,Kernel.gaussian((2,1)))
h=20;depths = h*(1:size(mncsd,1))
plotanalog(mncsd,fs=fs,y=depths,color=:RdBu)
foreach(i->savefig(joinpath(resultdir,"column_mean_normalized_csd$i")),[".png",".svg"])
save(joinpath(resultdir,"csd.jld2"),"csd",mncsd,"depth",depths,"fs",fs)


# Depth Power Spectrum
epochdur = 0.2
epoch = [0 epochdur]
epochs = condon.+epoch
ys=reshape2ref(subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),refchmask)
pys = cat((ys[:,i,:,:] for i in 1:2)...,dims=3)
ps,freq = powerspectrum(pys,fs)
mps = dropdims(mean(ps,dims=3),dims=3)
plotanalog(imfilter(mps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:db,color=:fire)

epoch = [-epochdur 0]
epochs = condon.+epoch
ys=reshape2ref(subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),refchmask)
pys = cat((ys[:,i,:,:] for i in 1:2)...,dims=3)
bps,freq = powerspectrum(pys,fs)
mbps =dropdims(mean(bps,dims=3),dims=3)
plotanalog(imfilter(mbps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:db,color=:fire)

rcps = ps./bps.-1
mrcps = imfilter(dropdims(mean(rcps,dims=3),dims=3),Kernel.gaussian((1,1)))
plotanalog(mrcps,x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:db,color=:fire)
foreach(i->savefig(joinpath(resultdir,"depth_rcpowerspectrum$i")),[".png",".svg"])
save(joinpath(resultdir,"powerspectrum.jld2"),"rcps",rcps,"depth",depths,"freq",freq)


# Unit Position
plotunitposition(spike)
foreach(i->savefig(joinpath(resultdir,"Unit_Position$i")),[".png",".svg"])
save(joinpath(resultdir,"spike.jld2"),"spike",spike)

# Unit Spike Trian
epochext = 0.1
@manipulate for u in 1:length(unitspike)
ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
plotspiketrain(ys,timeline=[0,minconddur],title="Unit_$(unitid[u])")
end

# Unit Depth PSTH
epochdur = 0.15
epoch = [0 epochdur]
epochs = condon.+epoch
bw = 0.002
psthbins = epoch[1]:bw:epoch[2]
h=20

baseindex = epoch2samplerange([0 0.015],1/bw) # 0-15ms of psth as baseline
unitpsth = map(ust->psth(subrv(ust,epochs,isminzero=true)[1],psthbins,israte=true,normfun=x->x.-mean(x[baseindex])),unitspike)
depthpsth,x,depths = spacepsth(unitpsth,unitposition,spacebinedges=h*(0:size(refchmask,1)+1))
depthpsth = imfilter(depthpsth,Kernel.gaussian((1,1)))
plotpsth(depthpsth,x,depths)
foreach(i->savefig(joinpath(resultdir,"depthpsth$i")),[".png",".svg"])
save(joinpath(resultdir,"depthpsth.jld2"),"depthpsth",depthpsth,"depth",depths,"x",x)








@manipulate for er1=0.03,erl=0.01,fpext=0.015, dr1=0:3000,drl=100
uid = dr1 .<= unitposition[:,2] .< dr1+drl
epochs = condon[1:2:end].+[er1 er1+erl]
ys = mapreduce(ust->flatrvv(subrv(ust,epochs)[1])[1],vcat,unitspike[uid])

sepochs = ys.+[-fpext fpext]
stfp = reshape2ref(subrm(mmlfp,fs,sepochs,chs=1:nch,meta=dataset["lf"]["meta"]),refchmask)


pys = cat((stfp[:,i,:,:] for i in 1:2)...,dims=3)
pcsd = csd(pys)
# baseline = pcsd[:,epoch2samplerange([fpext-0.003 fpext+0.003],fs),:]
baseline = pcsd[:,epoch2samplerange([0 0.015],fs),:]
ncsd=pcsd.-mean(baseline,dims=2)
#ncsd = pcsd
mncsd = dropdims(mean(ncsd,dims=3),dims=3)
mncsd[1,:].=0;mncsd[end,:].=0
mncsd = imfilter(mncsd,Kernel.gaussian((2,1)))
plotanalog(mncsd,fs=fs,color=:RdBu)
hline!([dr1,dr1+drl],linestyle=:dash)
vline!([fpext*1000],leg=false)
end


















# Single Unit Binary Spike Trian of Conditions
bepochext = -0.3
bepoch = [-bepochext minconddur]
bepochs = condon.+bepoch
spikebins = bepoch[1]:0.001:bepoch[2]
subst = map(ust->float.(histmatrix(subrv(ust,bepochs,isminzero=true,shift=bepochext)[1],spikebins)[1])',unitspike[unitgood])
# Single Unit Correlogram and Circuit
lag=50;suid = unitid[unitgood]
ccgs,x,ccgis,projs,eunits,iunits = circuitestimate(subst,lag=lag)

@manipulate for i in 1:length(ccgs)
    title = "Correlogram between Unit $(suid[ccgis[i][1]]) and $(suid[ccgis[i][2]])"
    bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
end

for i in 1:length(ccgs)
    title = "Correlogram between Unit $(suid[ccgis[i][1]]) and $(suid[ccgis[i][2]])"
    bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
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
