using NeuroAnalysis,Statistics,Plots,Mmap,Images,StatsBase,LightGraphs,Interact,Dierckx,PyCall


# Prepare Dataset
dataroot = "../Data"
dataexportroot = "../DataExport"
resultroot = "../Result"

subject = "AF5";recordsession = "HLV1";recordsite = "ODL3";test="Flash2Color_5"
siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
datadir = joinpath(dataroot,subject,siteid)
resultsitedir = joinpath(resultroot,subject,siteid)
testid = join(filter(!isempty,[siteid,test]),"_")
resultdir = joinpath(resultsitedir,testid)
isdir(resultdir) || mkpath(resultdir)

dataset = prepare(joinpath(dataexportroot,subject,"$testid.mat"))

# Condition Test
ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
spike = dataset["spike"]
eval.([:($(Symbol(k))=spike[$k]) for k in keys(spike)])
condon = ex["CondTest"]["CondOn"]
condoff = ex["CondTest"]["CondOff"]
minconddur=minimum(condoff-condon)
histogram(condoff-condon,nbins=20,title="Condition Duration(Set to $conddur)")

# Prepare LFP
lffile = matchfile(Regex("^$(testid)[A-Za-z0-9_]*.imec.lf.bin"),dir = datadir,adddir=true)[1]
nsavech=dataset["lf"]["meta"]["nSavedChans"]
nsample=dataset["lf"]["meta"]["nFileSamp"]
fs=dataset["lf"]["meta"]["fs"]
nch = dataset["lf"]["meta"]["snsApLfSy"][2]
hx,hy,hz = dataset["lf"]["meta"]["probespacing"]
badchmask = badchmaskim(dataset)
pnrow,pncol = size(badchmask)
mmlf = Mmap.mmap(lffile,Matrix{Int16},(nsavech,nsample),0)

# Laminar LFP and CSD
epochdur = 150
epoch = [0 epochdur]
epochs = condon.+epoch
ys=reshape2mask(subrm(mmlf,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),badchmask) # all LFP epochs for each channel of probe in shape

col=1
mcys=dropdims(mean(ys[:,col,:,:],dims=3),dims=3)
plotanalog(mcys,fs=fs,cunit=:uv)
foreach(i->savefig(joinpath(resultdir,"Column_$(col)_Mean_LFP$i")),[".png",".svg"])

mccsd = dropdims(mean(csd(ys[:,col,:,:],h=hy),dims=3),dims=3)
plotanalog(imfilter(mccsd,Kernel.gaussian((1,1))),fs=fs)
foreach(i->savefig(joinpath(resultdir,"Column_$(col)_Mean_CSD$i")),[".png",".svg"])

# Column Mean ΔCSD
pys = cat((ys[:,i,:,:] for i in 1:pncol)...,dims=3)
pcsd = csd(pys,h=hy)
basedur = 15
baseindex = epoch2samplerange([0 basedur],fs)
dcsd=stfilter(pcsd,temporaltype=:sub,ti=baseindex,spatialtype=:annulus,ir=3,or=10)
mdcsd = dropdims(mean(dcsd,dims=3),dims=3)
mdcsd[[1,end],:].=0
depths = hy*(1:size(mdcsd,1))
plotanalog(imfilter(mdcsd,Kernel.gaussian((2,1))),fs=fs,y=depths,color=:RdBu)
foreach(i->savefig(joinpath(resultdir,"Columns_Mean_dCSD$i")),[".png",".svg"])
save(joinpath(resultdir,"csd.jld2"),"csd",mdcsd,"depth",depths,"fs",fs,"log",ex["Log"],"color","$(ex["Param"]["ColorSpace"])_$(ex["Param"]["Color"])")

# Depth Power Spectrum
epochdur = 200
epoch = [0 epochdur]
epochs = condon.+epoch
ys=reshape2mask(subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),badchmask)
pys = cat((ys[:,i,:,:] for i in 1:2)...,dims=3)
ps,freq = powerspectrum(pys,fs)
mps = dropdims(mean(ps,dims=3),dims=3)
plotanalog(imfilter(mps,Kernel.gaussian((1,1))),x=freq,y=depths,xlabel="Freqency",xunit="Hz",timeline=[],cunit=:db,color=:fire)

epoch = [-epochdur 0]
epochs = condon.+epoch
ys=reshape2mask(subrm(mmlfp,fs,epochs,chs=1:nch,meta=dataset["lf"]["meta"]),badchmask)
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
epochext = 100
@manipulate for u in 1:length(unitspike)
ys,ns,ws,is = subrv(unitspike[u],condon.-epochext,condoff.+epochext,isminzero=true,shift=epochext)
plotspiketrain(ys,timeline=[0,minconddur],title="Unit_$(unitid[u])")
end

# Unit Depth PSTH
epochdur = 150
epoch = [0 epochdur]
epochs = condon.+epoch
bw = 2
psthbins = epoch[1]:bw:epoch[2]

baseindex = epoch2samplerange([0 basedur],1/(bw*SecondPerUnit))
normfun = x->x.-mean(x[baseindex])
unitpsth = map(ust->psth(subrv(ust,epochs,isminzero=true)[1],psthbins,israte=true,normfun=normfun),unitspike[unitgood])
depthpsth,x,depths,depthnunit = spacepsth(unitpsth,unitposition[unitgood,:],spacebinedges=hy*(0:pnrow+1))

plotpsth(imfilter(depthpsth,Kernel.gaussian((1,1))),x,depths,color=:minmax,n=depthnunit)
foreach(i->savefig(joinpath(resultdir,"DepthPSTH$i")),[".png",".svg"])
save(joinpath(resultdir,"depthpsth.jld2"),"depthpsth",depthpsth,"depth",depths,"x",x,"n",depthnunit)







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
bepochext = timetounit(-300)
bepoch = [-bepochext minconddur]
bepochs = ccondon.+bepoch
spikebins = bepoch[1]:timetounit(1):bepoch[2]
subst = map(ust->float.(histmatrix(subrv(ust,bepochs,isminzero=true,shift=bepochext)[1],spikebins)[1])',unitspike[unitgood])
# Single Unit Correlogram and Circuit
lag=50
ccgs,x,ccgis,projs,eunits,iunits,projweights = circuitestimate(subst,lag=lag,unitid=unitid[unitgood],condis=ccond.i,minepoch=5,minspike=5)



@manipulate for i in 1:length(ccgs)
    title = "Correlogram between Single Unit $(ccgis[i][1]) and $(ccgis[i][2])"
    bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
end

for i in 1:length(ccgs)
    title = "Correlogram between Single Unit $(ccgis[i][1]) and $(ccgis[i][2])"
    bar(x,ccgs[i],bar_width=1,legend=false,color=:gray15,linecolor=:match,title=title,xlabel="Time (ms)",ylabel="Coincidence/Spike",grid=(:x,0.4),xtick=[-lag,0,lag])
    foreach(i->savefig(joinpath(resultdir,"$title$i")),[".png",".svg"])
end
save(joinpath(resultdir,"circuit.jld2"),"projs",projs,"eunits",eunits,"iunits",iunits)
