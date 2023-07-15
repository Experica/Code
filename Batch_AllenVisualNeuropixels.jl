using CSV,DataFrames,NeuroAnalysis,Mmap,DSP,StatsBase,StatsPlots,JLD2,Images,ProgressMeter,PyCall

rootdir = raw"E:\SpikeGLXData\allen-brain-observatory\visual-coding-neuropixels"
cachedir = joinpath(rootdir,"ecephys-cache")
rawdatadir = joinpath(rootdir,"raw-data")
resultroot = "Z:/"

sessions = CSV.read(joinpath(cachedir,"sessions.csv"),DataFrame)
probes = CSV.read(joinpath(cachedir,"probes.csv"),DataFrame)
channels = CSV.read(joinpath(cachedir,"channels.csv"),DataFrame)
units = CSV.read(joinpath(cachedir,"units.csv"),DataFrame)

sp=innerjoin(rename(filter(i->i.phase=="PXI",probes),:id=>:probe_id,:ecephys_session_id=>:session_id,:lfp_sampling_rate=>:lffs,:sampling_rate=>:apfs),
    rename(sessions,:id=>:session_id),
    on=:session_id)
spc = innerjoin(transform(rename(channels,:ecephys_probe_id=>:probe_id,:id=>:channel_id),
                :ecephys_structure_acronym=>ByRow(i->ismissing(i) ? missing : Symbol(i))=>:neurostructure),
                sp,on=:probe_id)
sps = combine(groupby(dropmissing(spc,:neurostructure),[:probe_id,:session_id]),[:lffs,:apfs].=>first.=>identity,
            [:probe_vertical_position,:neurostructure]=>((p,s)->[Dict(k=>extrema(p[s.==k]).-20 for k in unique(s))])=>last)
spu = innerjoin(select(units,[:duration,:spread],:id=>:unit_id,:ecephys_channel_id=>:channel_id,:PT_ratio=>:ptr,:quality=>ByRow(i->i=="good")=>:unit_good),
                select(spc,[:channel_id,:probe_id,:session_id],:probe_vertical_position=>(i->i.-20)=>:chposy),on=:channel_id)


spudf = combine(groupby(filter(i->i.unit_good,spu),[:probe_id,:session_id]),[:chposy,:spread,:duration,:ptr,:unit_good].=>(i->[i]).=>identity)
select!(spudf,[:session_id,:probe_id],[:chposy,:unit_good]=>ByRow((y,g)->[unitdensity(y,w=replace(g,0=>1.2),bw=40,step=20,r=100,s=1.4)...].*[1e9,1])=>[:Density,:depth],
            [:chposy,:spread]=>ByRow((y,v)->unitdensity(y,w=v,wfun=mean,bw=40,step=20,s=1.4).n)=>:Spread,
            [:chposy,:duration]=>ByRow((y,v)->unitdensity(y,w=v,wfun=mean,bw=40,step=20,s=1.4).n)=>:Duration,
            [:chposy,:ptr]=>ByRow((y,v)->unitdensity(y,w=v,wfun=mean,bw=40,step=20,s=1.4).n)=>:PTR)
leftjoin!(spudf,sps[:,[:probe_id,:session_id,:neurostructure]],on=[:probe_id,:session_id])


ii = 1
session_id = sps.session_id[ii]
probe_id = sps.probe_id[ii]
lffs = sps.lffs[ii]
apfs = sps.apfs[ii]
neurostructure = sps.neurostructure[ii]
process_flash_allen(session_id,probe_id,lffs,apfs,neurostructure)



function process_flash_allen(session_id,probe_id,lffs,apfs,neurostructure;np=pyimport("numpy"),figfmt = [".png"],layer=nothing)
    sessionresultdir = joinpath(resultroot,"$session_id")
    sessionrawdir = joinpath(rawdatadir,"$session_id")
    ctc = CSV.read(joinpath(cachedir,"session_$(session_id)_flashes.csv"),DataFrame)
    condon = (ctc.start_time.-ctc.t0)*1000 # ms
    cond = condin(ctc[:,[:color]])

    resultdir = joinpath(sessionresultdir,"$probe_id")
    proberawdir = joinpath(sessionrawdir,"$probe_id")
    apbinfile = joinpath(proberawdir,"flash_g0_tcat.imec0.ap.kilosort3.whiten.dat")
    apkilosortdir = joinpath(proberawdir,"flash_g0_tcat.imec0.ap.kilosort3.phy")
    winvraw = np.load(joinpath(apkilosortdir,"whitening_mat_inv_raw.npy"))
    chmapraw = dropdims(np.load(joinpath(apkilosortdir,"channel_map_raw.npy")),dims=2).+1
    lfbinfile = joinpath(proberawdir,"lfp_band.dat")
    isdir(resultdir) || mkpath(resultdir)
    t0 = 1200*1000 # t0 ms for catgt flash ap segment

    # Prepare AP
    nch = 383 # kilosort whiten data
    nsample = Int(stat(apbinfile).size/2/nch)
    hy = 20 # for Neuropixels 1.0 (μm)
    pversion = 1
    pnrow = 192
    pncol = 2
    exch = [192] # ref ch 192 for Neuropixels 1.0
    probe_id in [848037568,841431754,841435557] && push!(exch,266) # bad channel for these probes
    exchmask = chmasknp(exch,pnrow,pncol)
    depths = hy*(0:(pnrow-1))
    mmapbinfile = mmap(apbinfile,Matrix{Int16},(nch,nsample),0)

    # epoch AP
    epoch = [-40 150]
    epochs = condon.-t0.+epoch
    # All AP epochs, unwhiten, bandpass filtered,
    # with all channels mapped in the shape of probe where excluded channels are replaced with rand local average
    ys = fill2mask(epochsamplenp(mmapbinfile,apfs,epochs,1:nch;bandpass=[500,3000],whiten=winvraw),exchmask,chmap=chmapraw,randreplace=true)
    # 3 fold downsampling(10kHz)
    ys = resample(ys,1//3,dims=3)
    fs = apfs/3
    baseindex = epoch2sampleindex([0 50],fs) # [-40, 10] ms

    # power spectrum of same depth
    @views pys = ys[.!exchmask,:,:] # exclude and flat channels
    chpos = vcat(chpositionnp(pversion)[.!exchmask]...) # exclude and flat channel positions 
    chgi = [findall(chpos[:,2].==up) for up in unique(chpos[:,2])] # group channels with same depth
    @views ppss,psfreqs = powerspectrum(pys[:,baseindex[end]+1:end,:],fs;freqrange=[500,3000])
    pss = Array{Float64}(undef,length(chgi),size(ppss)[2:end]...)
    @views foreach(i->pss[i,:,:] = mean(ppss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    tps = dropdims(mean(pss,dims=2),dims=2) # freq average
    @views cfps = Dict(condstring(r)=>
        dropdims(mean(pss[:,:,r.i],dims=3),dims=3) # trial average       
    for r in eachrow(cond))

    # power contrast of same depth
    prmss = mapwindow(x->rms(x)^2,pys,(1,101,1),border="symmetric") # 10ms rms window
    rmss = Array{Float64}(undef,length(chgi),size(prmss)[2:end]...)
    @views foreach(i->rmss[i,:,:] = mean(prmss[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    @views crms = Dict(condstring(r)=>
        stfilter(dropdims(mean(rmss[:,:,r.i],dims=3),dims=3),temporaltype=:rcb,ti=baseindex) # trial average
    for r in eachrow(cond))
    times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(pys,2))

    # local coherence
    @views lcs,lcfreqs = localcoherence(pys[:,baseindex[end]+1:end,:],chpos,fs;freqrange=[500,3000],lr=55,sigma=25,chgroupdim=2)
    tlc = dropdims(mean(lcs,dims=2),dims=2) # freq average
    @views cflc = Dict(condstring(r)=>
        dropdims(mean(lcs[:,:,r.i],dims=3),dims=3) # trial average       
    for r in eachrow(cond))


    plotanalog(tps;hy,color=:heat,n=mean(tps,dims=2),xlabel="Trial",xunit="",cunit=:fr,layer)
    foreach(ext->savefig(joinpath(resultdir,"tPower$ext")),figfmt)
    plotanalog(tlc;hy,xlabel="Trial",xunit="",color=:heat,cunit=:fr,n=mean(tlc,dims=2),layer)
    foreach(ext->savefig(joinpath(resultdir,"tCoherence$ext")),figfmt)
    fpslims = extrema(mapreduce(extrema,union,values(cfps)))
    rmslim = mapreduce(x->maximum(abs.(x)),max,values(crms))
    flclims = extrema(mapreduce(extrema,union,values(cflc)))
    for k in keys(crms)
        plotanalog(cfps[k];hy,x=psfreqs,xlabel="Frequency",xunit="Hz",clims=fpslims,color=:heat,cunit=:fr,n=mean(cfps[k],dims=2),layer)
        foreach(ext->savefig(joinpath(resultdir,"$(k)_fPower$ext")),figfmt)
        plotanalog(crms[k];x=times,hy,clims=(-rmslim,rmslim),layer)
        foreach(ext->savefig(joinpath(resultdir,"$(k)_dRMS$ext")),figfmt)
        plotanalog(cflc[k];hy,x=lcfreqs,xlabel="Frequency",xunit="Hz",clims=flclims,color=:heat,cunit=:fr,n=mean(cflc[k],dims=2),layer)
        foreach(ext->savefig(joinpath(resultdir,"$(k)_fCoherence$ext")),figfmt)
    end
    jldsave(joinpath(resultdir,"ap.jld2");tps,cfps,crms,tlc,cflc,times,psfreqs,lcfreqs,depths,baseindex,fs,session_id,probe_id,neurostructure)


    # Prepare LF
    nch = 384
    nsample = Int(stat(lfbinfile).size/2/nch)
    hy = 20 # for Neuropixels 1.0 (μm)
    pversion = 1
    pnrow = 192
    pncol = 2
    exch = [192] # ref ch 192 for Neuropixels 1.0
    probe_id in [848037568,841431754,841435557] && push!(exch,266) # bad channel for these probes
    exchmask = chmasknp(exch,pnrow,pncol) 
    depths = hy*(0:(pnrow-1))
    mmlfbinfile = mmap(lfbinfile,Matrix{Int16},(nch,nsample),0)

    # epoch LFP
    epoch = [-40 150]
    epochs = condon.+epoch
    # All LFP epochs, line noise(60,120,180Hz) removed, bandpass filtered,
    # with all channels mapped in the shape of probe where excluded channels are replaced with local average
    ys = fill2mask(epochsamplenp(mmlfbinfile,lffs,epochs,1:nch;bandpass=[1,100]),exchmask)
    # 2.5 fold downsampling(1kHz)
    ys = resample(ys,1/2.5,dims=3)
    fs = lffs/2.5
    baseindex = epoch2sampleindex([0 50],fs) # [-40, 10] ms

    # LFP of same depth
    @views pys = ys[.!exchmask,:,:] # exclude and flat channels
    chpos = vcat(chpositionnp(pversion)[.!exchmask]...) # exclude and flat channel positions 
    chgi = [findall(chpos[:,2].==up) for up in unique(chpos[:,2])] # group channels with same depth
    lfps = Array{Float64}(undef,length(chgi),size(pys)[2:end]...)
    @views foreach(i->lfps[i,:,:] = mean(pys[chgi[i],:,:],dims=1),eachindex(chgi)) # average of same depth
    @views clfp = Dict(condstring(r)=>
        dropdims(mean(lfps[:,:,r.i],dims=3),dims=3) # trial average
    for r in eachrow(cond))
    
    # CSD of same depth
    ccsd = Dict(k=>stfilter(csd(v,h=hy),temporaltype=:sub,ti=baseindex) for (k,v) in clfp)
    times = range(start=epoch[1],step=1/fs/SecondPerUnit,length=size(pys,2))

    lfplim = mapreduce(x->maximum(abs.(x)),max,values(clfp))
    csdlim = mapreduce(x->maximum(abs.(x)),max,values(ccsd))
    for k in keys(clfp)
        plotanalog(clfp[k];x=times,hy,clims=(-lfplim,lfplim))
        foreach(ext->savefig(joinpath(resultdir,"$(k)_LFP$ext")),figfmt)
        plotanalog(ccsd[k];x=times,color=:RdBu,hy,clims=(-csdlim,csdlim))
        foreach(ext->savefig(joinpath(resultdir,"$(k)_dCSD$ext")),figfmt)
    end
    jldsave(joinpath(resultdir,"lf.jld2");clfp,ccsd,times,depths,fs,baseindex,session_id,probe_id,neurostructure)
end


## batch each probe for each session
np = pyimport("numpy")
@showprogress "Batch Allen Session Probe ... " for r in eachrow(sps)
    process_flash_allen(r.session_id,r.probe_id,r.lffs,r.apfs,r.neurostructure;np)
end