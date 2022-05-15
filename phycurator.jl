## phy curation
using CSV,DataFrames,PyCall,Combinatorics,LinearAlgebra,ProgressMeter
np = pyimport("numpy")

function findoverlapspike(st,win)
    oi=Int[];n=length(st);i=1
    @label start
    for j in i+1:n
        if st[j]-st[i] > win
            i=j
            @goto start
        end
        push!(oi,j)
    end
    oi
end
function findoverlapspike(st,sa,win)
    oi=Int[];n=length(st);i=1
    @label start
    for j in i+1:n
        if st[j]-st[i] > win
            i=j
            @goto start
        end
        if sa[j] > sa[i]
            push!(oi,i)
            i=j
            @goto start
        end
        push!(oi,j)
    end
    oi
end


## load phy
phydir = raw"X:\AG1\AG1_V1_ODL17\F0F1F2F3C0C1H0H1H2H3O0O1O2O3.imec0.ap.kilosort3.phy"
backupdir = joinpath(phydir,"backup")
params = CSV.read(joinpath(phydir,"params.py"),DataFrame,header=0,delim='=',quotechar=''',stripwhitespace=true)
params = Dict(r.Column1=>r.Column2 for r in eachrow(params))
nsample = parse(Int,params["nsample"])
fs = parse(Float64,params["sample_rate"])
duration = nsample/fs


## remove putative double counted spikes
loaddir = ispath(backupdir) ? backupdir : phydir
spiketime = np.load(joinpath(loaddir,"spike_times.npy")) |> vec
spikecluster = np.load(joinpath(loaddir,"spike_clusters.npy")) |> vec
spiketemplate = np.load(joinpath(loaddir,"spike_templates.npy")) |> vec
spikeamplitude = np.load(joinpath(loaddir,"amplitudes.npy")) |> vec
templates = np.load(joinpath(phydir,"templates.npy"))
chpos = np.load(joinpath(phydir,"channel_positions.npy"))

overlapwinin = round(Int,0.4 * 1e-3 * fs)
clusterid = unique(spikecluster)
rmspikeidx = Int[]
minclustersize = 300
# remove overlap spikes within cluster
@showprogress for c in clusterid
    ci = findall(spikecluster.==c)
    length(ci) < minclustersize && continue
    ci = ci[sortperm(@view spiketime[ci])]
    @views crmidx = findoverlapspike(spiketime[ci],overlapwinin)
    append!(rmspikeidx,ci[crmidx])
end

sort!(unique!(rmspikeidx))
deleteat!(spiketime,rmspikeidx)
deleteat!(spikecluster,rmspikeidx)
deleteat!(spiketemplate,rmspikeidx)
deleteat!(spikeamplitude,rmspikeidx)


# remove overlap spikes between nearby clusters
templatechamp = map(i->i[2]-i[1],dropdims(extrema(templates,dims=2),dims=2))
templatechidx = map(i->i[2],dropdims(argmax(templatechamp,dims=2),dims=2))

overlapwinbetween = round(Int,0.35 * 1e-3 * fs)
clusterid = unique(spikecluster)
clusteridx = [findall(spikecluster.==c) for c in clusterid]
clustersize = length.(clusteridx)
clusterpos = chpos[templatechidx[clusterid.+1],:]
rmspikeidx = Int[]
radius = 60 # closest channel distance for Neuropixels1 = 25.6, dx=32, dy=20

@showprogress for (c1,c2) in combinations(1:length(clusterid),2)
    (clustersize[c1] < minclustersize || clustersize[c2] < minclustersize || norm(clusterpos[c2,:] .- clusterpos[c1,:]) > radius) && continue
    ci = [clusteridx[c1];clusteridx[c2]]
    ci = ci[sortperm(@view spiketime[ci])]
    @views crmidx = findoverlapspike(spiketime[ci],spikeamplitude[ci],overlapwinbetween)
    append!(rmspikeidx,ci[crmidx])
end

sort!(unique!(rmspikeidx))
deleteat!(spiketime,rmspikeidx)
deleteat!(spikecluster,rmspikeidx)
deleteat!(spiketemplate,rmspikeidx)
deleteat!(spikeamplitude,rmspikeidx)

# overwrite phy
if !ispath(backupdir)
    mkpath(backupdir)
    cp(joinpath(phydir,"spike_times.npy"),joinpath(backupdir,"spike_times.npy"))
    cp(joinpath(phydir,"spike_clusters.npy"),joinpath(backupdir,"spike_clusters.npy"))
    cp(joinpath(phydir,"spike_templates.npy"),joinpath(backupdir,"spike_templates.npy"))
    cp(joinpath(phydir,"amplitudes.npy"),joinpath(backupdir,"amplitudes.npy"))
end
np.save(joinpath(phydir,"spike_times.npy"),spiketime)
np.save(joinpath(phydir,"spike_clusters.npy"),spikecluster)
np.save(joinpath(phydir,"spike_templates.npy"),spiketemplate)
np.save(joinpath(phydir,"amplitudes.npy"),spikeamplitude)


### Now Do Phy Manual Curation ###


## for clusters not handled, assign low firing rate cluster as noise
clustergroup = CSV.read(joinpath(phydir,"cluster_group.tsv"),DataFrame)
clusterinfo = CSV.read(joinpath(phydir,"cluster_info.tsv"),DataFrame)
ti = ismissing.(clusterinfo.group)
minfr = 0.5
fi = clusterinfo.fr .<= minfr
vi = ti .& fi

clusterinfo.group[vi] .= "noise"
append!(clustergroup,DataFrame(cluster_id=clusterinfo.cluster_id[vi],group="noise")) |> sort!

CSV.write(joinpath(phydir,"cluster_info.tsv"),clusterinfo,delim='\t')
CSV.write(joinpath(phydir,"cluster_group.tsv"),clustergroup,delim='\t')



## finally clean up
rm(backupdir,force=true,recursive=true)

