## automatic phy curation
using CSV,DataFrames

phydir = raw"X:\AG2\AG2_V1_ODR1\F0F1F2F3C0C1C0C1C2C3C4C5H0H1H2H3I0I1I2I3O0O1O2O3O4.imec0.ap.kilosort3.phy"
params = CSV.read(joinpath(phydir,"params.py"),DataFrame,header=0,delim='=',quotechar=''',stripwhitespace=true)
params = Dict(r.Column1=>r.Column2 for r in eachrow(params))
nsample = parse(Int,params["nsample"])
fs = parse(Float64,params["sample_rate"])
duration = nsample/fs

clustergroup = CSV.read(joinpath(phydir,"cluster_group.tsv"),DataFrame)
clusterinfo = CSV.read(joinpath(phydir,"cluster_info.tsv"),DataFrame)

# clusters not assigned by phy operator
ti = ismissing.(clusterinfo.group)

## low firing rate cluster as noise
minfr = 0.15
fi = clusterinfo.fr .<= minfr
vi = ti .& fi

clusterinfo.group[vi] .= "noise"
append!(clustergroup,DataFrame(cluster_id=clusterinfo.cluster_id[vi],group="noise")) |> sort!


CSV.write(joinpath(phydir,"cluster_info.tsv"),clusterinfo,delim='\t')
CSV.write(joinpath(phydir,"cluster_group.tsv"),clustergroup,delim='\t')







