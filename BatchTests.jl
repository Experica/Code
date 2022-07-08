includet("Batch.jl")


## Prepare Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "X:/",
    :dataexportroot => "Y:/",
    :resultroot => "Z:/",
    :stimuliroot => "S:/",
    :spikesorter => "kilosort3")
# param = Dict{Any,Any}(
#     :dataroot => "I:\\AG1\\AG1_V1V2_Full",
#     :dataexportroot => "I:\\",
#     :resultroot => "I:\\ISI_analysis")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))


## Query Tests
tests = select!(filter(meta) do r
                    r.Subject_ID in ["AG1", "AG2"] &&
                    # r.RecordSession == "V1" &&
                    # r.RecordSite == "ODR1" &&
                    r.ID == "HartleySubspace" &&
                    r.sourceformat == "SpikeGLX"
                end,
                [:files,:ID,:UUID,:sourceformat])
## Setup Param

files = ".\\AG1\\AG1_V1_ODL1_HartleySubspace_0.mat"
pyplot()

param[:model] = [:STA]

param[:model] = [:ePPR]
param[:epprparam] = (ndelay=1, nft=[3,3,3], lambda=320000)
param[:batchunit] = batchunit


## Batch Tests
batchtests(tests,param,plot=true)
