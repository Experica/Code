includet("Batch/Batch.jl")


## Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "X:/",
    :dataexportroot => "Y:/",
    :resultroot => "Z:/",
    :stimuliroot => "S:/",
    :spikesorter => "kilosort3")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))


## Query Tests
tests = select!(filter(meta) do r
                    r.Subject_ID in ["AG5"] &&
                    # r.RecordSession == "V1" &&
                    # r.RecordSite == "1R" &&
                    r.ID in ["Color"] &&
                    r.sourceformat == "SpikeGLX"
                end,
                [:files,:ID,:UUID,:sourceformat])


## Setup Param
files = ".\\AG1\\AG1_V1_ODL1_Flash2Color_0.mat"
files = ".\\AG2\\AG2_V1_ODL3_OriSF_0.mat"
pythonplot()
pyplot()

# only process and save spike data
param[:onlyspike] = true
# images responses linear model
param[:model] = [:STA]
# images responses ePPR model
param[:model] = [:ePPR]
param[:epprparam] = (ndelay=1, nft=[3,3,3], lambda=320000)
param[:batchunit] = batchunit


## Batch Tests
batchtests(tests,param,plot=true)
