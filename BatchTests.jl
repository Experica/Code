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
                    r.Subject_ID in ["AG1", "AG2"] &&
                    # r.RecordSession == "V1" &&
                    # r.RecordSite == "ODR1" &&
                    r.ID in ["OriSF"] &&
                    r.sourceformat == "SpikeGLX"
                end,
                [:files,:ID,:UUID,:sourceformat])


## Setup Param
files = ".\\AG1\\AG1_V1V2_Full_ISICycleOri_0.mat"
tests[1:4,:].files
tests=tests[1:4,:]

files = ".\\AG1\\AG1_V1_ODL1_Flash2Color_0.mat"
files = ".\\AG2\\AG2_V1_ODL3_OriSF_0.mat"
pythonplot()

param[:model] = [:STA]

param[:model] = [:ePPR]
param[:epprparam] = (ndelay=1, nft=[3,3,3], lambda=320000)
param[:batchunit] = batchunit


## Batch Tests
batchtests(tests,param,plot=true)
