includet("Batch/Batch.jl")


## Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "I:/",
    :dataexportroot => "Y:/",
    :resultroot => "Y:/",
    :stimuliroot => "S:/",
    :spikesorter => "kilosort3")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))


## Query Tests
tests = select!(filter(meta) do r
                    r.Subject_ID in ["Test", "AG2"] &&
                    # r.RecordSession == "V1" &&
                    # r.RecordSite == "ODR1" &&
                    # r.ID in ["ISIEpochOri8"] &&
                    r.sourceformat == "Imager"
                end,
                [:files,:ID,:UUID,:sourceformat])


## Setup Param
files = ".\\AG1\\AG1_V1V2_Full_ISICycleOri_0.mat"
tests.files[2]
tests = tests[1:1,:]


## Batch Tests
batchtests(tests,param,plot=true)