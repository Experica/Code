includet("Batch.jl")


## Prepare Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "X:/",
    :dataexportroot => "Y:/",
    :resultroot => "Z:/",
    :stimuliroot => "D:/NaturalStimuli")
# param = Dict{Any,Any}(
#     :dataroot => "I:\\AG1\\AG1_V1V2_Full",
#     :dataexportroot => "I:\\",
#     :resultroot => "I:\\ISI_analysis")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))


## Query Tests
tests = select!(filter(meta) do r
                    r.Subject_ID == "AG2" &&
                    r.RecordSession == "V1" &&
                    r.RecordSite == "ODL3" &&
                    r.ID != "OriSF" &&
                    r.sourceformat == "SpikeGLX"
                    end,
                [:files,:ID,:UUID,:sourceformat])
# show(tests, allrows = true, allcols = true)
files = tests[9,:files]

## Setup Param
param[:model] = [:STA]

# param[:cell] = load(joinpath(param[:resultroot],"cell.jld2"),"cell")
# param[:model] = [:ePPR]
# param[:eppr_ndelay]=1
# param[:eppr_nft]=[3,3,3]
# param[:eppr_lambda]=64000


## Batch Tests
batchtests(tests,param,plot=true)
