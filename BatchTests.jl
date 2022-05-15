includet("Batch.jl")


## Prepare Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "X:/",
    :dataexportroot => "Y:/",
    :resultroot => "Z:/",
    :stimuliroot => "S:/")
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
                    r.ID == "Color" &&
                    r.sourceformat == "SpikeGLX"
                end,
                [:files,:ID,:UUID,:sourceformat])


## Setup Param

files = ".\\AG2\\AG2_V1_ODR1_Flash2Color_1.mat"
pyplot()

param[:spikesorter] = "kilosort3"
param[:model] = [:STA]

# param[:cell] = load(joinpath(param[:resultroot],"cell.jld2"),"cell")
# param[:model] = [:ePPR]
# param[:eppr_ndelay]=1
# param[:eppr_nft]=[3,3,3]
# param[:eppr_lambda]=64000


## Batch Tests
batchtests(tests,param,plot=true)
