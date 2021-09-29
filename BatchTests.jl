includet("Batch.jl")

## Prepare Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "X:/",
    :dataexportroot => "Y:/",
    :resultroot => "Z:/",
    :stimuliroot => "X:/NaturalStimuli")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))

## Query Tests
tests = select!(filter(meta) do r
                    r.Subject_ID == "AG2" &&
                    r.RecordSession == "V1" &&
                    r.RecordSite == "ODR2" &&
                    r.ID == "HartleySubspace" &&
                    r.sourceformat == "SpikeGLX"
                    end,
                [:ID,:files,:UUID,:sourceformat])

## Batch ISI Tests
batchtests(tests[3:6,:],param)



## Batch Condition Tests
batchtests(tests,param,plot=true)
## Batch HartleySubspace Tests
param[:model] = [:STA]
batchtests(tests,param,plot=true)



param[:cell] = load(joinpath(param[:resultroot],"cell.jld2"),"cell")
param[:model] = [:ePPR]
param[:eppr_ndelay]=1
param[:eppr_nft]=[3,3,3]
param[:eppr_lambda]=64000

batchtests(tests[1:1,:],param,plot=true)
