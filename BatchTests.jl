includet("Batch.jl")

## Prepare Param and Metadata
cd(@__DIR__)
param = Dict{Any,Any}(
    :dataroot => "../Data",
    :dataexportroot => "../DataExport",
    :resultroot => "../Result",
    :stimuliroot => "../NaturalStimuli")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))

param[:layer] = layer

## Query Tests
tests = @from i in meta begin
        @where startswith(get(i.Subject_ID), "AF5")
        @where i.RecordSite == "ODL2"
        @where i.ID == "HartleySubspace"
        @where i.sourceformat == "SpikeGLX"
        @select {i.ID,i.UUID,i.files,i.sourceformat}
        @collect DataFrame
        end

## Batch Condition Tests
batchtests(tests,param,plot=true)

## HartleySubspace Parametric and Image Response
param[:model]=[:STA]
param[:blank] = (:Ori_Final,NaN)

param[:epprndelay]=1
param[:epprnft]=[3]
param[:epprlambda]=100

batchtests(tests,param,plot=false)
