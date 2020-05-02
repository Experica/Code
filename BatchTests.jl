## Prepare Param and Metadata
includet("batch.jl")
param = Dict{Any,Any}(
    :dataroot => "../Data",
    :dataexportroot => "../DataExport",
    :resultroot => "../Result",
    :stimuliroot => "../NaturalStimuli")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))

## Query Tests
tests = @from i in meta begin
        @where startswith(get(i.Subject_ID), "AF5")
        @where i.RecordSite == "ODL1"
        @where i.ID == "OriSF"
        @where i.sourceformat == "SpikeGLX"
        @select {i.ID,i.UUID,i.files}
        @collect DataFrame
        end

## Condition Tests
batchtests(tests,param,plot=true)

## HartleySubspace Parametric and Image Response
param[:model]=[:STA]
param[:epprndelay]=1
param[:epprnft]=[3]
param[:epprlambda]=100

batchtests(tests,param,plot=false)
