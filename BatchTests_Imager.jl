includet("Batch/Batch.jl")


## Param and Metadata
param = Dict{Any,Any}(
    :dataroot => "D:/",
    :dataexportroot => "Y:/",
    :resultroot => "Z:/")
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))


## Query Tests
tests = select!(filter(meta) do r
                    r.Subject_ID in ["Test", "AG5"] &&
                    # r.RecordSession == "V1" &&
                    # r.RecordSite == "Full" &&
                    # r.ID in ["ISIEpochFlash2Color"] &&
                    r.sourceformat == "Imager"
                end,
                [:files,:ID,:UUID,:sourceformat])


## Setup Param
files = ".\\AG1\\AG1_V1V2_Full_ISICycleOri_3.mat"
files = ".\\AG1\\AG1_V1V2_Full_ISICycleColorPlane_0.mat"
files = ".\\AG1\\AG1_V1V2_Full_ISIEpochOri8_2.mat"
files = ".\\AG1\\AG1_V1V2_Full_ISICycle2Color_1.mat"
files = ".\\Test\\Test_Full_ISIEpoch2Color_0.mat"


## Batch Tests
batchtests(tests,param,plot=true)