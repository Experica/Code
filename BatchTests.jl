# Prepare Env and Metadata
includet("batch.jl")
env = Dict{Any,Any}(
    :dataroot => "../Data",
    :dataexportroot => "../DataExport",
    :resultroot => "../Result",
    :imageroot => "../NaturalStimuli")
meta = readmeta("$(env[:dataexportroot])/metadata.mat")

# Query Test
test = "Color"
tests=querymeta(meta,test=test,sourceformat="SpikeGLX",subject="AF5",recordsite="ODL1")



# Flash of Two Colors
param = copy(env)
r,c = batchtests(tests[5:5,:],param,isplot=true)

# Hartley Subspace Parametric and Image Response
param = copy(env)
param[:model]=[:ePPR]
param[:epprndelay]=1
param[:epprnft]=[3]
param[:epprlambda]=100
ur,uc = batchtests(tests[4:4,:],param,isplot=true);

# Drifting Grating with Ori and SpatialFreq
param = copy(env)
param[:responsedelay] = 0.015
param[:blank] = (:Ori, "blank")
param[:corrcond] = :allcond
r,c = batchtests(tests,param,isplot=true)

# Drifting Grating with Ori, SpatialFreq and Color
param = copy(env)
param[:responsedelay] = 0.015
param[:blank] = (:ColorID, 36)
param[:corrcond] = :allcond
r,c = batchtests(tests,param,isplot=true)

# Color
param = copy(env)
r,c = batchtests(tests,param,isplot=true)



# Save Result
resultdir = joinpath(rroot,testtype)
!isdir(resultdir) && mkpath(resultdir)
save(joinpath(resultdir,"batchresult.jld2"),"ur",ur,"uc",uc)
