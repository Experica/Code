using NeuroAnalysis,StatsBase,JLD2,YAML,DataFrames,Images,HypothesisTests,CairoMakie,AlgebraOfGraphics

function process_fixation(files,param;uuid="",log=nothing,plot=true)
    dataset = prepare(joinpath(param[:dataexportroot],files))

    ex = dataset["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    siteid = join(filter(!isempty,[subject,recordsession,recordsite]),"_")
    siteresultdir = joinpath(param[:resultroot],subject,siteid)
    testid = splitdir(splitext(ex["source"])[1])[2]
    testindex = parse(Int,split(testid,'_')[end])
    resultdir = joinpath(siteresultdir,testid)
    isdir(resultdir) || mkpath(resultdir)

    # test = splitext(splitdir(files)[end])[begin]
    # edata = load(Val(Experica),files)
    # ex = edata["ex"];subject = ex["Subject_ID"];recordsession = ex["RecordSession"];recordsite = ex["RecordSite"]
    # siteid = join(filter!(!isempty,[subject,recordsession,recordsite]),"_")
    # siteresultdir = joinpath(resultroot,subject,siteid)
    # resultdir = joinpath(siteresultdir,test);mkpath(resultdir)

    envparam = ex["EnvParam"];extparam = ex["ExtendParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"]
    exenv=Dict{Any,Any}("ID"=>ex["ID"],"eye"=>ex["Eye"])
    
    ct,ctc = ctctc(ex)
end


param = Dict{Any,Any}(
    :dataroot => "X:/",
    :dataexportroot => "Y:/",
    :resultroot => "Z:/",)
meta = readmeta(joinpath(param[:dataexportroot],"metadata.mat"))


files = ".\\Test\\Test_Full_Fixation_0.mat"
files = ".\\Test\\Test_Full_FixationRFMapping_12.mat"
files = ".\\Test\\Test_Full_FixationDrawBorder_8.mat"

ct,ctc = process_fixation(files,param)

ctd =  dropmissing(ct)
hr = round(count(i->i=="HIT", ctd.TaskResult)/nrow(ctd)*100,digits=1)
spec = data(ctd) * mapping(:TaskResult,color=:TaskResult ) * frequency()
draw(spec;figure=(;title="Hit Rate = $(hr)%"))
