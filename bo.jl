include("io.jl")
using NeuroAnalysis.NABase

function condfactor(cond::Dict)
  df = matdictarray2df(cond)
  df[:FigSide] = Array{Int}(df[:OFLIP])
  df[:EdgeContrast] = Array{Int}(df[:REVCOL])
  df[:xRFSize] = Array{Int}(df[:RFSIZ2WID])
  df[:FigType] = df[:Shape]
  delete!(df,[:AUTOOVERLAP,:NUMFIGS,:OFLIP,:RFSIZ2WID,:RIDGE,:RIDGE10,:Shape,:REVCOL])
  if any(names(df).==:OVERLAP);delete!(df,:OVERLAP);end
  return df
end

function condtest(ct1::Dict,cond::DataFrame)
  ks = ["condidx","condrepeat","status","figontime","figofftime","targetontime","fixontime","testofftime"]
  ct = DataFrame(Any[squeeze(v,2) for v in map(k->ct1[k],ks)],[Symbol(k) for k in ks])
  vi = Array{Bool}(map(x->!isempty(x),ct[:figontime])) & Array{Bool}(map(x->!isempty(x),ct[:figofftime]))
  vct = ct[vi,:]
  vctcond = cond[vct[:condidx],:]
  [vct vctcond]
end
condtest(ct1::Dict,cond::Dict) = condtest(ct1,condfactor(cond))

function bons(datafile;factors=["FigSide","EdgeContrast","xRFSize","FigType"])
  data,param = readmat("./data/$datafile")
  fixrad = param["SubjectParam"]["FPWinRad"]/param["SubjectParam"]["UnitsPerDeg"]
  # exclude test when the fixation radius was wrongly set to a huge number
  if fixrad > 1.5
    return nothing
  end

  spike = data["cellspike"][1][:]
  minconddur = Int(param["SubjectParam"]["MinCondDur"])
  ct = condtest(data["condtests"],param["Condition"])

  gct = ct[ct[:status] .!= "Early",:]
  gfigon = Array{Float64}(gct[:figontime])
  gfigoff = Array{Float64}(gct[:figofftime])
  dataset = deepcopy(gct[:,map(symbol,factors)])

  fl,fln = flfln(dataset,factors)
  ns = Dict("datafile"=>datafile,"minconddur"=>minconddur,"spike"=>spike,"figon"=>gfigon,"figoff"=>gfigoff,
            "dataset"=>dataset,"cellid"=>param["CellID"],"testtype"=>param["TestType"],
            "testrepeat"=>Int(param["TestRepeat"]),"fl"=>fl,"fln"=>fln)
end

function bobs(datafile;gminfixdur=200,factors=["FigSide","EdgeContrast","xRFSize","FigType"])
  data,param = readmat("./data/$datafile")
  fixrad = param["SubjectParam"]["FPWinRad"]/param["SubjectParam"]["UnitsPerDeg"]
  # exclude test when the fixation radius was wrongly set to a huge number
  if fixrad > 1.5
    return nothing
  end

  eye = data["eyepoint"]
  minconddur = Int(param["SubjectParam"]["MinCondDur"])
  ct = condtest(data["condtests"],param["Condition"])

  ti = find(Array{Bool}(map(x->!isempty(x),ct[:targetontime])))
  tgton = Array{Float64}(ct[ti[1:end-1],:targetontime])
  fixon = Array{Float64}(ct[ti[1:end-1],:fixontime])
  fixoff = Array{Float64}(ct[ti[2:end]-1,:testofftime])
  gidx = (fixoff-fixon) .> gminfixdur
  tgton = tgton[gidx];fixon = fixon[gidx];fixoff = fixoff[gidx]

  gi = ct[:status] .!= "Early"
  figon = Array{Float64}(ct[gi,:figontime])
  figoff = Array{Float64}(ct[gi,:figofftime])
  dataset = deepcopy(ct[gi,map(symbol,factors)])

  fl,fln = flfln(dataset,factors)
  bs = Dict("datafile"=>datafile,"minconddur"=>minconddur,"eye"=>eye,"figon"=>figon,"figoff"=>figoff,
            "tgton"=>tgton,"fixon"=>fixon,"fixoff"=>fixoff,
            "dataset"=>dataset,"cellid"=>param["CellID"],"testtype"=>param["TestType"],
            "testrepeat"=>Int(param["TestRepeat"]),"fl"=>fl,"fln"=>fln)
end
