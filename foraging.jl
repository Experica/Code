include("io.jl")
using NeuroAnalysis.NABase

"Prepare factors for all conditions in `DataFrame`"
function condfactor(cond::Dict)
  df = matdictarray2df(cond)
  # RFFIGIDX == 0:NA when no figure on RF, otherwise RFFIGIDX ∈ {1,2,3,...,10}
  t = DataArray{Int}(df[:RFFIGIDX]);t[t.==0]=NA;df[:RFFIGIDX]=t
  # RFFIGTYPE == -1:NA when no figure on RF, otherwise RFFIGTYPE ∈ {0:square,1:triangle,2:disk,...}
  t = DataArray{Int}(df[:RFFIGTYPE]);t[t.==-1]=NA;df[:RFFIGTYPE]=t
  df[:ORFLIP] = 1-df[:ORFLIP]
  df[:FigSide] = Array{Int}(df[:ORFLIP])
  df[:EdgeContrast] = Array{Int}(df[:ORFLIP] .!= df[:REVCOL])
  # FixTarget ∈ {0:Fix on non-reward figure type,1:Fix on reward figure type}
  df[:FixTarget] = DataArray{Int}(df[:FIXFIGTYPE] .== df[:REWARDTYPE])
  # RFTarget ∈ {0:RF on non-reward figure type,1:RF on reward figure type,NA:no figure on RF}
  df[:RFTarget] = DataArray{Int}(df[:RFFIGTYPE] .== df[:REWARDTYPE])
  df[:RFFigType] = df[:RFFIGTYPE]
  return df
end

"Prepare parameters for each condition test"
function condtest(ct1::Dict,cond::DataFrame,param)
  ct = matdictarray2df(ct1)
  ctcond = cond[ct[:condidx],:]
  ct = [ct ctcond]

  cond0 = param["Condition0"]
  x = cond0["POSX"];y=cond0["POSY"]
  ffi=map((ci,fi)->sub2ind(size(x),ci,Int(fi)),ct[:condidx0],ct[:FIXFIGIDX])
  ct[:FixFigPosX] = x[ffi]
  ct[:FixFigPosY] = y[ffi]
  return ct
end
condtest(ct1::Dict,cond::Dict,param) = condtest(ct1,condfactor(cond),param)

"Prepare factors for all condition tests"
function condtestfactor(ct::DataFrame)
  ct0i = ct[:condtestidx0]
  # There are no sufcond for the end of figurearray and no precond for the begining of figurearray
  ct0beginend = ct0i[1:end-1] .!= ct0i[2:end]
  # SufFIXFIGIDX ∈ {1,2,...10,NA:no sufcondtest}
  t = deepcopy(ct[:FIXFIGIDX][2:end]);t[ct0beginend]=NA;push!(t,NA);ct[:SufFIXFIGIDX] = t
  # SufFixTarget ∈ {0,1,NA:no sufcondtest}
  t = deepcopy(ct[:FixTarget][2:end]);t[ct0beginend]=NA;push!(t,NA);ct[:SufFixTarget] = t
  # SufRFTarget ∈ {0,1,NA:no figure on RF,-1:no sufcondtest}
  t = deepcopy(ct[:RFTarget][2:end]);t[ct0beginend]=-1;push!(t,-1);ct[:SufRFTarget] = t
  # SufRFFIGIDX ∈ {1,2,...,10,NA:no figure on RF,-1:no sufcondtest}
  t = deepcopy(ct[:RFFIGIDX][2:end]);t[ct0beginend]=-1;push!(t,-1);ct[:SufRFFIGIDX] = t
  # ToRFFIG ∈ {0:not saccade to RF,1:saccade to RF,NA:no figure on RF/no sufcondtest}
  ct[:ToRFFIG] = DataArray{Int}(ct[:RFFIGIDX] .== ct[:SufFIXFIGIDX])

  # PreFixTarget ∈ {0,1,NA:no precondtest}
  t = deepcopy(ct[:FixTarget][1:end-1]);t[ct0beginend]=NA;unshift!(t,NA);ct[:PreFixTarget] = t
  # PreRFTarget ∈ {0,1,NA:no figure on RF,-1:no precondtest}
  t = deepcopy(ct[:RFTarget][1:end-1]);t[ct0beginend]=-1;unshift!(t,-1);ct[:PreRFTarget] = t
  # PreRFFIGIDX ∈ {1,2,...,10,NA:no figure on RF,-1:no precondtest}
  t = deepcopy(ct[:RFFIGIDX][1:end-1]);t[ct0beginend]=-1;unshift!(t,-1);ct[:PreRFFIGIDX] = t
  return ct
end

"fixation sequence in each figure array"
function fixseq(ctx)
  ffseq=[];cti=[];ffrepn=[]
  ct0i = unique(ctx[:condtestidx0])
  for i in ct0i
    ctsofcts0 = find(ctx[:condtestidx0] .== i)
    cffseq = ctx[:FIXFIGIDX][ctsofcts0]
    ffseq = [ffseq;Any[Array{Int}(cffseq)]]
    cti = [cti;Any[ctsofcts0]]
    ffrepn = [ffrepn;(length(cffseq) - length(unique(cffseq)))]
  end
  return ct0i,cti,ffseq,ffrepn
end

"fixation distribution on different figures and types"
function fixdist(ctx,fign)
  figtype = unique(ctx[:FIXFIGTYPE])
  fixfigidx = [ctx[:FIXFIGIDX].==i for i=1:fign]
  fix1fig0 = map(x->countnz(x&(ctx[:FixTarget].==1)&(ctx[:FIXFIGTYPE].==figtype[1])),fixfigidx)
  fix1fig1 = map(x->countnz(x&(ctx[:FixTarget].==1)&(ctx[:FIXFIGTYPE].==figtype[2])),fixfigidx)
  fix0fig0 = map(x->countnz(x&(ctx[:FixTarget].==0)&(ctx[:FIXFIGTYPE].==figtype[1])),fixfigidx)
  fix0fig1 = map(x->countnz(x&(ctx[:FixTarget].==0)&(ctx[:FIXFIGTYPE].==figtype[2])),fixfigidx)
  return Array{Int}(fix1fig0), Array{Int}(fix1fig1), Array{Int}(fix0fig0), Array{Int}(fix0fig1)
end

function foragingns(datafile;minasfix=15,gminfigfixdur=200,factors=["FigSide","EdgeContrast","RFTarget","ToRFFIG","RFFigType"])
  data,param = readmat("./data/$datafile")
  fixrad = param["SubjectParam"]["FPWinRad"]/param["SubjectParam"]["UnitsPerDeg"]
  # exclude test when the fixation radius was wrongly set to a huge number
  if fixrad > 1.5
    return nothing
  end

  spike = data["cellspike"][1][:]
  minfigfixdur = Int(param["SubjectParam"]["MinFigFixDur"])
  ct = condtest(data["condtests"],param["Condition"],param)
  fixidx = (ct[:figofftime]-ct[:figontime]) .> minasfix
  ctx = condtestfactor(ct[fixidx,:])

  gidx = (ctx[:figofftime]-ctx[:figontime]) .> gminfigfixdur
  gctidx = !isna(ctx[:ToRFFIG]) & gidx
  gct = ctx[gctidx,:]
  gfigon = Array{Float64}(gct[:figontime])
  gfigoff = Array{Float64}(gct[:figofftime])

  nfgctidx = isna(ctx[:RFFIGIDX]) & gidx
  nfgct = ctx[nfgctidx,:]
  nfgfigon = Array{Float64}(nfgct[:figontime])
  nfgfigoff = Array{Float64}(nfgct[:figofftime])

  dataset = deepcopy(gct[:,map(symbol,factors)])
  fl,fln = flfln(dataset,factors)
  ns = Dict("datafile"=>datafile,"minfigfixdur"=>minfigfixdur,"spike"=>spike,"figon"=>gfigon,"figoff"=>gfigoff,
            "nfigon"=>nfgfigon,"nfigoff"=>nfgfigoff,"dataset"=>dataset,"cellid"=>param["CellID"],"testtype"=>param["TestType"],
            "testrepeat"=>Int(param["TestRepeat"]),"fl"=>fl,"fln"=>fln,"gminfigfixdur"=>gminfigfixdur)
end

function foragingbs(datafile;minasfix=15,gminfigfixdur=200)
  data,param = readmat("./data/$datafile")
  fixrad = param["SubjectParam"]["FPWinRad"]/param["SubjectParam"]["UnitsPerDeg"]
  # exclude test when the fixation radius was wrongly set to a huge number
  if fixrad > 1.5
    return nothing
  end

  ct0 = data["condtests0"]
  fign = Int(param["FigNum"])
  eye = data["eyepoint"]
  minfigfixdur = convert(Int,param["SubjectParam"]["MinFigFixDur"])
  ct=condtest(data["condtests"],param["Condition"],param)
  figon = convert(Array{Float64,1},ct[:figontime])
  figoff = convert(Array{Float64,1},ct[:figofftime])
  figfixdur = figoff-figon

  fixidx = (ct[:figofftime]-ct[:figontime]) .> minasfix
  ctx = condtestfactor(ct[fixidx,:])
  gidx = (ctx[:figofftime]-ctx[:figontime]) .> gminfigfixdur
  gct = ctx[gidx,:]
  gfigon = float64(gct[:figontime])
  gfigoff = float64(gct[:figofftime])
  dataset = deepcopy(gct[:,[:PosX,:PosY]])

  function foragingtaskbs(tminfigfixdur=0)
    tidx = figfixdur .> tminfigfixdur
    tctx= condtestfactor(ct)[tidx,:]
    tfigfixdur=figfixdur[tidx]
    fix1fig0, fix1fig1, fix0fig0, fix0fig1 = fixdist(tctx,fign)
    ct0i,cti,ffseq,ffrepn = fixseq(tctx)

    firsti = int(map(first,cti))
    ffix1fig0, ffix1fig1, ffix0fig0, ffix0fig1 = fixdist(tctx[firsti,:],fign)
    lasti = int(map(last,cti))
    lfix1fig0, lfix1fig1, lfix0fig0, lfix0fig1 = fixdist(tctx[lasti,:],fign)

    hr = countnz(ct0["status"][:][ct0i] .== "Hit")/length(ct0i)
    ffseqn = map(length,ffseq)
    ffrepr = ffrepn./ffseqn

    ffs = (x,y,z,w)->begin
      ff = [x y z w];ffn=sum(ff);ffr = cat(3,ff,ff/ffn);tr=sum([x,y])/ffn
      return ffr,tr
    end

    ffr,tr=ffs(fix1fig0, fix1fig1, fix0fig0, fix0fig1)
    fffr,ftr=ffs(ffix1fig0, ffix1fig1, ffix0fig0, ffix0fig1)
    lffr,ltr=ffs(lfix1fig0, lfix1fig1, lfix0fig0, lfix0fig1)

    bs = Dict("tminfigfixdur"=>tminfigfixdur,"tfigfixdur"=>tfigfixdur,"ffr"=>ffr,"tr"=>tr,"fffr"=>fffr,"ftr"=>ftr,"lffr"=>lffr,"ltr"=>ltr,
              "ct0i"=>ct0i,"cti"=>cti,"ffseq"=>ffseq,"ffrepn"=>ffrepn,"hr"=>hr,"ffseqn"=>ffseqn,"ffrepr"=>ffrepr)
  end

  aft = foragingtaskbs();gft = foragingtaskbs(gminfigfixdur)
  gr = length(gft["tfigfixdur"])/length(aft["tfigfixdur"])

  bs = Dict("datafile"=>datafile,"minfigfixdur"=>minfigfixdur,"dataset"=>dataset,"aft"=>aft,"gft"=>gft,"gr"=>gr,
            "eye"=>eye,"figon"=>gfigon,"figoff"=>gfigoff,"cellid"=>param["CellID"],"testtype"=>param["TestType"],"testrepeat"=>int(param["TestRepeat"]))
end
