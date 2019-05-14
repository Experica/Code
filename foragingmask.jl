include("io.jl")
using NeuroAnalysis.NABase

"""
extend equally from a point along an orientation to get a line segment,
then place an orthogonal square drifting grating at the origin and get distance between
grating half cycle line to one end point of the line segment. the second half cycle
of the grating has the same color as the line segment which would occlude the line segment periodically.
argument:
`x`,`y` coordinates of the point
`dir` drifting direction of grating
`l` length of the line segment
`gcw` cycle width of drifting grating
return:
`d` distance between occluding line to the start point of line segment
"""
function pgd(x,y,dir,l,gcw)
  # initial distance from point to grating center line when each presentation started.
  # d>0 indicates grating central line drifting towards point
  d = cos(deg2rad(dir))*y-sin(deg2rad(dir))*x
  # distance from occlusion line to start of line segment
  d-=(gcw+l)/2
end

"""
occlusion cycle by drifting grating
argument:
`gcw` cycle width of drifting grating
`gtf` temporal frequency of drifting grating
`l` length to be occluded
`t` time in ms. od,op=0 when t=0
return:
`od` occlusion distance
`op` occlusion phase
"""
function goc(gcw,gtf,l,t)
  T = 1/gtf;v = gcw/T;t = mod(t*0.001,T);d = v*t
  if d<l
    od=d
    op=(d/l)*0.5pi
  elseif d<gcw/2
    od = l
    op = (d-l)/(gcw/2-l)*0.5pi+0.5pi
  elseif d<gcw/2+l
    od = l-v*(t-T/2)
    op = (d-gcw/2)/l*0.5pi+pi
  else
    od = 0
    op = (d-gcw/2-l)/(gcw/2-l)*0.5pi+1.5pi
  end
  return od,op
end

"""
drifting grating occlusion function for a point, figure array and mask appear at t=0
"""
function got(x,y,dir,l,gcw,gtf,t)
  d = pgd(x,y,dir,l,gcw)
  dt = d/(gcw*gtf)
  goc(gcw,gtf,l,t-dt*1000)
end

"""
occlusion area ratio for different figure type
"""
function got(x,y,dir,l,gcw,gtf,t,ft)
  od,op = got(x,y,dir,l,gcw,gtf,t)
  if ft==0    # Rectangle/square
    or = od/l
  elseif ft==1    # Equilateral Triangle
    htr = od->2*(od/l)^2
    if od>l/2
      or = 1-htr(mod(-od,l/2))
    else
      or = htr(od)
    end
  elseif ft==2    # Circle
    or = nothing
  else
    or = nothing
  end
  return od,op,or
end

"""
drifting grating occlusion time for occluding phase, t=0 when p=0
"""
function gop(gcw,gtf,l,p)
  T = 1/gtf;v = gcw/T
  t1 = 1000(l/v);t2 = 1000((gcw/2-l)/v);T0=t1+t2
  n,p1=fldmod(p,pi)
  if p1<0.5pi
    t=(2p1/pi)*t1+n*T0
  else
    t=t1+(2p1/pi-1)*t2+n*T0
  end
  return t
end

"Prepare factors for all conditions in `DataFrame`"
function condfactor(cond::Dict)
  df = matdictarray2df(cond)
  t = DataArray{Int}(df[:RFFIGIDX]);t[t.==0]=NA;df[:RFFIGIDX]=t
  t = DataArray{Int}(df[:RFFIGTYPE]);t[t.==-1]=NA;df[:RFFIGTYPE]=t
  df[:ORFLIP] = 1-df[:ORFLIP]
  df[:FigSide] = Array{Int}(df[:ORFLIP])
  df[:RevCol] = Array{Int}(df[:REVCOL])
  df[:EdgeContrast] = Array{Int}(df[:ORFLIP] .!= df[:REVCOL])
  df[:FixTarget] = Array{Int}(df[:FIXFIGTYPE] .== df[:REWARDTYPE])
  df[:RFTarget] = DataArray{Int}(df[:RFFIGTYPE] .== df[:REWARDTYPE])
  df[:RFFigType] = df[:RFFIGTYPE]
  return df
end

"Prepare parameters for each condition test"
function condtest(ct1::Dict,cond::DataFrame,param,ct0)
  ct = matdictarray2df(ct1)
  ctcond = cond[ct[:condidx],:]
  ct = [ct ctcond]

  cond0 = param["Condition0"]
  mo = Array{Int}(cond0["MASKOR"][:])
  gcw = cond0["GRCYCLEWID"][:]
  gtf = cond0["GRDRIFTFREQ"][:]
  n0 = length(gtf)
  x = cond0["POSX"]
  x = [x[i,:][:] for i=1:n0]
  y = cond0["POSY"]
  y = [y[i,:][:] for i=1:n0]

  sim = param["SimulateParam"]
  l = [sim["FT0LEN"] for i=1:n0]
  upd = param["SubjectParam"]["UnitsPerDeg"]

  cond0 = DataFrame(Any[mo,gcw,gtf,x,y,l],[:MaskOri,:MaskCycleWidth,:MaskTemporalFreq,:FigPosX,:FigPosY,:FigLength])
  ct = [ct cond0[ct[:condidx0],:]]
  ct[:rffigofun] = map((fi,dir,gcw,gtf,x,y,l,ct0i,ft,ffi,fft)->begin
  if isna(fi)
    ix = (x[ffi]+(x[ffi]-x[ffi-1]))*upd;iy = (y[ffi]+(y[ffi]-y[ffi-1]))*upd
    t->got(ix,iy,dir,l,gcw,gtf,t-ct0["figontime"][:][ct0i],fft)
  else
    t->got(x[fi]*upd,y[fi]*upd,dir,l,gcw,gtf,t-ct0["figontime"][:][ct0i],ft)
  end
end,
ct[:RFFIGIDX],ct[:MaskOri],ct[:MaskCycleWidth],ct[:MaskTemporalFreq],ct[:FigPosX],ct[:FigPosY],ct[:FigLength],ct[:condtestidx0],ct[:RFFIGTYPE],ct[:FIXFIGIDX],ct[:FIXFIGTYPE])
ct[:fixfigofun] = map((fi,dir,gcw,gtf,x,y,l,ct0i,ft)->(t->got(x[fi],y[fi],dir,l,gcw,gtf,t-ct0["figontime"][:][ct0i],ft)),
ct[:FIXFIGIDX],ct[:MaskOri],ct[:MaskCycleWidth],ct[:MaskTemporalFreq],ct[:FigPosX],ct[:FigPosY],ct[:FigLength],ct[:condtestidx0],ct[:FIXFIGTYPE])
ct[:otfun] = map((gcw,gtf,l)->(p->gop(gcw,gtf,l,p)),ct[:MaskCycleWidth],ct[:MaskTemporalFreq],ct[:FigLength])
return ct
end
condtest(ct1::Dict,cond::Dict,param,ct0)=condtest(ct1,condfactor(cond),param,ct0)

"Prepare factors for all condition tests"
function condtestfactor(ct::DataFrame)
  ct0i = ct[:condtestidx0]
  ct0beginend = ct0i[1:end-1] .!= ct0i[2:end]
  t = deepcopy(ct[:FIXFIGIDX][2:end]);t[ct0beginend]=NA;push!(t,NA);ct[:SufFIXFIGIDX] = t
  ct[:ToRFFIG] = DataArray{Int}(ct[:RFFIGIDX] .== ct[:SufFIXFIGIDX])
  return ct
end

function fixseq(ctx)
  ffseq=[];cti=[];ffrepn=[]
  ct0i = convert(Array{Float64,1},unique(ctx[:condtestidx0]))
  for i in ct0i
    ctsofcts0 = find(ctx[:condtestidx0] .== i)
    cffseq = ctx[:FIXFIGIDX][ctsofcts0]
    ffseq = [ffseq;Any[Array{Int}(cffseq)]]
    cti = [cti;Any[ctsofcts0]]
    ffrepn = [ffrepn;(length(cffseq) - length(unique(cffseq)))]
  end
  return ct0i,cti,ffseq,ffrepn
end

function fixdist(ctx,fign)
  fixfigidx = [ctx[:FIXFIGIDX].==i for i=1:fign]
  fix1fig0 = map(x->countnz(x&(ctx[:FixTarget].==1)&(ctx[:FIXFIGTYPE].==0)),fixfigidx)
  fix1fig1 = map(x->countnz(x&(ctx[:FixTarget].==1)&(ctx[:FIXFIGTYPE].==1)),fixfigidx)
  fix0fig0 = map(x->countnz(x&(ctx[:FixTarget].==0)&(ctx[:FIXFIGTYPE].==0)),fixfigidx)
  fix0fig1 = map(x->countnz(x&(ctx[:FixTarget].==0)&(ctx[:FIXFIGTYPE].==1)),fixfigidx)
  return Array{Int}(fix1fig0), Array{Int}(fix1fig1), Array{Int}(fix0fig0), Array{Int}(fix0fig1)
end

function foragingmaskns(datafile;minasfix=15,gminfigfixdur=200,
  factors=["FigSide","EdgeContrast","RevCol","MaskOri","RFTarget","ToRFFIG","RFFigType"])
  data,param = readmat("./data/$datafile")
  ct0 = data["condtests0"]
  minfigfixdur = Int(param["SubjectParam"]["MinFigFixDur"])
  spike = data["cellspike"][1][:]
  ct=condtest(data["condtests"],param["Condition"],param,ct0)
  fixidx = (ct[:figofftime]-ct[:figontime]) .> minasfix
  ctx = condtestfactor(ct[fixidx,:])

  gidx = (ctx[:figofftime]-ctx[:figontime]) .> gminfigfixdur
  dataset = ctx[gidx,[map(symbol,factors);:figontime;:figofftime;:rffigofun;:fixfigofun;:otfun]]
  fl,fln = flfln(dataset,factors)
  ns = Dict("datafile"=>datafile,"minfigfixdur"=>minfigfixdur,"spike"=>spike,
  "dataset"=>dataset,"cellid"=>param["CellID"],"testtype"=>param["TestType"],
  "testrepeat"=>Int(param["TestRepeat"]),"fl"=>fl,"fln"=>fln,"gminfigfixdur"=>gminfigfixdur)
end
