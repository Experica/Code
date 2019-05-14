using NeuroAnalysis.NACore, NeuroAnalysis.NABase

function tooseq(wtranss::Vector,oos::Vector,otranss::Vector,orots::Vector;wo=Vec4())
  seqn = length(wtranss)
  if seqn == length(oos) == length(otranss) == length(orots) != 0
    toos = Array(Vec4{typeof(oos[1].x)},seqn+1)
    toos[1] = wo
    for i in 1:seqn
      too = translate(transrotatez(oos[i],otranss[i],orots[i]),wtranss[i])
      toos[i+1] = transform(too,translation(toos[i].x,toos[i].y,toos[i].z))
    end
    fti=[1:seqn+1]+1;fti[end]=0
    return toos,fti
  else
    error("Lengths of Argument Vector Don't Equal.")
  end
end
function tooseq(wtrans::Vec3,otranss::Vector,orot::Real;wo=Vec4())
  seqn = length(otranss)
  tooseq(fill(wtrans,seqn),fill(Vec4(),seqn),otranss,fill(orot,seqn),wo=wo)
end
function figarray343(wtrans::Vec3,otranss::Vector,orot::Real;linespace=1.0)
  ctooseq,cfti = tooseq(wtrans,otranss[2:4],orot)
  ltooseq,lfti = tooseq(wtrans,otranss[6:7],orot)
  rtooseq,rfti = tooseq(wtrans,otranss[9:10],orot)
  lfti += 4;lfti[end]=0
  rfti += 7;rfti[end]=0

  dv = [diff(ctooseq),diff(ltooseq),diff(rtooseq)]
  ml = linespace*mean(map(length,dv));md = -1*norm(sum(dv))
  ctooseq = translate(ctooseq,convert(Vec3,(ml/2)*md))
  ltooseq = translate(ltooseq,convert(Vec3,rotatexyz(ml*md,Vec3(0.0,0.0,pi/2))))
  rtooseq = translate(rtooseq,convert(Vec3,rotatexyz(ml*md,Vec3(0.0,0.0,-pi/2))))
  toos = translate([ctooseq,ltooseq,rtooseq],convert(Vec3,ml*md))
  return convert(Array{Vec3},toos),[cfti,lfti,rfti],ml
end
function figarray343(wtrans::Vec3,otranss::Vector,otransi::Vector,orot::Real;linespace=1.0)
  figarray343(wtrans,map(i->otranss[i],otransi),orot,linespace=linespace)
end

function randfp(vs::Vector,vsr,dismin;radmin=0.5)
  itermax = 100;iterscalestart=50
  vsrscale = 1.0;vsrscalemax=2.0
  vsrscalestep=(vsrscalemax-vsrscale)/(itermax-iterscalestart)
  for i=1:itermax
    if i>iterscalestart
      vsrscale += vsrscalestep
    end
    fp=vsrscale*vsr*randvec3(radmin,is2d=true)
    dis = map(x->length(fp-x),vs)
    if all(x->x>=dismin,dis)
      return fp
    end
  end
  error("Can't find valid fixation point.")
end
function rewardfig(figseq,rftype,rn=1)
  rfti = shuffle!(find(figseq.==rftype))
  rf = falses(length(figseq))
  rf[rfti[1:rn]]=true
  return rf
end
