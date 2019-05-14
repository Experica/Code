function f1psthf2(ds,bin,fl,f1,f2;timemark=[0],xmin=-100,xmax=200,ymax=100,align="",f1colors=["dodgerblue","orange"],ldf=[])
  lowlightfun = c -> begin
   t=LCHab(c);RGBA(LCHab(90,20,t.h),0.7)
  end
  df,ss = psth(ds,bin,flcond(fl,f1,f2),spike=symbol("spike$align"))
  df[:f1] = map(x->contains(x,"$(f1)=0")?"$(f1)=0":"$(f1)=1",df[:condition])
  df[:f2] = map(x->contains(x,"$(f2)=0")?"$(f2)=0":"$(f2)=1",df[:condition])
  if align == "off"
    mi=-xmax;xmax=-xmin;xmin=mi
  end
  ps=[]
  dfile = "$(ds[:cellid][1])$(ds[:testtype][1])$(lpad(ds[:testrepeat][1],2,0)).mat"
  t0 = "$dfile - $(f1)-$(f2)=0_$align"
  l=[layer(df[(df[:f2].=="$(f2)=0") & (df[:f1].=="$(f1)=0"),:],x=:x,y=:y,ymin=:ymin,ymax=:ymax,xintercept=timemark,Geom.line,Geom.ribbon,
  Geom.vline(color="gray",size=1pt),Theme(default_color=RGBA(color(f1colors[1]),0.9),lowlight_color=lowlightfun));
    layer(df[(df[:f2].=="$(f2)=0") & (df[:f1].=="$(f1)=1"),:],x=:x,y=:y,ymin=:ymin,ymax=:ymax,Geom.line,Geom.ribbon,
  Theme(default_color=RGBA(color(f1colors[2]),0.9),lowlight_color=lowlightfun))]
  if !isempty(ldf)
    l=[l;layer(ldf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,Geom.line,Geom.ribbon,
               Theme(default_color=RGBA(0.45,0.45,0.45,0.9),lowlight_color=c->RGBA(0.75,0.75,0.75,0.5)))]
  end
  p0=plot(l,Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=0,ymax=ymax),
          Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t0))
  ps=[ps,(p0,t0)]

  t1 = "$dfile - $(f1)-$(f2)=1_$align"
  l=[layer(df[(df[:f2].=="$(f2)=1") & (df[:f1].=="$(f1)=0"),:],x=:x,y=:y,ymin=:ymin,ymax=:ymax,xintercept=timemark,Geom.line,Geom.ribbon,
  Geom.vline(color="gray",size=1pt),Theme(default_color=RGBA(color(f1colors[1]),0.9),lowlight_color=lowlightfun));
    layer(df[(df[:f2].=="$(f2)=1") & (df[:f1].=="$(f1)=1"),:],x=:x,y=:y,ymin=:ymin,ymax=:ymax,Geom.line,Geom.ribbon,
  Theme(default_color=RGBA(color(f1colors[2]),0.9),lowlight_color=lowlightfun))]
  if !isempty(ldf)
    l=[l,layer(ldf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,Geom.line,Geom.ribbon,
               Theme(default_color=RGBA(0.45,0.45,0.45,0.9),lowlight_color=c->RGBA(0.75,0.75,0.75,0.5),lowlight_opacity=1))]
  end
  p1=plot(l,Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=0,ymax=ymax),
          Guide.xlabel("Time (ms)"),Guide.ylabel(""),Guide.title(t1),Guide.manual_color_key("f1",["$(f1)=0","$(f1)=0"],f1colors))
  ps=[ps,(p1,t1)]

  display(hstack(p0,p1))
  return ps
end

function savefigs(ps;path="",width=22cm,height=13cm,dpi=300)
  for p in ps
    NeuroAnalysis.NAVisualization.savefig(p[1],p[2],path=path,format="svg",width=width,height=height,dpi=dpi)
    NeuroAnalysis.NAVisualization.savefig(p[1],p[2],path=path,format="png",width=width,height=height,dpi=dpi)
  end
end

function phasecolor(colorspace::Symbol=:HSV,defaultcolor=RGB(0,0.05,0.15);hues=[],saturation=nothing)
  cm=c->begin
    rs=3;c==1?(qc=4;sc=0.25):(qc=fld(c,0.25)+1;sc=mod(c,0.25))
    if qc==1
      sc=sc*rs+1-0.25*rs
    elseif qc==2
      sc=1-sc*rs
    elseif qc==3
      sc=1-sc*rs
    else
      sc=1-sc*rs
    end
    return qc,sc
  end
  if colorspace==:HSV
    return c->begin
      cs=[240,0,240,120]
      if ~isempty(hues);cs=hues;end
      qc,sc=cm(c)
      if saturation!=nothing;sc=saturation;end
      Colors.HSV(cs[qc],sc,0.9)
    end
  elseif colorspace==:LCHab
    return c->begin
      cs=[280,10,280,150]
      if ~isempty(hues);cs=hues;end
      qc,sc=cm(c)
      if saturation!=nothing;sc=saturation;end
      Colors.LCHab(90,sc*100,cs[qc])
    end
  else
    return c->defaultcolor
  end
end

function gflatrvs(rvs1,l1,rvs2,l2;sortvar1=[],sortvar2=[])
  x1,y1,c1=flatrvs(rvs1,sortvar1);g1=fill(l1,length(x1))
  x2,y2,c2=flatrvs(rvs2,sortvar2);g2=fill(l2,length(x2))
  return [x1;x2],[y1;y2],[c1;c2],[g1;g2]
end

function phaseline(otfun,p0s,ps=-pi:0.5pi:3.5pi;colorfun=phasecolor())
    pts=map(p->map((f,p0)->[f(p)]-f(p0),otfun,p0s),ps)
    pls=map(pl->flatrvs(pl,p0s),pts)
    pcs=Array{RGB}(map(colorfun,mod(ps,2pi)/2pi))
    ls = map((pl,pc)->layer(x=pl[1],y=pl[2],Geom.path,Theme(line_width=1pt,default_color=RGBA(red(pc),green(pc),blue(pc),0.5)))[1],
    pls,pcs)
    cl = Any[["0","0.5\U003C0","\U003C0","1.5\U003C0"],map(colorfun,mod([0,0.5pi,pi,1.5pi],2pi)/2pi)]
    return ls,cl
end

function pstpl(x::Vector,y::Vector,c::Vector,pl=[];timemark=[0],theme=Theme(),xmin=minimum(x)-10,xmax=maximum(x)+10,
              ymin=minimum(y)-1,ymax=maximum(y)+1,colorkey="Phase")
    xl="Time (ms)";yl="Trial Sorted"
    ls = layer(x=x,y=y,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt))
    if ~isempty(pl)
pll = pl[1];cl=pl[2];ls=[ls;pll]
plot(ls,Guide.xlabel(xl),Guide.ylabel(yl),Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
     Guide.manual_color_key(colorkey,cl[1],cl[2]))
else
  plot(ls,Guide.xlabel(xl),Guide.ylabel(yl),Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))
    end
end

function pstpl(ston,pon,stoff,poff,otfun;theme=Theme())
  pstpl(flatrvs(ston,pon)...,phaseline(otfun,pon,colorfun=phasecolor(:HSV,hues=[50,0,240,120],saturation=0.7)),theme=theme,colorkey="RFFig Occluding Phase",colorfun=phasecolor(),colorminv=0,colormaxv=2pi)
  pstpl(flatrvs(stoff,poff)...,phaseline(otfun,poff,colorfun=phasecolor(:HSV,hues=[50,0,240,120],saturation=0.7)),theme=theme,colorkey="RFFig Occluding Phase",colorfun=phasecolor(),colorminv=0,colormaxv=2pi)


end

function pstp(x::Vector,y::Vector,c::Vector;xgroup::Vector=[],timemark=[0],theme=Theme(),colorkey="",colorfun=Scale.lab_gradient(color("white"),color("red")),colorminv=[],colormaxv=[])
  yl="Trial Sorted";ps=25
  if isempty(colorminv);colorminv=minimum(c);end
  if isempty(colormaxv);colormaxv=maximum(c);end
  pline = (x,y,c)->begin
    ti=unique(y);cc=Float64[c[find(y.==i)[1]] for i=ti]
    if abs(minimum(x)) < maximum(x)
      xx=minimum(x)-10-ps*cc
      xmin=minimum(xx)-50;xmax=maximum(x)+50
    else
      xx=maximum(x)+10+ps*cc
      xmin=minimum(x)-50;xmax=maximum(xx)+50
    end
    return xx,ti,cc,xmin,xmax
  end
  pvt = (p,v,t)->begin
    m,i=findmin(abs(p-v))
    return t[i]
  end
  if !isempty(xgroup)
    ug=unique(xgroup)
    xx1,yy1,cc1,xmin,xmax=pline(x[xgroup.==ug[1]],y[xgroup.==ug[1]],c[xgroup.==ug[1]])
    xx2,yy2,cc2,xmin,xmax=pline(x[xgroup.==ug[2]],y[xgroup.==ug[2]],c[xgroup.==ug[2]])
    gg1=fill(ug[1],length(xx1));gg2=fill(ug[2],length(xx2))
    X=[x;xx1;xx2];Y=[y;yy1;yy2];C=[c;cc1;cc2];G=[xgroup;gg1;gg2]
    plot(x=X,y=Y,color=C,xgroup=G,xintercept=fill(timemark[1],length(X)),theme,
         Geom.subplot_grid(Geom.point,Geom.vline(color="gray",size=1pt),free_x_axis=true),
         Guide.xlabel("Time (ms)"),Guide.ylabel(yl),Guide.colorkey(colorkey),
         Scale.ContinuousColorScale(colorfun,minvalue=colorminv,maxvalue=colormaxv))
  else
    xx,yy,cc,xmin,xmax=pline(x,y,c)
    plot(layer(x=x,y=y,color=c,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt)),
         layer(x=xx,y=yy,color=cc,yintercept=[pvt(cc,0.5pi,yy),pvt(cc,pi,yy),pvt(cc,1.5pi,yy)],Geom.hline(color="gray",size=1pt),theme,Geom.point),
         Coord.Cartesian(xmin=xmin,xmax=xmax),Guide.xlabel("Time (ms)"),Guide.ylabel(yl),Guide.colorkey(colorkey),
         Scale.ContinuousColorScale(colorfun,minvalue=colorminv,maxvalue=colormaxv))
  end
end
