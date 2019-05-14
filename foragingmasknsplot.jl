include("foragingmask.jl")
using Gadfly,Colors,NeuroAnalysis.NAVisualization
include("plot.jl")

function foragingmasknsplot(datafile;area="",minasfix=15,gminfigfixdur=200,isst=true,ispsth=false,fixex=100,
                            factors=["FigSide","EdgeContrast","RevCol","MaskOri","RFTarget","ToRFFIG","RFFigType"],
                            binwidth=10,path="./figure/foragingmaskns",mintrial=50,minrepeat=5)
  data,param = readmat("./data/$datafile")
  ct0 = data["condtests0"]
  minfigfixdur = Int(param["SubjectParam"]["MinFigFixDur"])
  spike = data["cellspike"][1][:]
  ct=condtest(data["condtests"],param["Condition"],param,ct0)
  fixidx = (ct[:figofftime]-ct[:figontime]) .> minasfix
  ctx = condtestfactor(ct[fixidx,:])
  gidx = (ctx[:figofftime]-ctx[:figontime]) .> gminfigfixdur
  gct = ctx[gidx,:]

  ps = []
  if isst
    dataset=gct[!isna(gct[:RFFIGIDX]),:]
    figon = Array{Float64}(dataset[:figontime])
    figoff = Array{Float64}(dataset[:figofftime])
    st,sn,ws,is = subrv(spike,figon-fixex,figon+gminfigfixdur+fixex,isminzero=true)
    ston = map(x->x-fixex,st)
    st,sn,ws,is = subrv(spike,figoff-gminfigfixdur-fixex,figoff+fixex,ismaxzero=true)
    stoff = map(x->x+fixex,st)
    pon=map((f,x)->f(x)[2],dataset[:rffigofun],figon)
    poff=map((f,x)->f(x)[2],dataset[:rffigofun],figoff)

    cdataset = gct[isna(gct[:RFFIGIDX]),:]
    cfigon = Array{Float64}(cdataset[:figontime])
    cfigoff = Array{Float64}(cdataset[:figofftime])
    st,sn,ws,is = subrv(spike,cfigon-fixex,cfigon+gminfigfixdur+fixex,isminzero=true)
    cston = map(x->x-fixex,st)
    st,sn,ws,is = subrv(spike,cfigoff-gminfigfixdur-fixex,cfigoff+fixex,ismaxzero=true)
    cstoff = map(x->x+fixex,st)
    cpon=map((f,x)->f(x)[2],cdataset[:rffigofun],cfigon)
    cpoff=map((f,x)->f(x)[2],cdataset[:rffigofun],cfigoff)


    title="spiketrain"
    sttheme = Theme(default_point_size=1.2pt,default_color=RGBA(0,0.2,0.5,0.7),continuous_highlight_color=c->nothing,highlight_width=0mm)
    t=gflatrvs(ston,"On",stoff,"Off")
    stn = length(unique(t[2]))
    if stn >= mintrial
      p = plotspiketrain(t[1:3]...,xgroup=t[4],theme=sttheme)
      ps = [ps;("$datafile - $title - $area",p)]

      p = pstpl(flatrvs(ston,pon)...,phaseline(dataset[:otfun],pon,colorfun=phasecolor(hues=[50,0,240,120],saturation=0.7)),theme=sttheme)
      ps = [ps;("$datafile - $(title)-phase - $area",p)]

      p = pstpl(flatrvs(cston,cpon)...,phaseline(cdataset[:otfun],cpon,colorfun=phasecolor(hues=[50,0,240,120],saturation=0.7)),theme=sttheme)
      ps = [ps;("$datafile - $(title)-phase-nofig - $area",p)]
    end
  end

  if ispsth
    fl,fln = flfln(dataset,factors)
    st,sn,ws,is = subrv(spike,figon-fixex,figon+gminfigfixdur+fixex,isminzero=true)
    ston = map(x->x-fixex,st)
    binon = -fixex:binwidth:gminfigfixdur
    st,sn,ws,is = subrv(spike,figoff-gminfigfixdur-fixex,figoff+fixex,ismaxzero=true)
    stoff = map(x->x+fixex,st)
    binoff = -gminfigfixdur:binwidth:fixex
    st={ston,stoff};bin={binon,binoff};phase={phaseon,phaseoff};align=["on","off"]



    st,sn,ws,is = subrv(spike,figon-fixex,figon+gminfigfixdur+fixex,isminzero=true)
    ston = map(x->x-fixex,st)
    binon = -fixex:binwidth:gminfigfixdur
    st,sn,ws,is = subrv(spike,figoff-gminfigfixdur-fixex,figoff+fixex,ismaxzero=true)
    stoff = map(x->x+fixex,st)
    binoff = -gminfigfixdur:binwidth:fixex
    st={ston,stoff};bin={binon,binoff};align=["on","off"]

    for f=keys(fl)
      is,ss = findcond(dataset,flcond(fl,f))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$datafile: level of \"$f\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      tdf=DataFrame()
      for i = 1:2
        df = psth(map(x->st[i][x],is),bin[i],ss)
        df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
        df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
        df[:align]=fill(align[i],size(df,1))
        tdf = [tdf,df]
      end
      t = "$datafile - $(f)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
      Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
      Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f1 = "FigSide"
    for f2=filter(x->x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f1,f2))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$datafile: level of \"$f1-$f2\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      tdf=DataFrame()
      for i = 1:2
        df = psth(map(x->st[i][x],is),bin[i],ss)
        df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
        df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
        df[:align]=fill(align[i],size(df,1))
        tdf = [tdf,df]
      end
      t = "$datafile - $(f1)-$(f2)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
      Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
      Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f2 = "EdgeContrast"
    for f3=filter(x->x!=f2 && x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f2,f3))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$datafile: level of \"$f2-$f3\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      tdf=DataFrame()
      for i = 1:2
        df = psth(map(x->st[i][x],is),bin[i],ss)
        df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
        df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
        df[:align]=fill(align[i],size(df,1))
        tdf = [tdf,df]
      end
      t = "$datafile - $(f2)-$(f3)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
      Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
      Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f3 = "RFTarget"
    for f4=filter(x->x!=f3 && x!=f2 && x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f3,f4))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$datafile: level of \"$f3-$f4\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      tdf=DataFrame()
      for i = 1:2
        df = psth(map(x->st[i][x],is),bin[i],ss)
        df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
        df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
        df[:align]=fill(align[i],size(df,1))
        tdf = [tdf,df]
      end
      t = "$datafile - $(f3)-$(f4)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
      Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
      Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f4 = "ToRFFIG"
    for f5=filter(x->x!=f4 && x!=f3 && x!=f2 && x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f4,f5))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$datafile: level of \"$f4-$f5\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      tdf=DataFrame()
      for i = 1:2
        df = psth(map(x->st[i][x],is),bin[i],ss)
        df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
        df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
        df[:align]=fill(align[i],size(df,1))
        tdf = [tdf,df]
      end
      t = "$datafile - $(f4)-$(f5)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
      Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
      Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end
  end

  if !isempty(ps)
    for p in ps
      NeuroAnalysis.NAVisualization.savefig(p[2],p[1],path=path,format="svg",width=30cm)
      NeuroAnalysis.NAVisualization.savefig(p[2],p[1],path=path,format="png",width=30cm)
    end
  end
end
