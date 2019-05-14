using Gadfly,Color,NeuroAnalysis.NABase

function bonsplot(bons;path="./figure/bons")
  dfile = bons["datafile"]
  #   if !contains(dfile,"29j")
  #     return
  #   end
  spike = bons["spike"]
  fl = bons["fl"]
  factors = collect(keys(fl))
  minconddur = bons["minconddur"]

  dataset = bons["dataset"]
  figon = bons["figon"]
  figoff = bons["figoff"]
  binwidth = 10 # ms
  figonex = 100 # ms
  figoffex = 200 # ms
  minrepeat = 5
  ps = []

  ispsth=true
  if ispsth
    st,sn,ws,is = subrv(spike,figon-figonex,figoff+figoffex,isminzero=true)
    st = map(x->x-figonex,st)
    bin = -figonex:binwidth:minconddur+figoffex
    onoff=[0,minconddur]

    for f=keys(fl)
      is,ss = findcond(dataset,flcond(fl,f))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$dfile: level of \"$f\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      df = psth(map(x->st[x],is),bin,ss)
      df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
      df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
      t = "$dfile - $(f)"
      p = plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=onoff,
               Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
               Coord.Cartesian(xmin=bin[1],xmax=bin[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f1 = "FigSide"
    for f2=filter(x->x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f1,f2))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$dfile: level of \"$f1-$f2\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      df = psth(map(x->st[x],is),bin,ss)
      df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
      df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
      t = "$dfile - $(f1)-$(f2)"
      p = plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=onoff,
               Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
               Coord.Cartesian(xmin=bin[1],xmax=bin[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f2 = "EdgeContrast"
    for f3=filter(x->x!=f2 && x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f2,f3))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$dfile: level of \"$f2-$f3\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      df = psth(map(x->st[x],is),bin,ss)
      df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
      df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
      t = "$dfile - $(f2)-$(f3)"
      p = plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=onoff,
               Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
               Coord.Cartesian(xmin=bin[1],xmax=bin[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f3 = "xRFSize"
    for f4=filter(x->x!=f3 && x!=f2 && x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f3,f4))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$dfile: level of \"$f3-$f4\" do not exist or repeats < $minrepeat.\n")
        continue
      end
      df = psth(map(x->st[x],is),bin,ss)
      df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
      df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
      t = "$dfile - $(f3)-$(f4)"
      p = plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=onoff,
               Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
               Coord.Cartesian(xmin=bin[1],xmax=bin[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

  end
  if !isempty(ps)
    for p in ps
      NeuroAnalysis.NAVisualization.savefig(p[2],p[1],path=path,format="svg")
      NeuroAnalysis.NAVisualization.savefig(p[2],p[1],path=path,format="png")
    end
  end
end
