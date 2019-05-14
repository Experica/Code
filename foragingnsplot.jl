using Gadfly,Color,DataFrames,NeuroAnalysis.NABase

function foragingnsplot(foragingns;path="./figure/foragingns")
  dfile = foragingns["datafile"]
  #   if !contains(dfile,"29j")
  #     return
  #   end
  spike = foragingns["spike"]
  fl = foragingns["fl"]
  factors = collect(keys(fl))
  minfigfixdur = foragingns["minfigfixdur"]
  gminfigfixdur = foragingns["gminfigfixdur"]

  dataset = foragingns["dataset"]
  figon = foragingns["figon"]
  figoff = foragingns["figoff"]
  fixdur = figoff-figon
  binwidth = 10 # ms
  fixex = 100 # ms
  minrepeat = 5
  ps = []

  isfixdur=true;ispsth=true
  if isfixdur
    t = "$dfile - Longer than $(gminfigfixdur)ms Fixation with Valid Factors"
    p = plot(x=fixdur,Geom.histogram(bincount=100),Coord.Cartesian(xmin=gminfigfixdur,xmax=maximum(fixdur)),
             Guide.xlabel("Fixation Duraion (ms)"),Guide.ylabel("Counts of Fixation"),Guide.title(t))
    ps = [ps,(t,p)]
  end

  if ispsth
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
        print("$dfile: level of \"$f\" do not exist or repeats < $minrepeat.\n")
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
      t = "$dfile - $(f)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
               Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
               Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f1 = "FigSide"
    for f2=filter(x->x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f1,f2))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$dfile: level of \"$f1-$f2\" do not exist or repeats < $minrepeat.\n")
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
      t = "$dfile - $(f1)-$(f2)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
               Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
               Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f2 = "EdgeContrast"
    for f3=filter(x->x!=f2 && x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f2,f3))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$dfile: level of \"$f2-$f3\" do not exist or repeats < $minrepeat.\n")
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
      t = "$dfile - $(f2)-$(f3)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
               Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
               Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f3 = "RFTarget"
    for f4=filter(x->x!=f3 && x!=f2 && x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f3,f4))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$dfile: level of \"$f3-$f4\" do not exist or repeats < $minrepeat.\n")
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
      t = "$dfile - $(f3)-$(f4)"
      p = plot(tdf,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xgroup=:align,xintercept=fill(0,size(tdf,1)),
               Geom.subplot_grid(Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),free_x_axis=true),
               Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"),Guide.title(t))
      ps = [ps,(t,p)]
    end

    f4 = "ToRFFIG"
    for f5=filter(x->x!=f4 && x!=f3 && x!=f2 && x!=f1,factors)
      is,ss = findcond(dataset,flcond(fl,f4,f5))
      if isempty(is) || any(map(length,is).<minrepeat)
        print("$dfile: level of \"$f4-$f5\" do not exist or repeats < $minrepeat.\n")
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
      t = "$dfile - $(f4)-$(f5)"
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
