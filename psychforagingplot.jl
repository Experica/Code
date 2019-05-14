include("foraging.jl")
using PyPlot,Color

function psychforagingplot(dfname)
  if !contains(dfname,"29j")
    return
  end
  showffseq = false;showffdist=true;showffrep=true;showffdur=true

  data,param = readmat("./data/$dfname")
  cts0 = data["condtests0"]
  status0 = cts0["status"][:]
  cts0n = length(status0)
  fign = param["FigNum"]
  minfigfixdur = convert(Int,param["SubjectParam"]["MinFigFixDur"])
  ctscondex=ctscondfactorextend(ctscondextend(data["condtests"],factorextend(param["Condition"])));
  figon = convert(Array{Float64,1},ctscondex[:figontime])
  figoff = convert(Array{Float64,1},ctscondex[:figofftime])
  figfixdur = figoff-figon;

  PyPlot.svg(true);nc=10
  pygui(false)
  cmb = colormap("Blues",25)[14:23]
  cmo = colormap("Oranges",19)[19-nc+1:end]
  cmg = colormap("Greens",5)[3:4]
  cmr = colormap("Reds",5)[3:4]
  ft0 = "Rectangle";ft1 = "Triangle"
  ps=[];dpi=5000;fs=9;ptoi=1/72;nbins=100
  gidx = ctscondex[:status] .!= "Early"
  for i=["All Fixation","Longer than $(minfigfixdur)ms Fixation"]
    if i=="All Fixation"
      cctscond = ctscondex
      fixdur = figfixdur
    else
      cctscond = ctscondex[gidx,:]
      fixdur = figfixdur[gidx]
    end
    binwidth = maximum(fixdur)/nbins
    binedge = 0:binwidth:maximum(fixdur)+2binwidth

    if showffdur
      p=figure(figsize=(12,6),dpi=dpi)
      dns,wins,subs,sis = histtps(fixdur,binedge)
      bar(binedge[1:end-1],dns,color=[cmb[6].r, cmb[6].g, cmb[6].b],width=binwidth,linewidth=0.5)
      xlim(binedge[1]-2binwidth,binedge[end])
      t = "$dfname - $i Duration"
      title(t);xlabel("Fixation Duration (ms)");ylabel("Counts of Fixation")
      ps=[ps,(t,p)]
    end

    fix1fig0, fix1fig1, fix0fig0, fix0fig1 = fixdist(cctscond,fign)
    cts0idx,ffrepn,ffseq,cctscondi = fixseq(cctscond)
    if showffdist
      p=figure(figsize=(11,6),dpi=dpi)
      f1f0 = bar(1:fign,fix1fig0,color=[cmg[2].r, cmg[2].g, cmg[2].b],align="center",linewidth=0)
      f1f1 = bar(1:fign,fix1fig1,bottom=fix1fig0,color=[cmg[1].r, cmg[1].g, cmg[1].b],align="center",linewidth=0)
      f0f0 = bar(1:fign,fix0fig0,bottom=fix1fig0+fix1fig1,color=[cmr[2].r, cmr[2].g, cmr[2].b],align="center",linewidth=0)
      f0f1 = bar(1:fign,fix0fig1,bottom=fix1fig0+fix1fig1+fix0fig0,color=[cmr[1].r, cmr[1].g, cmr[1].b],align="center",linewidth=0)
      xlim(0,11)
      xticks(1:fign)
      legend((f1f0[1],f1f1[1],f0f0[1],f0f1[1]),("FixOnTarget & FigType_$ft0","FixOnTarget & FigType_$ft1",
                                                "FixOffTarget & FigType_$ft0","FixOffTarget & FigType_$ft1"),
             bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0,fontsize="small",frameon=false)
      t = "$dfname - $i"
      title(t);xlabel("Figure");ylabel("Counts of Fixation")
      ps=[ps,(t,p)]
    end

    if showffseq
      x=[];y=[];m=[];d=[];ymax=1
      for t=1:length(cts0idx)
        x = [x,fill(cts0idx[t],length(ffseq[t]))]
        y = [y,[1:length(ffseq[t])]]
        m = [m,convert(Array{Int},ffseq[t])]
        d = [d, fixdur[cctscondi[t]]]
        ymax=maximum([ymax,length(ffseq[t])])
      end
      sfw = maximum([5, fs*ptoi*(cts0n+2)])
      sfh = maximum([3, fs*ptoi*(ymax+2)])
      p=figure(figsize=(sfw,sfh),dpi=dpi)
      #scatter(x,y)
      xlim(0,cts0n+1)
      ylim(0,ymax+1)
      xn = length(x)
      for j=1:xn
        if d[j] >= minfigfixdur
          ci = floor((d[j]-minfigfixdur)/((minfigfixdur+400)/nc))+1
          ci=minimum([nc,ci])
          c = [cmo[ci].r, cmo[ci].g, cmo[ci].b]
        else
          ci = floor(d[j]/(minfigfixdur/nc))+1
          ci=minimum([nc,ci])
          c = [cmb[ci].r, cmb[ci].g, cmb[ci].b]
        end
        ((j+1 > xn) || (x[j] != x[j+1])) && (status0[x[j]] == "Hit") && (c = [cmr[2].r, cmr[2].g, cmr[2].b])
        annotate(string(m[j]-1),[x[j], y[j]],ha="center",va="center",size=fs,color=c)
      end
      t = "$dfname - $i Sequence in Trial"
      title(t);ylabel("Sequence");xlabel("FigArray Trial")
      ps = [ps,(t,p)]
    end

    if showffdist
      firsti = int(map(x->x[1],cctscondi))
      fix1fig0, fix1fig1, fix0fig0, fix0fig1 = fixdist(cctscond[firsti,:],fign)
      p=figure(figsize=(11,6),dpi=dpi)
      f1f0 = bar(1:fign,fix1fig0,color=[cmg[2].r, cmg[2].g, cmg[2].b],align="center",linewidth=0)
      f1f1 = bar(1:fign,fix1fig1,bottom=fix1fig0,color=[cmg[1].r, cmg[1].g, cmg[1].b],align="center",linewidth=0)
      f0f0 = bar(1:fign,fix0fig0,bottom=fix1fig0+fix1fig1,color=[cmr[2].r, cmr[2].g, cmr[2].b],align="center",linewidth=0)
      f0f1 = bar(1:fign,fix0fig1,bottom=fix1fig0+fix1fig1+fix0fig0,color=[cmr[1].r, cmr[1].g, cmr[1].b],align="center",linewidth=0)
      xlim(0,11)
      xticks(1:fign)
      legend((f1f0[1],f1f1[1],f0f0[1],f0f1[1]),("FixOnTarget & FigType_$ft0","FixOnTarget & FigType_$ft1",
                                                "FixOffTarget & FigType_$ft0","FixOffTarget & FigType_$ft1"),
             bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0,fontsize="small",frameon=false)
      t="$dfname - First Fixation in FigArray Trial($i)"
      title(t);xlabel("Figure");ylabel("Counts of Fixation")
      ps=[ps,(t,p)]

      lasti = int(map(x->x[end],cctscondi))
      fix1fig0, fix1fig1, fix0fig0, fix0fig1 = fixdist(cctscond[lasti,:],fign)
      p=figure(figsize=(11,6),dpi=dpi)
      f1f0 = bar(1:fign,fix1fig0,color=[cmg[2].r, cmg[2].g, cmg[2].b],align="center",linewidth=0)
      f1f1 = bar(1:fign,fix1fig1,bottom=fix1fig0,color=[cmg[1].r, cmg[1].g, cmg[1].b],align="center",linewidth=0)
      f0f0 = bar(1:fign,fix0fig0,bottom=fix1fig0+fix1fig1,color=[cmr[2].r, cmr[2].g, cmr[2].b],align="center",linewidth=0)
      f0f1 = bar(1:fign,fix0fig1,bottom=fix1fig0+fix1fig1+fix0fig0,color=[cmr[1].r, cmr[1].g, cmr[1].b],align="center",linewidth=0)
      xlim(0,11)
      xticks(1:fign)
      legend((f1f0[1],f1f1[1],f0f0[1],f0f1[1]),("FixOnTarget & FigType_$ft0","FixOnTarget & FigType_$ft1",
                                                "FixOffTarget & FigType_$ft0","FixOffTarget & FigType_$ft1"),
             bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0,fontsize="small",frameon=false)
      t="$dfname - Last Fixation in FigArray Trial($i)"
      title(t);xlabel("Figure");ylabel("Counts of Fixation")
      ps=[ps,(t,p)]
    end

    if showffrep
      p=figure(figsize=(12,6),dpi=dpi)
      plot(cts0idx,ffrepn,aa=true,color=[cmb[6].r, cmb[6].g, cmb[6].b],lw=1)
      xlim(0,cts0n+1)
      t="$dfname - $i in Trial"
      title(t);ylabel("Repetition");xlabel("FigArray Trial")
      ps=[ps,(t,p)]
    end
  end
  figpath="./figure/psychforaging"
  for p in ps
    NeuroAnalysis.NAVisualization.savefig(p[2],p[1],path=figpath)
  end
end
