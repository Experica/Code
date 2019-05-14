function cvformula(p,mf)
  x,i = findmin(p.meanloss)
  vti = find(p.path.betas[2:end,i])
  vt = mf.terms.terms[vti]
  rt = mf.terms.eterms[1]
  f = Formula(rt,Expr(:call,:+,vt...))
end

function showmodelpath(p)
  pn = length(p.lambda)
  tn = size(p.path.betas,1)
  minloss,i = findmin(p.meanloss)
  x=1:pn
  ps = Layer[]
  for t =1:tn
    y = p.path.betas[t,:]
    push!(ps,layer(y=y,x=x,Geom.line)[1])
  end
  display(plot(ps,Guide.xlabel("Model Path"),Guide.ylabel("Coeffients")))
  plot(y=p.meanloss,x=1:pn,xintercept=[i],Geom.line,Geom.vline(color=colorant"red"),
       Guide.xlabel("Model Path"),Guide.ylabel("Cross-Validation Meanloss"))
end
