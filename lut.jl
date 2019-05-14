
using Images
rg=0:0.05:1
lut=[RGBA(r,g,b,1) for g=reverse(rg),r=rg,b=rg]

##
 clut=cat(map(i->lut[:,:,i],1:21)...,dims=2)
##
using Plots
plot(sin)
