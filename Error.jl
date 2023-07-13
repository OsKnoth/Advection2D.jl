function Error!(e,c,dx)
#  c_0     c_1     c_2      c_n
#      dx_1    dx_2     dx_n

n = size(c,1)
@inbounds for i = 2 : n - 1
  e[i-1] = ErrorLoc(c[i-1],c[i],c[i+1],dx[i-1],dx[i],.2)
end
end

@inline function ErrorLoc(cm,c,cp,dxm,dxp,d)
  eLoc = abs(dxm*cp-(dxm+dxp)*c+dxp*cm) /
    (dxm*abs(cp-cm) + dxp*abs(c-cm) +
    d*(dxm*abs(cp)+(dxm+dxp)*abs(c)+dxp*abs(cm))+1.e-14);
end
