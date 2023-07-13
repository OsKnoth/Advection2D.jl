function Lagrange(x,xP,iP)
  f = 1.0
  for i = 1 : size(xP,1)
     if i != iP
        f = f * (x - xP[i]) / (xP[iP] - xP[i])
     end
  end
  return f
end
