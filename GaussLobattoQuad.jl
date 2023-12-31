function GaussLobattoQuad(n)
  w = zeros(n+1)
  x = zeros(n+1)
  if n == 1
    w[1]=1
    w[2]=1
    x[1]=-1
    x[2]=1
      
  elseif n == 2
    w[2]=4/3
    w[1]=1/3
    w[3]=1/3
      
    x[2]=0
    x[1]=-1
    x[3]=1
      
  elseif n == 3
    w[2]=5/6
    w[3]=5/6
    w[1]=1/6
    w[4]=1/6
      
    x[2]=-1/sqrt(5)
    x[3]=1/sqrt(5)
    x[1]=-1
    x[4]=1
      
  elseif n == 4
    w[3]=32/45
    w[2]=49/90
    w[4]=49/90
    w[1]=1/10
    w[5]=1/10
    
    x[3]=0
    x[2]=-sqrt(3/7)
    x[4]=sqrt(3/7)
    x[1]=-1
    x[5]=1
      
  elseif n == 5
    w[3]=(14+sqrt(7))/30
    w[4]=(14+sqrt(7))/30
    w[2]=(14-sqrt(7))/30
    w[5]=(14-sqrt(7))/30
    w[1]=1/15
    w[6]=1/15
    
    x[3]=-sqrt(1/3-2*sqrt(7)/21)
    x[4]=sqrt(1/3-2*sqrt(7)/21)
    x[2]=-sqrt(1/3+2*sqrt(7)/21)
    x[5]=sqrt(1/3+2*sqrt(7)/21)
    x[1]=-1
    x[6]=1
      
  elseif n == 6
    w[4]=256/525
    w[3]=(124+7*sqrt(15))/350
    w[5]=(124+7*sqrt(15))/350
    w[2]=(124-7*sqrt(15))/350
    w[6]=(124-7*sqrt(15))/350
    w[1]=1/21
    w[7]=1/21
      
    x[4]=0
    x[3]=-sqrt(5/11-2/11*sqrt(5/3))
    x[5]=sqrt(5/11-2/11*sqrt(5/3))
    x[2]=-sqrt(5/11+2/11*sqrt(5/3))
    x[6]=sqrt(5/11+2/11*sqrt(5/3))
    x[1]=-1
    x[7]=1
      
  else
    println("ord1+ord2 zu gross oder zu klein")
  end
  return w,x
end

