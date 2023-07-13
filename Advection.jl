function fF = Advection(cF,uF,dXdxI,J,DW,S,w)
N = size(cF,2);
n = size(cF,1);
fF = zeros(size(cF));
temp = zeros(n,N);
for i = 1 : N
  temp[:,i] = dXdxI[:,i] .* uF[:,i] .* cF[:,i];
end

for i = 1 : N
  fc = Flux[temp[:,i],S];
  fF[1,i] = fc[1] ;
  for j = 1 : n - 2
    fF[j+1,i] = [fc[j+1] - fc[j]] ;
  end
  fF[n,i] = -fc[n-1] ;
end


for i = 1 : N - 1
  fM = fF[end,i] + fF[1,i+1];
  JM = J[end,i] + J[1,i+1];
  fF[end,i] = fM;
  fF[1,i+1] = fM;
  J[end,i] = JM;
  J[1,i+1] = JM;
end
fM = fF[end,N] + fF[1,1];
JM = J[end,N] + J[1,1];
fF[end,N] = fM;
fF[1,1] = fM;
J[end,N] = JM;
J[1,1] = JM;
for i = 1 : N
  for j = 1 : n
    fF[j,i] = fF[j,i] / [w[j] * J[j,i]];
  end
end
end



