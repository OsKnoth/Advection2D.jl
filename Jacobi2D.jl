function Jacobi2D(P1,P2,P3,P4,ksi,D)
  
  n = size(D,1)
  X = zeros(n,n,2)
  dXdx=zeros(n,n,2,2)
  J = zeros(n,n)
  dXdxI=zeros(n,n,2,2)
  for i = 1 : n
    for j = 1 : n
      X[i,j,:] = 0.25 *( 
        (1-ksi[i])*(1-ksi[j]) * P1[:] + 
        (1+ksi[i])*(1-ksi[j]) * P2[:] + 
        (1+ksi[i])*(1+ksi[j]) * P3[:] + 
        (1-ksi[i])*(1+ksi[j]) * P4[:])
    end
  end
  for j = 1 : n
    for i = 1 : n
      for l = 1 : n
        dXdx[i,j,:,1] = dXdx[i,j,:,1] + D[i,l] * X[l,j,:]
      end
    end
  end
  for i = 1 : n
    for j = 1 : n
      for l = 1 : n
        dXdx[i,j,:,2] = dXdx[i,j,:,2] + D[j,l] * X[i,l,:]
      end
    end
  end
  for i = 1 : n
    for j = 1 : n
      J[i,j] = det(reshape(dXdx[i,j,:,:],2,2))
      dXdxI[i,j,:,:] = inv(reshape(dXdx[i,j,:,:],2,2)) * J[i,j]
    end
  end
  return X,dXdxI,J
  end
