function DSS!(fF,J,w)  
  n = size(fF,1)
  Nx = size(fF,3)
  Ny = size(fF,4)
  JJ = deepcopy(J)
  fM = zeros(n)
  JM = zeros(n)
  for iy = 1 : Ny
    for ix = 1 : Nx - 1
      @views @. fM = fF[end,:,ix,iy] + fF[1,:,ix+1,iy]
      @views @. JM = JJ[end,:,ix,iy] + JJ[1,:,ix+1,iy]
      @views @. fF[end,:,ix,iy] = fM
      @views @. fF[1,:,ix+1,iy] = fM
      @views @. JJ[end,:,ix,iy] = JM
      @views @. JJ[1,:,ix+1,iy] = JM
    end
    @views @. fM = fF[end,:,Nx,iy] + fF[1,:,1,iy]
    @views @. JM = JJ[end,:,Nx,iy] + JJ[1,:,1,iy]
    @views @. fF[end,:,Nx,iy] = fM
    @views @. fF[1,:,1,iy] = fM
    @views @. JJ[end,:,Nx,iy] = JM
    @views @. JJ[1,:,1,iy] = JM
  end
  for ix = 1 : Nx
    for iy = 1 : Ny - 1
      @views @. fM = fF[:,end,ix,iy] + fF[:,1,ix,iy+1]
      @views @. JM = JJ[:,end,ix,iy] + JJ[:,1,ix,iy+1]
      @views @. fF[:,end,ix,iy] = fM
      @views @. fF[:,1,ix,iy+1] = fM
      @views @. JJ[:,end,ix,iy] = JM
      @views @. JJ[:,1,ix,iy+1] = JM
    end
    @views @. fM = fF[:,end,ix,Ny] + fF[:,1,ix,1]
    @views @. JM = JJ[:,end,ix,Ny] + JJ[:,1,ix,1]
    @views @. fF[:,end,ix,Ny] = fM
    @views @. fF[:,1,ix,1] = fM
    @views @. JJ[:,end,ix,Ny] = JM
    @views @. JJ[:,1,ix,1] = JM
  end
  for ix = 1 : Nx
    for iy = 1 : Ny
      for j = 1 : n
        for i = 1 : n
          fF[i,j,ix,iy] = fF[i,j,ix,iy] / (w[i] * w[j] * JJ[i,j,ix,iy])
        end
      end
    end
  end
end
  
