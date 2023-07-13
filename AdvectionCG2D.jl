function AdvectionCG2D!(fF,cF,uF,vF,dXdxI,J,S,w)
  Nx = size(cF,3)
  Ny = size(cF,4)
  n = size(cF,1)
  
  TempH = zeros(n,n) 
  FluxHighX = zeros(n-1,n)
  FluxHighY = zeros(n,n-1)
  uC = zeros(n,n)
  vC = zeros(n,n)
  @. fF = 0.0
  for ix = 1 : Nx
    for iy = 1 : Ny
      @views @. uC = dXdxI[:,:,1,1,ix,iy] * uF[:,:,ix,iy] +
         dXdxI[:,:,1,2,ix,iy] * vF[:,:,ix,iy]
      @views @. vC = dXdxI[:,:,2,1,ix,iy] * uF[:,:,ix,iy] +
         dXdxI[:,:,2,2,ix,iy] * vF[:,:,ix,iy]
      @views @. TempH = uC * cF[:,:,ix,iy]    
      FluxX!(FluxHighX,TempH,S)
      FluxHighX = FluxHighX * diagm(w)
      @views @. fF[1,:,ix,iy] += FluxHighX[1,:] 
      for i = 1 : n - 2
        @views @. fF[i+1,:,ix,iy] += FluxHighX[i+1,:] - FluxHighX[i,:] 
      end
      @views @. fF[n,:,ix,iy] += -FluxHighX[n-1,:]

      @views @. TempH = vC * cF[:,:,ix,iy]    
      FluxY!(FluxHighY,TempH,S)
      FluxHighY = diagm(w) * FluxHighY 
      @views @. fF[:,1,ix,iy] += FluxHighY[:,1] 
      for j = 1 : n - 2
        @views @. fF[:,j+1,ix,iy] += FluxHighY[:,j+1] - FluxHighY[:,j] 
      end
      @views @. fF[:,n,ix,iy] += -FluxHighY[:,n-1]
    end
  end  
  DSS!(fF,J,w)
end
