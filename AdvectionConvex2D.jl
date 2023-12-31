function AdvectionConvex2D!(fF,cF,uF,vF,dXdxI,J,S,w)
  Nx = size(cF,3)
  Ny = size(cF,4)
  n = size(cF,1)
  
  TempH = zeros(n,n) 
  FluxHighX = zeros(n-1,n)
  FluxLowX = zeros(n-1,n)
  FluxX = zeros(n-1,n)
  uCL = zeros(n-1,n)
  alphaX = zeros(n-1,n)
  FluxHighY = zeros(n,n-1)
  FluxLowY = zeros(n,n-1)
  FluxY = zeros(n,n-1)
  vCL = zeros(n,n-1)
  alphaY = zeros(n,n-1)
  uC = zeros(n,n)
  vC = zeros(n,n)
  dxF = zeros(n)
  dxC = zeros(n-1)
  dyF = zeros(n)
  dyC = zeros(n-1)
  e = zeros(n-2)
  @. fF = 0.0
  @inbounds for ix = 1 : Nx
    @inbounds for iy = 1 : Ny
      @views @. uC = dXdxI[:,:,1,1,ix,iy] * uF[:,:,ix,iy] +
         dXdxI[:,:,1,2,ix,iy] * vF[:,:,ix,iy]
      @views @. vC = dXdxI[:,:,2,1,ix,iy] * uF[:,:,ix,iy] +
         dXdxI[:,:,2,2,ix,iy] * vF[:,:,ix,iy]

      @views @. TempH = uC * cF[:,:,ix,iy]    
      FluxX!(FluxHighX,TempH,S)
      FluxHighX = FluxHighX * diagm(w)

      @views @. uCL = 0.5 * (uC[1:end-1,:] + uC[2:end,:])
      @views @. FluxLowX = -0.5 *  ((abs(uCL) + uCL) * cF[1:end-1,:,ix,iy] +
            (uCL - abs(uCL)) * cF[2:end,:,ix,iy])
      FluxLowX = FluxLowX * diagm(w)

      @inbounds for j = 1 : n
        @views @. dxF = w*J[:,j,ix,iy]
        dxF[1] = 2 * dxF[1]
        dxF[end] = 2 * dxF[end]
        dxC = 0.5*(dxF[1:end-1] + dxF[2:end])
        @views Error!(e,cF[:,j,ix,iy,],dxC)
        @. e = min(e,1.0)
        alphaX[1,j] = e[1]
        @. alphaX[2:end-1,j] = max(e[1:end-1],e[2:end])
        alphaX[end,j] = e[end]
      end
      @. FluxX = alphaX * FluxLowX + (1.0 - alphaX) * FluxHighX
      @views @. fF[1,:,ix,iy] += FluxX[1,:] 
      @views @. fF[2:end-1,:,ix,iy] += FluxX[2:end,:] - FluxX[1:end-1,:] 
      @views @. fF[end,:,ix,iy] += -FluxX[end,:]

      @views @. TempH = vC * cF[:,:,ix,iy]    
      FluxY!(FluxHighY,TempH,S)
      FluxHighY = diagm(w) * FluxHighY 

      @views @. vCL = 0.5 * (vC[:,1:end-1] + vC[:,2:end])
      @views @. FluxLowY = -0.5 *  ((abs(vCL) + vCL) * cF[:,1:end-1,ix,iy] +
            (vCL - abs(vCL)) * cF[:,2:end,ix,iy])
      FluxLowY = diagm(w) * FluxLowY 
      @inbounds for i = 1 : n
        @views @. dyF = w*J[i,:,ix,iy]
        dyF[1] = 2 * dyF[1]
        dyF[end] = 2 * dyF[end]
        dyC = 0.5*(dyF[1:end-1] + dyF[2:end])
        @views Error!(e,cF[i,:,ix,iy,],dyC)
        @. e = min(e,1.0)
        alphaY[i,1] = e[1]
        @. alphaY[i,2:end-1] = max(e[1:end-1],e[2:end])
        alphaY[i,end] = e[end]
      end

      @. FluxY = alphaY * FluxLowY + (1.0 - alphaY) * FluxHighY
      @views @. fF[:,1,ix,iy] += FluxY[:,1] 
      @inbounds for j = 1 : n - 2
        @views @. fF[:,j+1,ix,iy] += FluxY[:,j+1] - FluxY[:,j] 
      end
      @views @. fF[:,n,ix,iy] += -FluxY[:,n-1]
    end
  end  
  DSS!(fF,J,w)
end
