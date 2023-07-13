function AdvectionFV2D!(fF,cF,uF,vF,dXdxI,J,w)
  Nx = size(cF,3)
  Ny = size(cF,4)
  n = size(cF,1)
  
  
  @. fF = 0.0
  uC = zeros(n,n)
  vC = zeros(n,n)
  for ix = 1 : Nx
    for iy = 1 : Ny
      @views @. uC = dXdxI[:,:,1,1,ix,iy] * uF[:,:,ix,iy] + 
           dXdxI[:,:,1,2,ix,iy] * vF[:,:,ix,iy]
      @views @. vC = dXdxI[:,:,2,1,ix,iy] * uF[:,:,ix,iy] + 
           dXdxI[:,:,2,2,ix,iy] * vF[:,:,ix,iy]
      for j = 1 : n 
        for i = 1 : n - 1
          uCL = 0.5 * (uC[i,j] + uC[i+1,j])
          if uCL > 0.0
            fluxLoc = uCL * cF[i,j,ix,iy] * w[j]
          else
            fluxLoc = uCL * cF[i+1,j,ix,iy] * w[j]
          end
          fF[i,j,ix,iy] = fF[i,j,ix,iy] - fluxLoc
          fF[i+1,j,ix,iy] = fF[i+1,j,ix,iy] + fluxLoc
        end
      end
      for i = 1 : n 
        for j = 1 : n - 1
          vCL = 0.5 * (vC[i,j] + vC[i,j+1])
          if vCL > 0.0
            fluxLoc = vCL * cF[i,j,ix,iy] * w[i]
          else
            fluxLoc = vCL * cF[i,j+1,ix,iy] * w[i]
          end
          fF[i,j,ix,iy] = fF[i,j,ix,iy] - fluxLoc
          fF[i,j+1,ix,iy] = fF[i,j+1,ix,iy] + fluxLoc
        end
      end
    end
  end
  DSS!(fF,J,w)
end
  
