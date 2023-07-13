function Plot2DC(c,IntXF2E,Name)
  Nx = size(c,3)
  Ny = size(c,4)
  OrdPoly = size(c,1) - 1
  cPlot = zeros(Nx*OrdPoly,Ny*OrdPoly)

  for ix = 1 : Nx
    for iy = 1 : Ny  
      cPlot[(ix-1)*OrdPoly + 1 : ix*OrdPoly, (iy-1)*OrdPoly + 1 : iy*OrdPoly] =
#     cPlot[(iy-1)*OrdPoly + 1 : iy*OrdPoly, (ix-1)*OrdPoly + 1 : ix*OrdPoly] =
        (IntXF2E*reshape(c[:,:,ix,iy],OrdPoly+1,OrdPoly+1)*IntXF2E')
    end
  end  
  fig = Figure()
  ax = Axis(fig[1, 1])
  hm = heatmap!(ax, cPlot)
  Colorbar(fig[1, 2], hm)
  fig
  save(Name*".png", fig)
end  

