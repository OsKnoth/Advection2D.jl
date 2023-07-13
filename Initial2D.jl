function Initial2D(x,y,Case)
  #f = exp(-(x-3)^2)
  if Case == "SmoothExp"
    f = exp(-0.1*((x-5)^2+(y-5)^2))
  elseif Case == "Cube"  
    if x >= 400 && x <= 600 && y >= 400 && y <= 600
      f = 1
    else
      f = 0
    end  
  else
    f = 1.0  
  end  
  return f
end
