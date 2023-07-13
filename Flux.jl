function FluxX!(fc,c,S)
  n = size(c,1)
  @inbounds for i = 1 : n -1
    @views @. fc[i,:] = 0.0
    @inbounds for l = 1 : i
      @inbounds for k = 1 : n
        @views @. fc[i,:] = fc[i,:] - 0.5 * S[l,k] * (c[l,:] + c[k,:])
      end
    end
  end
end

function FluxY!(fc,c,S)
  n = size(c,1)
  @inbounds for i = 1 : n -1
    @views @. fc[:,i] = 0.0
    @inbounds for l = 1 : i
      @inbounds for k = 1 : n
        @views @. fc[:,i] = fc[:,i] - 0.5 * S[l,k] * (c[:,l] + c[:,k])
      end
    end
  end
end
