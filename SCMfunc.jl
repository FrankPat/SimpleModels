# Function definitions


function CO2dQ(CO2, CO2o)
    return 5.35 * log.(CO2 ./ CO2o)
end

function load_data(filename)
    return readdlm(filename)
end

function savitzkyGolay(x::Vector, windowSize::Int, polyOrder::Int; deriv::Int=0)

  isodd(windowSize) || throw("Window size must be an odd integer.")
  polyOrder < windowSize || throw("Polynomial order must me less than window size.")

  halfWindow = Int( ceil((windowSize-1)/2) )

  # Setup the S matrix of basis vectors
  S = zeros.(windowSize, polyOrder+1)
  for ct = 0:polyOrder
    S[:,ct+1] = (-halfWindow:halfWindow).^(ct)
  end

  ## Compute the filter coefficients for all orders

  # From the scipy code it seems pinv(S) and taking rows should be enough
  G = S * pinv(S' * S)

  # Slice out the derivative order we want
  filterCoeffs = G[:, deriv+1] * factorial(deriv)

  # Pad the signal with the endpoints and convolve with filter
  paddedX = [x[1]*ones(halfWindow); x; x[end]*ones(halfWindow)]
  y = conv(filterCoeffs[end:-1:1], paddedX)

  # Return the valid midsection
  return y[2*halfWindow+1:end-2*halfWindow]

end

function moving_average(data, window_size)
    n = length(data)
    ma = zeros(n)
    for i in 1:n
        if i < window_size
            ma[i] = mean(data[1:i])
        else
            ma[i] = mean(data[i-window_size+1:i])
        end
    end
    return ma
end

