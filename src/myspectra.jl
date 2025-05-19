"this are the functions for all that stuff on the seismic moment estimation!"

using ScatteredInterpolation, RollingFunctions

# FUNCTION TO GET THE FREQUENCY SPECTRUM WITH ZERO-PADDING!
function get_spectrum(signal, fs)
    # Get size and time sampling interval
    n_samples = size(signal, 1)
    Δt = 1 / fs  # Time sampling interval

    # Determine the next power of 2 for zero-padding
    padded_length = nextpow(2, n_samples)

    # Zero-pad the signal to the next power of 2
    padded_signal = vcat(signal, zeros(padded_length - n_samples))

    # Compute the FFT
    itrace = CartTrace(0.0, Δt, padded_signal)
    taper!(itrace)
    padded_signal = trace(itrace)
    fft_result = fft(padded_signal)  # Perform FFT

    # Compute frequency bins for the padded signal
    freq_bins = (0:padded_length-1) * (fs / padded_length)  # Frequency vector

    # Get the positive frequencies and corresponding magnitudes
    half_n = div(padded_length, 2) + 1
    freq_bins = freq_bins[1:half_n]
    magnitudes = abs.(fft_result[1:half_n]) .* Δt  # Normalize by time sampling interval Δt

    # Avoid 0 Hz component
    fv = freq_bins[2:end]       # Frequency vector (excluding 0 Hz)
    smv = magnitudes[2:end]      # Spectrum magnitude vector (excluding 0 Hz)

    return fv, smv  # Return frequency vector and spectrum magnitude vector
end

# FUCNTION TO SMOOTH THE FREQUENCY SPECTRUM!

function smooth_spectrum(fv, smv)
    x = fv
    y = smv
    y2 = [y[1:2]; rollmean(y, 5); y[end-1:end]]
    #y2=[y[1]; rollmean(y, 3); y[end]]
    y3 = [y[1]; rollmean(y2, 3); y[end]]
    newx = 10 .^ range(log10(x[1]), stop=log10(x[end]), length=20 * ceil(Int, log10(x[end])))
    itp = interpolate(Polyharmonic(1), x', y3)
    newy = evaluate(itp, newx')
    sfv = newx
    ssmv = newy
    return sfv, ssmv #smooth frequency vector(x), smooth spectrum magnitude vector(y)
end

# FUCNTION TO INTEGRATE THE FREQUENCY SPECTRUM TWICE!
function divomega2(fv, smv)
    f = fv
    y = smv
    ω = 2 * pi * f
    y2 = y ./ (ω .^ 2)
    return y2
end

#FUNCTION FOR THE SPECTRUM
function displacement_spectrum(f, M0, fc, Q, r, c, ρ, Rθϕ=0.62, F=1, G=1)
    # Parameters:
    # f  - frequency (Hz)
    # M0 - seismic moment (N·m)
    # fc - corner frequency (Hz)
    # Q  - quality factor (dimensionless)
    # r  - distance to source (m)
    # c  -  wave velocity (m/s)
    # ρ  - rock density (kg/m^3)
    # Rθφ - radiation pattern coefficient (default: 0.62)
    # F  - free surface coefficient (default: 1)
    # G  - geometrical spreading (default: 1)
    A = 1 / r^G
    B = F * Rθϕ / (4 * π * ρ * c^3)
    C = M0
    D = 1 ./ (1 .+ (f ./ fc) .^ 2)
    E = exp.(-π .* f .* r ./ (Q * c))
    # Calculate the spectral amplitude A(f)
    A_f = A .* B .* C .* D .* E

    return A_f
end


export get_spectrum, smooth_spectrum, divomega2
