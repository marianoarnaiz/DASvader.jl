using FFTW, GLMakie

# Define a structure to store spectrogram components
struct Spec
    amp::Matrix{Float64}     # Amplitude matrix
    time::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64}   # Time vector
    freq::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int64}    # Frequency vector
end

# Gaussian window function
function gaussian_window(n, alpha=2.5)
    # alpha controls the width of the Gaussian window
    k = -(n - 1)/2:(n-1)/2
    return exp.(-(alpha * k / ((n - 1) / 2)) .^ 2)
end


# Function to compute a spectrogram with Gaussian window
function gspectrogram(signal, fs; window_size=2^ceil(Int, log2(size(signal, 1) / 8)), overlap=0.5, alpha=2.5)
    stepw = Int(round(window_size * (1 - overlap)))                         # stepw size between windows
    n_windows = div(length(signal) - window_size, stepw)             # Number of windows
    half_size = div(window_size, 2) + 1                             # Half-size for FFT (Nyquist)

    # Preallocate spectrogram matrix for only the positive frequencies
    spectrogram = zeros(half_size, n_windows)
    window = gaussian_window(window_size, alpha)

    for i in 1:n_windows
        start_idx = (i - 1) * stepw + 1
        segment = signal[start_idx:start_idx+window_size-1] .* window
        spectrum = abs2.(fft(segment)[1:half_size])                  # Take only positive half
        spectrogram[:, i] = spectrum
    end

    time_vector = range(start=0.0, step=stepw / fs, length=n_windows)# (0:n_windows - 1) * (stepw / fs)
    freq_vector = range(start=0.0, step=fs / window_size, length=half_size) #(0:half_size - 1) * (fs / window_size)

    return Spec(spectrogram, time_vector, freq_vector)
end

# Plotting the spectrogram
function plot_spectrogram(spec::Spec; plot_type=:contour)
    fig = Figure(; size=(800, 600))
    ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Frequency (Hz)", title="Spectrogram")

    if plot_type == :heatmap
        heatmap!(ax, spec.time, spec.freq, spec.amp', colormap=:plasma)
    elseif plot_type == :contour
        contourf!(ax, spec.time, spec.freq, spec.amp', levels=25, colormap=:plasma)
    end

    fig
end

# Compute and plot the spectrogram with Gaussian window
#spec = gspectrogram(ch, 10000, window_size=512, overlap=0.5, alpha=2.5)
#plot_spectrogram(spec)
export gaussian_window, plot_spectrogram, gspectrogram
