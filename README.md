
<p align="center">
  <img src="https://raw.githubusercontent.com/marianoarnaiz/DASvader.jl/main/Documents/Logo.png" alt="DASVader Logo" width="200" />
</p>

# DASVader.jl

**DASVader.jl** is an open, fast, and flexible package for analyzing **Distributed Acoustic Sensing (DAS)** data in [Julia](https://julialang.org). It is designed for ease of use, speed, and adaptability, making it ideal for processing large DAS datasets.

This `README` provides a brief overview of installing and using DASVader. Comprehensive documentation and examples will be available soon.

## Overview

DASVader is a framework designed to read, process, and visualize Distributed Acoustic Sensing (DAS) data, similar to how software like SAC ([SAC - IRIS](https://ds.iris.edu/ds/nodes/dmc/software/downloads/sac/)), PQL ([PQL II](https://epic.earthscope.org/content/pql-ii-program-viewing-data)), and SeisGram ([SeisGram2K](http://alomax.free.fr/seisgram/SeisGram2K.html)) handle more general seismological data.

The framework provides functionality for many typical signal processing steps in both the frequency and wavelength domains. DASVader leverages the excellent Seis.jl ([Seis.jl GitHub](https://github.com/anowacki/Seis.jl)) package for processing, with additional support from other packages like FFTW ([FFTW.jl GitHub](https://github.com/JuliaMath/FFTW.jl)) and FourierAnalysis ([FourierAnalysis.jl GitHub](https://github.com/Marco-Congedo/FourierAnalysis.jl)).

Plotting is a critical component of seismic data analysis, and DASVader enhances this by offering dynamic, interactive visualizations. This is accomplished using a customized version of InteractiveViz.jl ([InteractiveViz.jl GitHub](https://github.com/org-arl/InteractiveViz.jl)), allowing users to explore large datasets without worrying about performance.

Although DASVader is not intended for highly advanced processing techniques such as machine learning-based denoising, we are continuously working on improvements and new features. We welcome contributions, feedback, and feature requests from the community.

## Key Features

- **DAS data processing**: First package in Julia dedicated to DAS processing.
- **Dynamic visualization**: Interactive plots for real-time data exploration.
- **Open-source**: Contributions are welcome!

Feel free to contribute or request new features, and help improve DASVader for the seismological community.


---

## Features
- Process large DAS datasets efficiently.
- Flexible tools for data visualization, transformation, and analysis.
- Designed with Julia’s high-performance capabilities.

---

## Installation

At present, **DASVader** is unregistered, and both it and its dependencies must be installed manually. Follow the steps below to get started:

1. Launch Julia from a terminal or your favorite IDE.
2. Enter **Pkg mode** by pressing `]` in the Julia REPL.
3. Run the following command to add DASVader and its required dependencies:

   ```julia
   (v1.11) pkg> add https://github.com/marianoarnaiz/DASvader.jl https://github.com/anowacki/Geodesics.jl https://github.com/anowacki/Seis.jl
   ```

4. Once the installation is complete, you can start using DASVader by loading it into your Julia session:

   ```julia
   julia> using DASVader
   ```

---
## Usage

### Step 0: Load the Package
```julia
using DASVader
```

### Step 1: Read a raw Febus HDF5 file to memory. 
**Note:** Only FEBUS A1 DAS is readable.
```julia
dDAS = rdas("SR_DS_2023-08-24_14-06-17_UTC.h5")
```

### Step 2: View the data in the file. 
You can change the colormap (`cm`) (e.g., `:grays`, `:viridis`, `:RdBu_9`) and adjust the color limits (`climit`) to your preference.
```julia
fig = viewdas(dDAS; cm=:RdBu_9, climit=10000)
```

### Step 3: Write the figure to a PDF
```julia
savefig(fig, "Matrix.pdf")
```
---
## Some data

### Notice that this data belongs to the FIMOPTIC projec!

If you need some data to test the code you can download this files. The run with the examples provided:

[A Noisy blast file.](https://www.dropbox.com/scl/fi/c6ui9cxcb1gxawm0qqkcp/SR_DS_2023-08-24_14-06-17_UTC_Noisy_Blast.h5?rlkey=e3yn1likn3mkrhesyx85pcpw8&st=6zh6w1i0&dl=0)

[A file with a micro event.](https://www.dropbox.com/scl/fi/xxrd8rlw8kwthmfgwamwx/SR_DS_2023-10-30_12-01-40_UTC_Microevent.h5?rlkey=3zjvn706s46grco4gzhzrqmu5&st=0di6y3x8&dl=0)

[A file with a big blast.](https://www.dropbox.com/scl/fi/abcb1zphevctfkzetn7ql/SR_DS_2024-10-22_14-08-02_UTC_Big_Blast.h5?rlkey=mi9mvj3ynzptgxj3e1khkre8l&st=46eo4s9h&dl=0)

[A file with something that might be an event.](https://www.dropbox.com/scl/fi/n5czzuez7lq2yt3j2s5k5/SR_DS_2024-10-22_21-27-02_UTC_Hidden_Event.h5?rlkey=fsb9tq7wuxgbek5av3uc8rpzd&st=ncrqbwhx&dl=0)

---

## Contribution

Contributions, suggestions, and bug reports are welcome! Please feel free to contact me via email, open an issue or submit a pull request on the [GitHub repository](https://github.com/marianoarnaiz/DASvader.jl).

---

## License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/marianoarnaiz/DASvader.jl/blob/main/LICENSE) file for details.

---

*“Use the force of fast and efficient DAS analysis with DASVader.jl.”*
