# J0FFT - open and run 'J0FFT_vX.ipynb' as a Jupyter Notebook
A python based code to extract the dark recombination current density at metal contacts (J0c) in half metallised solar cells.

# Extracting Contact Recombination from FFT-filtered Photoluminescence Imaging

## Photoluminescence Imaging (PLI)

Photoluminescence imaging (PLI) is a powerful tool for evaluating carrier recombination processes in photovoltaic materials. In this study, PLI was used to measure luminescence intensity patterns, which correlate with recombination at the metal-silicon interface. Samples with a tunneling oxide passivated contact (TOPCon) structure and industrial screen-printed metallisation were used.

## Fourier Transform Filtering

The primary analytical method involves converting spatial photoluminescence data into the frequency domain using the Fast Fourier Transform (FFT).

### Key Steps:
1. **Transform to Frequency Domain**:
   - PL line profiles are converted into the frequency domain using FFT:
    
     $X_k = \sum_{n=0}^{N-1} x_n e^{-i 2\pi k n / N}$
     
     where:
     - $X_k$: Frequency domain representation.
     - $x_n$: Spatial domain data.
     - $N$: Total number of data points.

2. **Bandpass Filtering**:
   - A selective filter isolates the fundamental spatial frequency of the metallisation pattern.
   - Noise and higher harmonics are suppressed by applying a rectangular window:
  
     $\text{Filtered sim\\_PL\\_signal} = \text{FFT}^{-1}(X_k \cdot W_k)$
  
     where $W_k$ is the filter window.

3. **Reconstruction in Spatial Domain**:
   - The inverse FFT reconstructs the filtered sim_PL_signal back into the spatial domain:

     $x_n = \frac{1}{N} \sum_{k=0}^{N-1} X_k e^{i 2\pi k n / N}$

4. **Metal Contrast Calculation**:
   - Metal contrast $C_{\text{metal}}$ is defined as the average amplitude of oscillations:

     $C_{\text{metal}} = \frac{A_{\text{ave}}}{S_{\text{ave}}} \times 100\%$

     where:
     - $A_{\text{ave}}$: Average oscillation amplitude.
     - $S_{\text{ave}}$: Average sim_PL_signal level.




## Numerical Simulations
Numerical device simulations link the observed metal contrast to the metal recombination parameter $J_0,c$.

### Simulation Tools:
1. **PC3D**:
   - Solves 2D carrier transport equations.
   - Models carrier density distributions and luminescence profiles.

2. **Quokka3**:
   - Simulates spatial PL emission patterns with optical filtering effects.
   - Incorporates material-specific properties for higher accuracy.

### Key Simulation Outputs:
- **Carrier Profiles**:
  - 2D maps of excess carrier densities as a cuntion of depth and lateral position in the cell, under steady-state illumination.
- **PL Profiles**:
  - Simulated luminescence intensity as a function of lateral position, for different metal contact recombination $J_0,c$.

### Relationship Between $J_{0,c}$ and Metal Contrast:
A linear correlation was established between simulated PL metal contrast and $J_0,c$ values:
$C_{\text{metal}} \propto (J_{0,c} - J_{0,\text{TOPCon}})$
where $J_0,\text{TOPCon}$ is the baseline recombination at non-contacted surfaces.


