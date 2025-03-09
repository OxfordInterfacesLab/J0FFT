from PIL import Image
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.fft import fft, ifft, fftfreq
from scipy.signal import find_peaks, peak_prominences
import h5py
from scipy.interpolate import interp1d


def pl_line_contrast(pl_line, plot_option,image_width_mm, center_frequency):
    """
    Process a single line of the PL image, apply FFT filtering, and plot results.

    Parameters:
        pl_line (numpy array): One line of the PL image.

    Returns:
        contrast_percentage (float): FFT contrast in percent.
    """
    # Input line from the image
    original_line = pl_line

    # Line length
    line_length = len(original_line)

    # Copy of the original line for processing
    processed_line = original_line.copy()

    # Define image length parameters for calibration to mm
    
    spatial_resolution = image_width_mm / line_length
    spatial_axis = np.arange(0, image_width_mm, spatial_resolution)  # x vector for plotting

    # Frequency vector for k-space
    sampling_rate = line_length / image_width_mm  # Sampling frequency [samples/mm]
    frequency_axis = fftfreq(len(processed_line), 1 / sampling_rate)

    # Normalize using the mean value
    normalized_line = processed_line / np.mean(processed_line)  # Normalized PL to its mean

    # Perform FFT
    fft_result = fft(normalized_line)
    fft_magnitude = np.abs(fft_result)  # FFT integrals are complex; this is the modulus

    # Frequency filtering parameters
    delta_bandwidth = 0.02  # Half bandwidth
    
    band_start, band_end = center_frequency - delta_bandwidth, center_frequency + delta_bandwidth
    normalization_bandwidth = 0.02  # Normalization bandwidth

    # Create frequency-selective filter (band-pass)
    bandpass_filter = (np.abs(frequency_axis) >= band_start) & (np.abs(frequency_axis) <= band_end)
    filtered_fft = fft_result * bandpass_filter

    # Create normalization filter to remove noise
    lower_band_filter = (np.abs(frequency_axis) >= band_start - normalization_bandwidth) & (np.abs(frequency_axis) <= band_start)
    upper_band_filter = (np.abs(frequency_axis) >= band_end) & (np.abs(frequency_axis) <= band_end + normalization_bandwidth)
    normalization_filter = lower_band_filter | upper_band_filter

    # Remove the base level from the FFT
    idx_norm = np.nonzero(normalization_filter)[0]
    idx_filter = np.nonzero(bandpass_filter)[0]
    adjusted_fft = fft_result - np.mean(np.real(fft_result[idx_norm]))
    adjusted_fft[idx_filter[:len(idx_filter) // 2]] -= 1j * np.mean(np.imag(fft_result[idx_norm[:len(idx_norm) // 2]]))
    adjusted_fft[idx_filter[len(idx_filter) // 2:]] -= 1j * np.mean(np.imag(fft_result[idx_norm[len(idx_norm) // 2:]]))

    # Apply the filter and perform inverse FFT
    enhanced_fft_sim_PL_signal = adjusted_fft * bandpass_filter
    enhanced_sim_PL_signal = ifft(enhanced_fft_sim_PL_signal, n=len(processed_line))

    if plot_option == 1:
        
        # Plot the original and enhanced sim_PL_signals
        fig, axs = plt.subplots(2, 2, figsize=(12, 6))

        # Subplot 1: Original line
        axs[0, 0].plot(spatial_axis, normalized_line, 'royalblue')
        axs[0, 0].set_xlabel('Distance (mm)', fontsize=20)
        axs[0, 0].set_ylabel('Normalized PL (a.u.)', fontsize=20)
        axs[0, 0].set_xlim(0, image_width_mm)
        axs[0, 0].set_ylim(0.75, 1.1)
        axs[0, 0].tick_params(axis='both', which='major', labelsize=18)

        # Subplot 2: Full k-spectrum
        half_freq = frequency_axis[:len(frequency_axis) // 2]
        half_fft_magnitude = fft_magnitude[:len(fft_magnitude) // 2]
        axs[0, 1].plot(half_freq, np.log(half_fft_magnitude), 'royalblue')
        axs[0, 1].set_xlim(0, 1 / 3 * max(frequency_axis))
        axs[0, 1].set_ylim(-4, 8)
        axs[0, 1].set_xlabel('Spatial frequency k (1/mm)', fontsize=20)
        axs[0, 1].set_ylabel('FFT of PL profile (logscale)', fontsize=20)
        axs[0, 1].tick_params(axis='both', which='major', labelsize=18)

        # Subplot 3: Zoomed k-spectrum with filters
        ax1 = axs[1, 0]
        ax2 = ax1.twinx()
        ax1.plot(half_freq, np.log(half_fft_magnitude), 'royalblue')
        ax1.set_xlim(0.4, 0.92)
        ax1.set_ylim(-2.0, 1.5)
        ax1.set_xlabel('Spatial frequency k (1/mm)', fontsize=20)
        ax1.set_ylabel('Log.abs.fft(PL counts)', color='royalblue', fontsize=20)
        ax1.tick_params(axis='y', labelcolor='royalblue', labelsize=18)

        # Plot filters
        def filter_bandpass(freq):
            return 1 if band_start <= freq <= band_end else 0

        def filter_normalization(freq):
            return 1 if (band_start - normalization_bandwidth <= freq <= band_start) or \
                        (band_end <= freq <= band_end + normalization_bandwidth) else 0

        bandpass_plot = np.vectorize(filter_bandpass)(half_freq)
        normalization_plot = np.vectorize(filter_normalization)(half_freq)
        ax2.plot(half_freq, normalization_plot, color='gold', label='Normalization Filter')
        ax2.plot(half_freq, bandpass_plot, color='forestgreen', label='Bandpass Filter')
        ax2.set_ylim(-0.5, 1.5)
        ax2.set_ylabel('Filter Value', color='forestgreen', fontsize=20)
        ax2.tick_params(axis='y', labelcolor='forestgreen', labelsize=18)

        # Subplot 4: Reconstructed sim_PL_signal
        axs[1, 1].plot(spatial_axis, np.real(enhanced_sim_PL_signal), 'royalblue')
        axs[1, 1].set_ylim(-0.02, 0.02)
        axs[1, 1].set_xlim(0, image_width_mm)
        axs[1, 1].set_xlabel('Distance (mm)', fontsize=20)
        axs[1, 1].set_ylabel('Normalized Amplitude', fontsize=20)
        axs[1, 1].tick_params(axis='both', which='major', labelsize=18)

        plt.tight_layout()
        plt.show()

    # Analyze peaks and calculate contrast
    peaks, _ = find_peaks(np.real(enhanced_sim_PL_signal))
    peak_prominence_values = peak_prominences(np.real(enhanced_sim_PL_signal), peaks)[0]
    contrast_percentage = np.mean(peak_prominence_values) / np.mean(normalized_line) * 100  # Contrast in percentage

    return contrast_percentage

    
## For running Quokka we need the function below, that extracts the data.

def get_Q3_PL_signal(Q3_resultsfile, Nx):
    # gets luminescence profile along X-axis from Quokka3 resultsfile
    # returns list of dictionary with keys 'X' [cm] and 'Y' [a.u.] being an equidistant luminescence profile
    # if a sweep over J0cont is defined the list has an entry for each J0cont, otherwise it is a single-entry list
    
    def read_interpolated_lumi_profile(h5file,sweepName,Nx):
        lumi_profile = dict()
        X = np.array(h5file['Spatial' + sweepName + '/X'][:])
        if 'LumiIntensity(X)' in h5file['Curves' + sweepName]:
            # full version, directly read luminescence profile which includes detailed reabsorption modeling
            I = h5file['Curves' + sweepName + '/LumiIntensity(X)']
        else:
            # free version, calculate luminescence from 2D carrier profile with simple reasborption model
            alpha = 0.76625 # absorption coefficient for single-wavelength reabsorption model as in paper / PC3D
            Bradni2 = 4.853e5 # Brad * ni^2 from Nguyen at 300K (does not influence luminescence contrast)
            Vt = 0.02586 # thermal voltage at 300 K
            Z = np.array(h5file['Spatial' + sweepName + '/Z'][:])
            QFPsplit = np.array(h5file['Spatial' + sweepName + '/QFPsplit'][:]) # 2D profile of quasi Fermi-level split
            I = np.zeros(X.shape)
            for ix in np.arange(0,len(X)):
                for iz in np.arange(0,len(Z)):
                    I[ix] += Bradni2 * np.exp(QFPsplit[0][ix][iz]/Vt - 1) * np.exp(-alpha*Z[iz])
                    
        # now interpolate to equidistant grid
        f = interp1d(X, I, kind='cubic', fill_value="extrapolate")
        lumi_profile['X'] = np.linspace(0, h5file['Spatial' + sweepName].attrs['Wx'], Nx)
        lumi_profile['Y'] = f(lumi_profile['X'])
        return lumi_profile
    
    lumi_profiles = []
    try:
        with h5py.File(Q3_resultsfile, 'r') as h5file:
            if h5file.attrs['isSweep'] == 1:
                J0cont = np.array(h5file['Sweep']['GroupA(1)']['Values'][:])
                for group in h5file['Curves']:
                    lumi_profiles.append(read_interpolated_lumi_profile(h5file, '/' + group, Nx))
            else:
                J0cont = []
                lumi_profiles.append(read_interpolated_lumi_profile(h5file, '', Nx))
        return lumi_profiles, J0cont
                
    except Exception as e:
        print(f"An error occurred getting luminescence profile from Quokka3 resultsfile: {e}")
	

