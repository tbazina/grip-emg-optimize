#!/usr/bin/python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numba as nb
import numpy as np
import pandas as pd
import scipy.fft
import scipy.signal


@nb.jit(nb.float64[:, :](nb.float64[:], nb.int64), nopython=True)
def rolling_window_to_array(emg_arr, window_size):
    """Apply rolling window to 1D numpy array and return with shape
    [length - window_size + 1, window_size]

    Args:
        emg_arr (np.array): emg signal to apply rolling window on
        window (int): size of the window

    Returns:
        2D np.array: array with window_size in each row
    """
    shape = emg_arr.shape[:-1] + (emg_arr.shape[-1] - window_size + 1, window_size)
    strides = emg_arr.strides + (emg_arr.strides[-1],)
    return np.lib.stride_tricks.as_strided(emg_arr, shape=shape, strides=strides)


class SignalProcessingParams:
    """Problem Class"""

    def __init__(
        self,
        measure_position,
        in_data_file="measurements_june_2024/csv/emg_grip_dat_complete.zip",
    ):
        self.in_data_file = in_data_file
        self.measure_position = measure_position
        self.file_names = []
        self.time_dat = []
        self.emg_fft = []
        self.grip_dat = []
        self.emg_dat = []
        self.sampling_rate = None
        self.window_size = None
        self.emg_fft_freq = None
        self.fft_dim = None

        self.data_loader()

    def data_loader(self):
        file_names = []
        time_dat = []
        emg_dat = []
        grip_dat = []
        window_size_relative = 0.5
        # Smoothing dimension (decay factor for exponential smoothing + window size)
        self.smooth_dim = 2

        # Load data from aligned and upsampled compressed csv
        filename = self.in_data_file
        data = pd.read_csv(filename)
        # Convert timestamp to float seconds from epoch
        data.measure_ts = data.measure_ts.astype("datetime64[ns]").astype(int) / 1e9
        # Filter data based on selected measurement position
        data = data[data.position == self.measure_position]

        # Group data by initials, position and repetition
        data_grouped = data.groupby(by=["initials", "position", "repetition"])
        # Loop over grouped df and store data as a list of arrays
        for _, group_df in data_grouped:
            # Timestamps
            time_dat.append(group_df.measure_ts.to_numpy())
            # EMG data
            emg_dat.append(group_df.emg.to_numpy())
            # Grip force data
            grip_dat.append(group_df.grip_force.to_numpy())
            # File names
            file_names.append(np.unique(group_df.filename.to_numpy()))

        for file_name, time_t, emg, grip in zip(
            file_names, time_dat, emg_dat, grip_dat
        ):
            # Data length
            dat_len = time_t.shape[0]
            # print(f'{dat_len=}')
            # Calculate sampling rate from time (remove first and last 10 data points)
            # Using median instead of mean - measurement imperfections
            self.sampling_rate = 1.0 / np.median(np.diff(time_t[10:-10]))
            # print(f'{self.sampling_rate=}')
            # Calculate window size for FFT - always round to nearest even number
            self.window_size = round(window_size_relative * self.sampling_rate)
            self.window_size += self.window_size % 2
            # print(f'{self.window_size=}')
            # FFT resolution
            self.fft_resolution = self.sampling_rate / self.window_size
            # print(f'{self.fft_resolution=}')

            # Reshape data to window_size × n array (disregard starting few points)
            time_t = time_t[dat_len % self.window_size :].reshape(
                (-1, self.window_size)
            )
            emg = emg[dat_len % self.window_size :].reshape((-1, self.window_size))
            grip = grip[dat_len % self.window_size :].reshape((-1, self.window_size))

            # Store filenames
            self.file_names.append(file_name)
            # Store emg data
            self.emg_dat.append(emg)
            # Calculate RFFT on emg_data
            self.emg_fft.append(scipy.fft.rfft(emg))
            # Store grip data
            self.grip_dat.append(grip)
            # Store time data
            self.time_dat.append(time_t)
            # Calculate RFFT frequencies and stack them
            self.emg_fft_freq = scipy.fft.rfftfreq(
                n=self.window_size, d=1.0 / self.sampling_rate
            )
            # print(f'{self.emg_fft_freq=}')
            # FFT problem dimension
            self.fft_dim = self.emg_fft_freq.shape[0]
            # Problem dimension
            self.dim = self.fft_dim + self.smooth_dim

            # Initialize mask with ones - no change in amplitudes
            self.mask = np.ones(self.fft_dim)
            # Always remove first frequency - DC component
            self.mask[0] = 0
            # Optimal indices of frequencies to optimize in decision vector
            # self.optim_ind = [0, 1, 24, 59, 61, 67, 87, 206, 235]
            # Shift indices by 1 to account for not including DC component in optimization
            # self.optim_ind = [i + 1 for i in self.optim_ind]
            # TODO: select all indices (except DC offset) for sensitivity analysis
            # and truncate later
            self.optim_ind = [i + 1 for i in range(self.fft_dim - 1)]

    def get_bounds(self):
        # TODO: Calculate rolling window size for smoothing (must be <= FFT window size)
        # Lower bounds for decision vector
        # lower_bounds = (
        #     [0.0] * 2 + [0.0] + [0.0] * 2 + [0.0] * 2 + [0.0] * 2 + [0.0] + [300]
        # )
        # Upper bounds for decision vector
        # upper_bounds = (
        #     [0.0] * 2
        #     + [1.0]
        #     + [3.0] * 2
        #     + [5.0] * 2
        #     + [0.0] * 2
        #     + [0.004]
        #     + [self.window_size - 1]
        upper_bounds = [5.0] * (self.fft_dim - 1) + [0.05] + [self.window_size - 1]
        # Lower bounds set to zero for all frequencies, 0 for decay factor and 2 for window size
        lower_bounds = [0.0] * (self.fft_dim - 1) + [0.0] + [2]
        # Upper bounds set to five for all frequencies, 0.05 for decay factor and
        # self.window_size - 1 for smoothing window size
        upper_bounds = [5.0] * (self.fft_dim - 1) + [0.05] + [self.window_size - 1]
        return (lower_bounds, upper_bounds)

    def get_nix(self) -> int:
        return 1

    def get_name(self):
        return "Optimize EMG signal processing parameters to cross-corelate with grip force"

    def fitness(self, x):
        """Fitness function

        Args:
          x (np.float64): Decision vector
            x[0:len(self.optim_ind)] - Custom maks values for fft spectrum of signal (0-5)
            x[-2] - Exponential moving average decay (0-0.05)
            x[-1] - Rolling window size for smoothing in sec (1-495)

            - amplification - maybe only in frequency domain 0 - 2 instead 0-1
        """
        max_corrs = []
        # Set decision vector values only to index of frequencies necessary to optimize
        # determined by self.optim_ind
        self.mask[self.optim_ind] = [x[: len(self.optim_ind)]]

        for emg_fft, grip_dat in zip(self.emg_fft, self.grip_dat):
            # Compute inverse FFT to reconstruct filtered signal as n-dim array
            emg_ifft = scipy.fft.irfft(emg_fft * self.mask[: self.fft_dim])
            rolling_window_size = round(x[-1])
            # Rectify flattened iFFT, construct window for smoothing and
            # calculate exponential moving average using FFT convolution
            emg_abs = np.abs(emg_ifft.ravel())
            window_ema = np.array(
                [(1 - x[-2]) ** i for i in range(rolling_window_size)]
            )
            window_ema = window_ema / window_ema.sum()
            emg_ema = scipy.signal.fftconvolve(emg_abs, window_ema[::-1], mode="valid")
            grip_flat = grip_dat.ravel()[rolling_window_size - 1 :]
            emg_ema = emg_ema - emg_ema.mean()
            grip_flat = grip_flat - grip_flat.mean()

            corr_ema = scipy.signal.fftconvolve(emg_ema, grip_flat[::-1], mode="full")
            corr_ema /= len(grip_flat) * emg_ema.std() * grip_flat.std()

            max_corr = np.abs(corr_ema).max()
            max_corrs.append(max_corr)

        avg_corr = np.mean(max_corrs)

        return [1 - avg_corr]

    def get_extra_info(self):
        return (
            f"\tSampling rate: {self.sampling_rate} Hz\n"
            f"\tFFT window size: {self.window_size} samples\n"
            f"\tFFT resolution: {self.fft_resolution} Hz\n"
        )

    def plot(self, x):
        fig, axs = plt.subplots(10, 6, figsize=(16, 30), dpi=320)
        self.mask[self.optim_ind] = [x[:9]]
        i = 0
        row = 0
        col = 0

        for file_name, emg_fft, grip_dat, time_dat in zip(
            self.file_names, self.emg_fft, self.grip_dat, self.time_dat
        ):
            # Compute inverse FFT to reconstruct filtered signal as n-dim array
            emg_ifft = scipy.fft.irfft(emg_fft * self.mask[: self.fft_dim])
            rolling_window_size = round(x[-1])
            # Rectify flattened iFFT, construct window for smoothing and
            # calculate exponential moving average using FFT convolution
            emg_abs = np.abs(emg_ifft.ravel())
            window_ema = np.array(
                [(1 - x[-2]) ** i for i in range(rolling_window_size)]
            )
            window_ema = window_ema / window_ema.sum()
            emg_ema = scipy.signal.fftconvolve(emg_abs, window_ema[::-1], mode="valid")
            grip_flat = grip_dat.ravel()[rolling_window_size - 1 :]

            axs[row, col].plot(
                time_dat.ravel()[rolling_window_size - 1 :], emg_ema, color="black"
            )
            axs[row, col].set_xticks([])  # Remove x-axis labels
            axs[row, col].set_yticks([])  # Remove y-axis labels
            axs[row, col].set_title("Clean EMG Data", fontsize=9)

            axs[row + 1, col].plot(
                time_dat.ravel()[rolling_window_size - 1 :], grip_flat, color="black"
            )
            axs[row + 1, col].set_xticks([])  # Remove x-axis labels
            axs[row + 1, col].set_yticks([])  # Remove y-axis labels
            axs[row + 1, col].set_title("Grip Force Data", fontsize=9)

            # Calculate max_corr
            emg_ema = emg_ema - emg_ema.mean()
            grip_flat = grip_flat - grip_flat.mean()

            corr_ema = scipy.signal.fftconvolve(emg_ema, grip_flat[::-1], mode="full")
            corr_ema /= len(grip_flat) * emg_ema.std() * grip_flat.std()

            max_corr = np.abs(corr_ema).max()

            # Add max_corr as a label to the second subplot
            axs[row, col].text(
                0.5,
                -0.07,
                f"{file_name} - Max Corr: {max_corr*100:.2f}%",
                ha="center",
                va="top",
                transform=axs[row, col].transAxes,
                fontsize=9,
            )

            i += 1
            col += 1
            if col >= 6:
                col = 0
                row += 2
            if row >= 10:
                break

        # Adjust spacing between subplots
        plt.subplots_adjust(hspace=0.3)
        plt.show()

    def correlations(self, x):
        self.mask[self.optim_ind] = [x[:9]]
        corrs = []

        for emg_fft, grip_dat, time_dat in zip(
            self.emg_fft, self.grip_dat, self.time_dat
        ):
            # Compute inverse FFT to reconstruct filtered signal as n-dim array
            emg_ifft = scipy.fft.irfft(emg_fft * self.mask[: self.fft_dim])
            rolling_window_size = round(x[-1])
            # Rectify flattened iFFT, construct window for smoothing and
            # calculate exponential moving average using FFT convolution
            emg_abs = np.abs(emg_ifft.ravel())
            window_ema = np.array(
                [(1 - x[-2]) ** i for i in range(rolling_window_size)]
            )
            window_ema = window_ema / window_ema.sum()
            emg_ema = scipy.signal.fftconvolve(emg_abs, window_ema[::-1], mode="valid")
            grip_flat = grip_dat.ravel()[rolling_window_size - 1 :]

            # Calculate max_corr
            emg_ema = emg_ema - emg_ema.mean()
            grip_flat = grip_flat - grip_flat.mean()

            corr_ema = scipy.signal.fftconvolve(emg_ema, grip_flat[::-1], mode="full")
            corr_ema /= len(grip_flat) * emg_ema.std() * grip_flat.std()

            max_corr = np.abs(corr_ema).max()

            corrs.append(max_corr)

        return corrs
