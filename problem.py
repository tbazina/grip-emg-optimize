#!/usr/bin/python

import pygmo as pg
import numpy as np
import numba as nb
import scipy.fft
import scipy.signal
import glob



@nb.jit(nb.float64[:,:](nb.float64[:], nb.int64), nopython=True)
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


# Učitavanje podataka s razdvajanjem na dvije grupe od 30 mjerenja s obzirom na broj grupe
# Vraća array time_dat, emg_dat i grip_dat
def data_loader(group_number):
  time_dat = []
  emg_dat = []
  grip_dat = []
        
  pattern = f"???{group_number}??_sl.csv"
        
  for filename in glob.glob(pattern):
    data = np.genfromtxt(filename, delimiter=",", skip_header=1)
    time_dat.append(data[:, 0])
    emg_dat.append(data[:, 1])
    grip_dat.append(data[:, 2])
        
  return np.array(time_dat), np.array(emg_dat), np.array(grip_dat)


class SignalProcessingParams:
  """Problem Class
  """
  def __init__(self, group_number):
    time_dat, emg_dat, grip_dat = data_loader(group_number)
    self.time_dat = time_dat
    self.emg_dat = emg_dat
    self.grip_dat = grip_dat

    # Data length
    self.dat_len = self.time_dat.shape[0]

    # Calculate sampling rate from time (remove first and last 10 data points) 
    # Using median instead of mean - measurement imperfections
    self.sampling_rate = 1. / np.median( np.diff(self.time_dat[10:-10]))
    # Calculate window size for FFT - always round to nearest even number
    self.window_size_relative = 0.5                                                  #window_size_relative uvijek 0.5
    self.window_size = round(self.window_size_relative * self.sampling_rate)
    self.window_size += (self.window_size % 2)

    # FFT resolution
    self.fft_resolution = self.sampling_rate / self.window_size

    # Reshape data to window_size × n array (disregard starting few points)
    self.time_dat = self.time_dat[self.dat_len%self.window_size:].reshape(
      (-1, self.window_size)) 
    self.emg_dat = self.emg_dat[self.dat_len%self.window_size:].reshape(
      (-1, self.window_size)) 
    self.grip_dat = self.grip_dat[self.dat_len%self.window_size:].reshape(
      (-1, self.window_size)) 
    
    # Calculate RFFT on emg_data
    self.emg_fft = scipy.fft.rfft(self.emg_dat)

    # Calculate RFFT frequencies and stack them
    self.emg_fft_freq = scipy.fft.rfftfreq(
      n=self.window_size, d=1./self.sampling_rate
      )
    # self.emg_fft_freq = np.vstack(
    #   [self.emg_fft_freq for _ in range(self.dat_len//self.window_size)])
    
    # FFT problem dimension
    self.fft_dim = self.emg_fft_freq.shape[0]
    # Smoothing window dimension (size + decay facotr for exponential smoothing)
    self.smooth_dim = 2
    # Problem dimension
    self.dim = self.fft_dim + self.smooth_dim
    return

  def get_bounds(self) -> int:
    return ([0.] * self.dim, [2.] * self.dim)

  def get_name(self):
    return (
      'Optimize EMG signal processing parameters to cross-corelate with grip force'
      )

  def fitness(self, x):
    """Fitness function

    Args:
      x (np.float64): Decision vector
        [0:self.fft_dim] - Custom window for fft spectrum of signal (0-2)
        [self.fft_dim] - Rolling window size for smoothing in sec (0-0.5)
        [self.fft_dim+1] - Exponential moving average decay (0-0.05)
        - amplification - maybe only in frequency domain 0 - 2 instead 0-1

    """
    # Shallow copy of emg fft
    emg_fft = np.copy(self.emg_fft)
    # Filtering with custom filter
    emg_fft = emg_fft * x[:self.fft_dim]
    # Compute inverse FFT to reconstruct filtered signal as n-dim array
    emg_ifft = scipy.fft.irfft(emg_fft)
    # Calculate rolling window size for smoothing (must be <= FFT window size)
    rolling_window_size = round(x[self.fft_dim] * self.sampling_rate)
    if rolling_window_size > self.window_size:
      rolling_window_size  = self.window_size
    # Construct 2-dimensional array using rolling window on iFFT
    # emg_rolling = rolling_window_to_array(emg_ifft.ravel(), rolling_window_size)
    # Rectify flattened iFFT, construct window for smoothing and
    # calculate exponential moving average using FFT convolution
    emg_abs = np.abs(emg_ifft.ravel())
    window_ema = np.array(
      [(1-x[self.fft_dim+1])**i for i in range(rolling_window_size)]
      )
    window_ema = window_ema / window_ema.sum()
    emg_ema = scipy.signal.fftconvolve(emg_abs, window_ema[::-1], mode='valid')
    
    
    # Računanje korelacije između emg_ema i grip_dat
    
    # Calculate the mean of the two arrays
    emg_ema_mean = np.mean(emg_ema)
    grip_dat_mean = np.mean(self.grip_dat)
    
    # Subtract the means from the arrays
    emg_ema_centered = emg_ema - emg_ema_mean
    grip_dat_centered = self.grip_dat - grip_dat_mean
    
    # Calculate the correlation coefficient
    corr_coef = np.sum(emg_ema_centered * grip_dat_centered) / np.sqrt(
        np.sum(emg_ema_centered**2) * np.sum(grip_dat_centered**2)
        )
    
    return (1 - corr_coef)                                                      # Fitness vraća 1 minus korelacija
  
  
  def get_extra_info(self):
    return (
      f'\tSampling rate: {self.sampling_rate} Hz\n'
      f'\tFFT window size: {self.window_size} samples\n'
      f'\tFFT resolution: {self.fft_resolution} Hz\n'
            )
