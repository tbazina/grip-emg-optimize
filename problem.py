#!/usr/bin/python

import pygmo as pg
import numpy as np
import numba as nb
import scipy.fft
import scipy.signal
import glob
import os
import itertools


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


class SignalProcessingParams:
  """Problem Class
  """
  def __init__(self, group_number):
    self.group_number = group_number
    self.emg_fft = None
    self.fft_dim = None
    self.sampling_rate = None
    self.window_size = None
    self.grip_dat = None
    
    self.data_loader()


  def data_loader(self, group_number):
    time_dat = []
    emg_dat = []
    grip_dat = []
        
    directory = "new_data"
    pattern = os.path.join(directory, f"???{self.group_number}??_sl.csv")
            
    for filename in glob.glob(pattern):
      data = np.genfromtxt(filename, delimiter=",", skip_header=1)
      time_dat.append(data[:, 0])
      emg_dat.append(data[:, 1])
      grip_dat.append(data[:, 2])
      
    # Combine data from multiple files into a single array
    time_dat = np.concatenate(time_dat)
    emg_dat = np.concatenate(emg_dat)
    grip_dat = np.concatenate(grip_dat)
    
    # Data length
    dat_len = time_dat.shape[0]

    # Calculate sampling rate from time (remove first and last 10 data points) 
    # Using median instead of mean - measurement imperfections
    self.sampling_rate = 1. / np.median( np.diff(time_dat[10:-10]))

    # Calculate window size for FFT - always round to nearest even number
    window_size_relative = 0.5                                     
    self.window_size = round(window_size_relative * self.sampling_rate)
    self.window_size += (self.window_size % 2)

    # FFT resolution
    fft_resolution = self.sampling_rate / self.window_size

    # Reshape data to window_size × n array (disregard starting few points)
    time_dat = time_dat[dat_len%self.window_size:].reshape((-1, self.window_size))
    emg_dat = emg_dat[dat_len%self.window_size:].reshape((-1, self.window_size)) 
    grip_dat = grip_dat[dat_len%self.window_size:].reshape((-1, self.window_size)) 
    
    # Calculate RFFT on emg_data
    self.emg_fft = scipy.fft.rfft(emg_dat)

    # Calculate RFFT frequencies and stack them
    emg_fft_freq = scipy.fft.rfftfreq(n=self.window_size, d=1./self.sampling_rate)
    # self.emg_fft_freq = np.vstack(
    #   [self.emg_fft_freq for _ in range(self.dat_len//self.window_size)])
    
    # FFT problem dimension
    self.fft_dim = emg_fft_freq.shape[0]

    # Smoothing window dimension (size + decay factor for exponential smoothing)
    smooth_dim = 2

    # Problem dimension
    self.dim = self.fft_dim + smooth_dim
    
    # Store grip data
    self.grip_dat = grip_dat
  
  
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
    rolling_window_relative_size = [i/100 for i in range(1, 50)]
    smoothing_factor = [i/1000 for i in range(50)]
    max_corrs = []

    for params in itertools.product(rolling_window_relative_size, smoothing_factor):
      
      # Random 0 - 1, 250 kom
      optimization_params = np.random.random(self.fft_dim)
      optimization_params = np.append(optimization_params, params)
      optimization_params[0] = 0 #zanemariti prvu vrijednost
            
      # Shallow copy of emg fft
      emg_fft = np.copy(self.emg_fft)

      # Filtering with custom filter
      emg_fft = emg_fft * x[:self.fft_dim]

      # Compute inverse FFT to reconstruct filtered signal as n-dim array
      emg_ifft = scipy.fft.irfft(emg_fft * optimization_params[:self.fft_dim])

      # Calculate rolling window size for smoothing (must be <= FFT window size)
      rolling_window_size = round(optimization_params[self.fft_dim] * self.sampling_rate)
      if rolling_window_size > self.window_size:
        rolling_window_size  = self.window_size

      # Construct 2-dimensional array using rolling window on iFFT
      # emg_rolling = rolling_window_to_array(emg_ifft.ravel(), rolling_window_size)
      # Rectify flattened iFFT, construct window for smoothing and
      # calculate exponential moving average using FFT convolution
      emg_abs = np.abs(emg_ifft.ravel())
      window_ema = np.array([(1-optimization_params[self.fft_dim+1])**i for i in range(rolling_window_size)]) 
      window_ema = window_ema / window_ema.sum()
      emg_ema = scipy.signal.fftconvolve(emg_abs, window_ema[::-1], mode='valid')
      
      grip_flat = self.grip_dat.ravel()[rolling_window_size-1:]
      emg_ema = emg_ema - emg_ema.mean()
      grip_flat = grip_flat - grip_flat.mean()
      
      corr_ema = scipy.signal.fftconvolve(emg_ema, grip_flat[::-1], mode='full')
      corr_ema /= (len(grip_flat) * emg_ema.std() * grip_flat.std())
      
      max_corr = np.abs(corr_ema).max()
      max_corrs.append(max_corr)
    
    avg_corr = np.mean(max_corrs)
      
    return 1 - avg_corr
  
  def get_extra_info(self):
    return (
      f'\tSampling rate: {self.sampling_rate} Hz\n'
      f'\tFFT window size: {self.window_size} samples\n'
      f'\tFFT resolution: {self.fft_resolution} Hz\n'
            )
