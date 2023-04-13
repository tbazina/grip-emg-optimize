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
    self.emg_fft = []
    self.grip_dat = []
    self.sampling_rate = None
    self.window_size = None
    self.emg_fft_freq = None
    self.fft_dim = None
    
    self.data_loader()

  def data_loader(self):
    time_dat = []
    emg_dat = []
    grip_dat = []
    window_size_relative = 0.5                                     
    # Smoothing dimension (decay factor for exponential smoothing + window size)
    self.smooth_dim = 2
        
    directory = "new_data"
    pattern = os.path.join(directory, f"???{self.group_number}??_sl.csv")
            
    for filename in glob.glob(pattern):
      data = np.genfromtxt(filename, delimiter=",", skip_header=1)
      # print(f'{data.shape=}')
      # print(f'{filename=}')
      time_dat.append(data[:, 0])
      emg_dat.append(data[:, 1])
      grip_dat.append(data[:, 2])
      
    # print(f'{len(time_dat)=}')
    # print(f'{len(emg_dat)=}')
    # print(f'{len(grip_dat)=}')
      
    for time_t, emg, grip in zip(time_dat, emg_dat, grip_dat):
      # Data length
      dat_len = time_t.shape[0]
      # print(f'{dat_len=}')
      # Calculate sampling rate from time (remove first and last 10 data points) 
      # Using median instead of mean - measurement imperfections
      self.sampling_rate = 1. / np.median( np.diff(time_t[10:-10]))
      # print(f'{self.sampling_rate=}')
      # Calculate window size for FFT - always round to nearest even number
      self.window_size = round(window_size_relative * self.sampling_rate)
      self.window_size += (self.window_size % 2)
      # print(f'{self.window_size=}')
      # FFT resolution
      self.fft_resolution = self.sampling_rate / self.window_size
      # print(f'{self.fft_resolution=}')

      # Reshape data to window_size × n array (disregard starting few points)
      time_t = time_t[dat_len%self.window_size:].reshape((-1, self.window_size))
      emg = emg[dat_len%self.window_size:].reshape((-1, self.window_size)) 
      grip = grip[dat_len%self.window_size:].reshape((-1, self.window_size)) 
      
      # Calculate RFFT on emg_data
      self.emg_fft.append(scipy.fft.rfft(emg))
      # Store grip data
      self.grip_dat.append(grip)
      # Calculate RFFT frequencies and stack them
      self.emg_fft_freq = scipy.fft.rfftfreq(n=self.window_size, d=1./self.sampling_rate)
      # print(f'{self.emg_fft_freq=}')
      # FFT problem dimension
      self.fft_dim = self.emg_fft_freq.shape[0]
      # Problem dimension
      self.dim = self.fft_dim + self.smooth_dim
      
  def get_bounds(self):
    #TODO: procjeriti ako uvijek ukloniti DC offset ili ne
    # Calculate rolling window size for smoothing (must be <= FFT window size)
    return (
      [0.] * self.fft_dim + [0.] + [1], 
      [2.] * self.fft_dim + [0.05] + [self.window_size-1]
      )

  def get_nix(self) -> int:
    return 1

  def get_name(self):
    return (
      'Optimize EMG signal processing parameters to cross-corelate with grip force'
      )

  def fitness(self, x):
    """Fitness function

    Args:
      x (np.float64): Decision vector
        [0:self.fft_dim] - Custom window for fft spectrum of signal (0-2)
        [self.fft_dim] - Rolling window size for smoothing in sec (1-0.5)
        [self.fft_dim+1] - Exponential moving average decay (0-0.05)
        - amplification - maybe only in frequency domain 0 - 2 instead 0-1
    """
    max_corrs = []

    for emg_fft, grip_dat in zip(self.emg_fft, self.grip_dat):
      # print(f'{emg_fft.shape=}')
      # Compute inverse FFT to reconstruct filtered signal as n-dim array
      emg_ifft = scipy.fft.irfft(emg_fft * x[:self.fft_dim])
      # print(f'{emg_ifft.shape=}')
      rolling_window_size = int(x[-1])
      # print(f'{rolling_window_size=}')
      # Rectify flattened iFFT, construct window for smoothing and
      # calculate exponential moving average using FFT convolution
      emg_abs = np.abs(emg_ifft.ravel())
      # print(f'{emg_abs.shape=}')
      window_ema = np.array([(1-x[-2])**i for i in range(rolling_window_size)]) 
      window_ema = window_ema / window_ema.sum()
      # print(f'{window_ema.shape=}')
      emg_ema = scipy.signal.fftconvolve(emg_abs, window_ema[::-1], mode='valid')
      # print(f'{emg_ema.shape=}')
      
      #TODO: što radi ovaj dio
      grip_flat = grip_dat.ravel()[rolling_window_size-1:]
      # print(f'{grip_flat.shape=}')
      emg_ema = emg_ema - emg_ema.mean()
      grip_flat = grip_flat - grip_flat.mean()
      
      corr_ema = scipy.signal.fftconvolve(emg_ema, grip_flat[::-1], mode='full')
      # print(f'{corr_ema.shape=}')
      corr_ema /= (len(grip_flat) * emg_ema.std() * grip_flat.std())
      # print(f'{corr_ema.shape=}')
      
      max_corr = np.abs(corr_ema).max()
      # print(f'{max_corr=}')
      max_corrs.append(max_corr)
      # print(f'{max_corrs=}')
    
    avg_corr = np.mean(max_corrs)
    # print(f'{avg_corr=}')
      
    return [1-avg_corr]
  
  def get_extra_info(self):
    return (
      f'\tSampling rate: {self.sampling_rate} Hz\n'
      f'\tFFT window size: {self.window_size} samples\n'
      f'\tFFT resolution: {self.fft_resolution} Hz\n'
            )
