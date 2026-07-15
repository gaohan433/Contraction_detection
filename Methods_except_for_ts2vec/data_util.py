import numpy as np
from scipy.fftpack import fft
from scipy.signal import butter,filtfilt,lfilter
import wfdb
import os
from brainflow.data_filter import DataFilter, FilterTypes, AggOperations, NoiseTypes
import math
import librosa





def get_signal(record_name, channels=None):
    '''
    record_name为文件位置，文件名不需要后缀
    channels为数组，表明取哪些通道的信号
    return:signals为 N*L数组，N为信号数量，L为信号长度
    '''
  #  record_name = file_dir+ name
    signals, fields=wfdb.rdsamp(record_name,channels=channels)
    return signals.T

    
def butter_bandpass(signal,lowcut, highcut, fs = 20, order=4, butterworth = True):
    signal_ = signal.copy()
    if butterworth:
        filter_type = FilterTypes.BUTTERWORTH.value
    else:
        filter_type = FilterTypes.BESSEL.value
    DataFilter.perform_bandpass(signal_, fs, (lowcut+highcut)/2, highcut-lowcut,order,
                                        filter_type, 0)
    
    
    return signal_


    

def re_sample(signal, origin_fs, target_fs):
    return librosa.resample(signal, orig_sr= origin_fs, target_sr = target_fs, fix=False, scale=False)  # fix=True
    
class Time:
    def __init__(self, hour, minute, second):
        self.hour = hour
        self.minute = minute
        self.second = second
    
    def get_remain_second(self):
        seconds = (59-self.minute)*60 + 60 - self.second
        return seconds
    
    def add_minute(self, add):
        self.minute += add
        while self.minute>=60:
            self.minute = self.minute -60
            self.hour += 1
    
    def add_second(self, add):
        self.second += add
        while self.second >= 60:
            self.second -= 60
            self.add_minute(1)
            
    def get_x_label(self, L, fs):
      #  hours = np.zeros((L))
        minutes = np.zeros((L))
       # seconds = np.zeros((L))
        move_time = Time(self.hour, self.minute, self.second)
        for i in range(L):
            if i != 0 and i % fs == 0:
                move_time.second +=1
                if move_time.second == 60:
                    move_time.second=0
                    move_time.minute += 1
                    if move_time.minute == 60:
                        move_time.minute=0
                        move_time.hour += 1
            minutes[i] = move_time.minute + move_time.second/60

        minutes = np.array(np.round(minutes, 0),dtype=np.int32)
        return minutes
    
    def copy(self):
        t = Time(self.hour, self.minute,self.second)
        return t
        
    
def get_total_minute(L, fs):
    return np.round(L//(fs*60) + (L%(fs*60))/(60*fs),1)




def normalization(data):
    _range = np.max(data) - np.min(data)
    return (data - np.min(data)) / _range
 

def standardization(data):
    mu = np.mean(data)#, axis=0
    sigma = np.std(data)#, axis=0
    return (data - mu) / sigma


def mean_filter(signal,size):
    signal_ = signal.copy()
    DataFilter.perform_rolling_filter(signal_, size, AggOperations.MEAN.value)
    return signal_
