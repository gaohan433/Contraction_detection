



from librosa.core import resample
import numpy as np
import math


def get_basal_tones(signal, window_size=20*60, shift=60, percentile=0.2):
    basal_tones = []

    for i in range(0, len(signal) - window_size + 1, shift):
        window = signal[i:i+window_size]
        window = sorted(window)
        base_value = np.mean(window[:int(percentile * len(window))])
        
        basal_tones.append(base_value)
       
    basal_tones = np.array(basal_tones)

    basal_tones = np.interp(list(range(len(signal))),
                            list(range(window_size//2, len(signal) - window_size//2 + 1, shift)),
                            basal_tones)
    return basal_tones


def remove_base(y, base_window_size, shift):
    base_line = get_basal_tones(y,window_size=base_window_size, shift=shift)
    y_remove = []
    for i in range(len(y)):
        tmp = y[i] - base_line[i]
        if tmp >= 0:
            y_remove.append(tmp)
        else:
            y_remove.append(0)
    y_ = np.asarray(y_remove)
    
    return y_


def get_curve(signal, window_size=60*20, shift = 20, method = None):
    y = []
    
    for i in range(0, len(signal) - window_size + 1, shift):
        window = signal[i:i+window_size]       
        y.append(method(window)) 
    y = np.interp(list(range(len(signal))),
                            list(range(window_size//2, len(signal) - window_size//2 + 1, shift)),
                            y)
    return y
