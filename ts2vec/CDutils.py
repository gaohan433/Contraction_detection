import scipy.io as sio

import matplotlib.pyplot as plt
import bottleneck as bn
import numpy as np
import ruptures as rpt

from scipy import signal
import warnings

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y


def rms_envelope(sig, frame_length):
    """
    Calculates the RMS envelope of a signal.

    Args:
        signal (numpy.ndarray): The input signal.
        frame_length (int): The length of the RMS window in samples.

    Returns:
        numpy.ndarray: The RMS envelope of the signal.
    """

    # Calculate the analytic signal using the Hilbert transform
    analytic_signal = signal.hilbert(sig)

    # Calculate the magnitude of the analytic signal (envelope)
    envelope = np.abs(analytic_signal)

    # Calculate the RMS envelope
    rms_envelope = np.sqrt(np.convolve(envelope ** 2, np.ones(frame_length), 'same') / frame_length)

    return rms_envelope


def tkeo_computation(emg_ch):
    t1 = emg_ch[1:-1]
    t2 = emg_ch[:-2]
    t3 = emg_ch[2:]
    tkeo = t1 ** 2 - t2 * t3
    return tkeo


def find_regions(binary_sequence, value):
    # Convert the binary sequence string to a numpy array of integers
    binary_array = np.array(list(map(int, binary_sequence)))

    # Identify transitions between different values
    transition_points = np.where(np.diff(binary_array) != 0)[0] + 1

    # Add starting and end points to transition points
    transitions = np.concatenate(([0], transition_points, [len(binary_array)]))

    regions = []
    for i in range(len(transitions) - 1):
        start = transitions[i]
        end = transitions[i + 1]
        if binary_array[start] == value:
            regions.append((start, end))

    return regions

def find_cpd(binary_sequence, value):
    # Convert the binary sequence string to a numpy array of integers
    binary_array = np.array(list(map(int, binary_sequence)))

    # Identify transitions between different values
    transition_points = np.where(np.diff(binary_array) != 0)[0] + 1

    # Add starting and end points to transition points
    transitions = np.concatenate(([0], transition_points, [len(binary_array)]))

    cpds = []
    for i in range(len(transitions) - 1):
        start = transitions[i]
        end = transitions[i + 1]
        if binary_array[start] == value:
            cpds.append(start)
            cpds.append(end)
    if len(cpds) ==0:
        cpds.insert(0, 0)
        cpds.append(binary_sequence.shape[0])
    if cpds[0]!=0:
        cpds.insert(0, 0)
    if cpds[-1]!=binary_sequence.shape[0]:
        cpds.append(binary_sequence.shape[0])

    return cpds

def find_closest_value(lst, target):
  return min(lst, key=lambda x: abs(x - target))

def cp_score(gt_cpds,pre_cpds,time_series_l):
    sum_diff = 0
    for i in range(len(pre_cpds)):
        value = find_closest_value(gt_cpds, pre_cpds[i])
        sum_diff = sum_diff+ np.abs(value-pre_cpds[i])

    return sum_diff/(len(pre_cpds)*time_series_l)

def cp_score1(gt_cpds,pre_cpds,time_series_l):
    sum_diff = 0
    for i in range(len(pre_cpds)):
        value = find_closest_value(gt_cpds, pre_cpds[i])
        sum_diff = sum_diff+ np.abs(value-pre_cpds[i])
    ratio_diff = np.abs(len(pre_cpds)-len(gt_cpds))/len(gt_cpds)
    # print(ratio_diff)
    return (1+ratio_diff)*sum_diff/(time_series_l)

def artifact_detection(emg_ch, global_th=6, local_th=4, window_size=41):
    thr_global = np.mean(emg_ch) + global_th * np.std(emg_ch)
    index_global = emg_ch > thr_global
    tkeo = tkeo_computation(emg_ch)
    m_sq_mean = bn.move_mean(tkeo, window_size)
    ma1 = m_sq_mean[window_size - 1:]
    hf_wind = int((window_size + 1) / 2)
    tkeo_err_adj = (tkeo[hf_wind:hf_wind + ma1.shape[0]] - ma1) / ma1
    thr_adj = np.mean(tkeo_err_adj) + local_th * np.std(tkeo_err_adj)
    index_adj = (tkeo_err_adj > thr_adj)

    index_local = np.zeros(emg_ch.shape[0],)

    index_local[hf_wind:hf_wind + ma1.shape[0]] = index_adj
    index_all = index_local + index_global
    # EMG_single_channel = EMG_ch
    index_expand = np.copy(index_all)
    af_loc = np.where(index_all >= 1)[0]

    for i in af_loc:
        l = np.max([0,i-5])
        r = np.min([i+5,emg_ch.shape[0]])
        index_expand[l:r] = 1

    # index_expand1 = np.append(index_expand, index_expand[-1])
    # index_expand2 = np.insert(index_expand, 0, index_expand[0])
    # dif_index = index_expand1 - index_expand2
    # index_left = np.where(dif_index == 1)
    # index_right = np.where(dif_index == -1)
    ar_region = find_regions(index_expand, 1)

    return ar_region, index_expand


def artifact_removal(emg_ch, ar_region, index_expand, s_rate=2.56, min_f=0.34, max_f=1, degree=50,
                     in_radius=20, l_init=0.5):
    emg_single_channel = np.copy(emg_ch)

    for i in range(len(ar_region)):
        left = ar_region[i][0]
        right = ar_region[i][1]

        signal_ori = emg_ch[left - in_radius:right + in_radius]  # Original signal

        af_index = np.where(index_expand[left - in_radius:right + in_radius] == 1)
        signal_masked = np.copy(signal_ori)
        signal_pre = np.copy(signal_ori)
        x = np.arange(1, signal_masked.shape[0] + 1) / s_rate
        x1 = np.copy(x)

        sig_masked = np.delete(signal_masked, af_index, 0)  # Masked signal
        en_masked = np.sum(np.square(sig_masked)) / (sig_masked.shape[0])


        x1 = np.delete(x1, af_index, 0)
        step = (max_f - min_f) / (degree - 1)

        y_reg = np.zeros(degree, )
        A = np.ones((x1.shape[0], degree))
        for j in range(degree - 1):
            A[:, j] = np.sin(2 * np.pi * (min_f + step * j) * x1)

        sim = np.inf
        new_sim = 0
        l = l_init
        while new_sim <= sim:
            A_reg = l * np.identity(degree)

            A_add = np.concatenate((A, A_reg), axis=0)
            y_add = np.concatenate((sig_masked, y_reg), axis=0)
            A_pinv = np.linalg.pinv(A_add)

            coef = np.dot(A_pinv, y_add)

            A_full = np.ones((x.shape[0], degree))
            for j in range(degree - 1):
                A_full[:, j] = np.sin(2 * np.pi * (min_f + step * j) * x)
            prediction = np.matmul(A_full, coef)
            pre = prediction[af_index]

            en_pre = np.sum(np.square(pre)) / (pre.shape[0])

            new_sim = abs(en_pre - en_masked)
            if new_sim < sim:
                signal_pre[af_index] = prediction[af_index]  # New signal
                emg_single_channel[left - in_radius:right + in_radius] = signal_pre
                sim = new_sim
                l = l + 0.5

    return emg_single_channel


def ar_remove_patient(emg_sig):
    emg_ar_removal = np.zeros_like(emg_sig)
    for ch in range(emg_sig.shape[0]):
        emg_ch = emg_sig[ch, :]
        ar_region, index_all = artifact_detection(emg_ch)
        emg_ar_single_ch = artifact_removal(emg_ch, ar_region, index_all)
        emg_ar_removal[ch, :] = emg_ar_single_ch
    return emg_ar_removal


def contraction_detection(sig, p):
    algo = rpt.KernelCPD(kernel="rbf").fit(sig)
    result = algo.predict(n_bkps=p)
    # rpt.display(sig, result)
    # plt.show()
    c = np.zeros(sig.shape[0])
    result.insert(0, 0)

    first_sub = sig[result[0]:result[1], :]
    ef = np.sum(np.square(first_sub)) / first_sub.shape[0]
    second_sub = sig[result[1]:result[2], :]
    es = np.sum(np.square(second_sub)) / second_sub.shape[0]

    last_sub = sig[result[-2]:result[-1], :]
    el = np.sum(np.square(last_sub)) / last_sub.shape[0]
    second_last_sub = sig[result[-3]:result[-2], :]
    e_sl = np.sum(np.square(second_last_sub)) / second_last_sub.shape[0]

    # print(result)
    if ef > es:
        c[result[0]:result[1]] = 1
    if el > e_sl:
        c[result[-2]:result[-1]] = 1
    for j in range(len(result) - 3):
        sub1 = sig[result[j]:result[j + 1], :]
        sub2 = sig[result[j + 1]:result[j + 2], :]
        sub3 = sig[result[j + 2]:result[j + 3], :]
        e1 = np.sum(np.square(sub1)) / sub1.shape[0]
        e2 = np.sum(np.square(sub2)) / sub2.shape[0]
        e3 = np.sum(np.square(sub3)) / sub3.shape[0]
        if e2 > e1 and e2 > e3:
            c[result[j + 1]:result[j + 2]] = 1
    return c

def contraction_detection_single_ch(sig, p):
    algo = rpt.KernelCPD(kernel="rbf").fit(sig)
    result = algo.predict(n_bkps=p)
    # rpt.display(sig, result)
    # plt.show()
    c = np.zeros(sig.shape[0])
    result.insert(0, 0)

    first_sub = sig[result[0]:result[1]]
    ef = np.sum(np.square(first_sub)) / first_sub.shape[0]
    second_sub = sig[result[1]:result[2]]
    es = np.sum(np.square(second_sub)) / second_sub.shape[0]

    last_sub = sig[result[-2]:result[-1]]
    el = np.sum(np.square(last_sub)) / last_sub.shape[0]
    second_last_sub = sig[result[-3]:result[-2]]
    e_sl = np.sum(np.square(second_last_sub)) / second_last_sub.shape[0]

    # print(result)
    if ef > es:
        c[result[0]:result[1]] = 1
    if el > e_sl:
        c[result[-2]:result[-1]] = 1
    for j in range(len(result) - 3):
        sub1 = sig[result[j]:result[j + 1]]
        sub2 = sig[result[j + 1]:result[j + 2]]
        sub3 = sig[result[j + 2]:result[j + 3]]
        e1 = np.sum(np.square(sub1)) / sub1.shape[0]
        e2 = np.sum(np.square(sub2)) / sub2.shape[0]
        e3 = np.sum(np.square(sub3)) / sub3.shape[0]
        if e2 > e1 and e2 > e3:
            c[result[j + 1]:result[j + 2]] = 1
    return c

def find_detection_single_ch(sig, result):

    result = result[result!= 0]

    c = np.zeros(sig.shape[0])
    if result.shape[0] == 0:
        return c
    else:
        result = np.concatenate(([0],result))
        first_sub = sig[result[0]:result[1]]
        ef = np.sum(np.square(first_sub)) / first_sub.shape[0]
        if result.shape[0] == 2:
            second_sub = sig[result[1]:]
            es = np.sum(np.square(second_sub)) / second_sub.shape[0]
            if ef > es:
                c[result[0]:result[1]] = 1
            return c
        else:
            second_sub = sig[result[1]:result[2]]
            es = np.sum(np.square(second_sub)) / second_sub.shape[0]

            last_sub = sig[result[-2]:result[-1]]
            el = np.sum(np.square(last_sub)) / last_sub.shape[0]
            second_last_sub = sig[result[-3]:result[-2]]
            e_sl = np.sum(np.square(second_last_sub)) / second_last_sub.shape[0]

            # print(result)
            if ef > es:
                c[result[0]:result[1]] = 1
            if el > e_sl:
                c[result[-2]:result[-1]] = 1
            for j in range(len(result) - 3):
                sub1 = sig[result[j]:result[j + 1]]
                sub2 = sig[result[j + 1]:result[j + 2]]
                sub3 = sig[result[j + 2]:result[j + 3]]
                e1 = np.sum(np.square(sub1)) / sub1.shape[0]
                e2 = np.sum(np.square(sub2)) / sub2.shape[0]
                e3 = np.sum(np.square(sub3)) / sub3.shape[0]
                if e2 > e1 and e2 > e3:
                    c[result[j + 1]:result[j + 2]] = 1
            return c


def zero_crossing_rate(sig):
    # Ensure numpy array for easy computation

    # Compute the sign of the data points
    sign_data = np.sign(sig)

    # Count zero crossings by comparing consecutive elements
    zero_crossings = np.where(np.diff(sign_data))[0]

    # Calculate ZCR
    zcr = len(zero_crossings) / len(sig)

    return zcr
# def calculate_metrics(y_true, y_pred):
#     """
#     Calculate Precision, Recall, and F1 Score for two binary sequences.
#
#     :param y_true: List of true binary values
#     :param y_pred: List of predicted binary values
#
#     :return: Dictionary containing Precision, Recall, and F1 Score
#     """
#     tp = sum((1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 1))  # True Positives
#     fp = sum((1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 1))  # False Positives
#     fn = sum((1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 0))  # False Negatives
#
#     # Precision: tp / (tp + fp)
#     precision = tp / (tp + fp) if (tp + fp) > 0 else 0
#
#     # Recall: tp / (tp + fn)
#     recall = tp / (tp + fn) if (tp + fn) > 0 else 0
#
#     # F1 Score: 2 * (precision * recall) / (precision + recall)
#     f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
#
#     # return {
#     #     "Precision": precision,
#     #     "Recall": recall,
#     #     "F1 Score"
#     return f1_score

# def pca_computation(sig,seg_length,overlap,seg_num):
#     '''
#     Precompute the PCA for each subsequence
#     '''
#     sig_dim = sig.shape[0]
#
#     U_array = np.zeros((seg_num,sig_dim,seg_length))
#
#     S_array = np.zeros((seg_num,seg_length))
#     for i in range(seg_num):
#         l = i * overlap
#         r = i * overlap + seg_length
#         seg = sig[:, l:r]
#         d_mean = np.mean(seg, axis=1, keepdims=True)
#         c_seg = seg - d_mean
#         U, S, VT = np.linalg.svd(c_seg/seg_length, full_matrices=False)
#         U_array[i,:,:]=U
#         S_array[i,:]=S
#     return U_array,S_array
#
#
# def pca_dist(U1, U2, S1, S2, dim=10):
#
#     dist_array = np.zeros((S1.shape[0], 1))
#     for i in range(dim):  # S1.shape[0]
#         cos_theta = np.dot(U1[:, i], U2[:, i])
#         dist_array[i, :] = np.abs(S1[i] * cos_theta - S2[i]) + np.abs(S2[i] * cos_theta - S1[i])
#     return np.sum(dist_array)
#
#
# def matrix_profile(sig, seg_length, overlap):
#     sig_length = sig.shape[1]
#     seg_num = int(np.floor((sig_length - seg_length) / overlap + 1))
#     U_array,S_array = pca_computation(sig,seg_length,overlap,seg_num)
#
#     dis_matrix = np.full((seg_num, seg_num), np.inf)
#     for i in range(seg_num):
#         l1 = i * overlap
#         r1 = i * overlap + seg_length
#         U1 = U_array[i,:,:]
#         S1 = S_array[i,:]
#
#         for j in range(seg_num):
#             l2 = j * overlap
#             r2 = j * overlap + seg_length
#             b1 = l1 <= l2 and l2 <= r1
#             b2 = l2 <= l1 and l1 <= r2
#             if b1 or b2:
#                 continue
#             else:
#                 U2 = U_array[j, :, :]
#                 S2 = S_array[j, :]
#                 dis_matrix[i, j] = pca_dist(U1, U2, S1, S2)
#     return np.min(dis_matrix, axis=1)

def matrix_profile(sig, seg_length):
    ch_num = sig.shape[0]
    rms_array = np.zeros_like(sig)
    for i in range(ch_num):
        sig_ch = sig[i, :]
        rms_array[i, :] = rms_envelope(sig_ch, seg_length)
    return np.sum(rms_array,axis=0)

def change_point_num(profile,t=0.5,d=10):
    profile_extend = np.zeros(profile.shape[0] + 2)
    profile_extend[1:profile.shape[0] + 1] = profile
    peaks, _ = signal.find_peaks(profile_extend, height=np.max(profile) * t, distance = d)
    # fig, ax = plt.subplots(figsize=(10, 4),layout='constrained')
    # ax.plot(profile_extend)
    # ax.plot(peaks, profile_extend[peaks], "x")
    # plt.yticks([0, 7])
    # plt.xticks([0, 3500])
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # plt.savefig("contraction_profile.svg", format="svg")
    cp_num = 2*len(peaks)
    #If the peak is on the boundary
    if peaks[0]==1:
        cp_num = cp_num-1
    if peaks[-1]==profile.shape[0]:
        cp_num = cp_num-1
    return cp_num

def remove_fake_contraction(contraction,signal):
    '''
    contraction: ch*signal_length
    '''
    for ch in range(contraction.shape[0]):
        one_regions = find_regions(contraction[ch, :], 1)
        baseline_signal = np.sum(np.square(-1 * (contraction[ch, :] - 1) * signal[ch, :]))
        baseline_energy = baseline_signal / (contraction.shape[1] - np.sum(contraction[ch, :]))
        for r in range(len(one_regions)):
            c_duration = one_regions[r][1] - one_regions[r][0]
            c_signal = np.sum(np.square(signal[ch, one_regions[r][0]:one_regions[r][1]]))
            c_energy = c_signal / c_duration
            if c_duration < 40 or c_duration > 300 or c_energy < 2 * baseline_energy:
                contraction[ch,one_regions[r][0]:one_regions[r][1]] = 0
    return contraction
# def pca_dist(seg1,seg2,dim=10):
#     npoints = seg1.shape[1]
#     d1_mean = np.mean(seg1, axis=1, keepdims=True)
#
#     d2_mean = np.mean(seg2, axis=1, keepdims=True)
#     c_seg1 = seg1 - d1_mean
#     c_seg2 = seg2 - d2_mean
#
#     U1, S1, VT1 = np.linalg.svd(c_seg1 / np.sqrt(npoints), full_matrices=False)
#
#     U2, S2, VT2 = np.linalg.svd(c_seg2 / np.sqrt(npoints), full_matrices=False)
#     dist_array = np.zeros([S1.shape[0], 1])
#     for i in range(dim):#S1.shape[0]
#         cos_theta = np.dot(U1[:, i], U2[:, i])
#         dist_array[i, :] = np.abs(S1[i] * cos_theta - S2[i]) + np.abs(S2[i] * cos_theta - S1[i])
#     return np.sum(dist_array)
#
# def matrix_profile(signal,seg_length,overlap):
#     sig_length = signal.shape[1]
#     seg_num = int(np.floor((sig_length-seg_length)/overlap+1))
#     dis_matrix = np.full((seg_num,seg_num),np.inf)
#     for i in range(seg_num):
#         l1 = i*overlap
#         r1 = i*overlap+seg_length
#         s1 = signal[:,l1:r1]
#         for j in range(seg_num):
#             l2=j*overlap
#             r2 = j*overlap+seg_length
#             b1 = l1<=l2 and l2<=r1
#             b2 = l2<=l1 and l1<=r2
#             if b1 or b2:
#                 continue
#             else:
#                 s2 = signal[:,l2:r2]
#                 dis_matrix[i,j] = pca_dist(s1,s2)
#     return np.min(dis_matrix,axis = 1)