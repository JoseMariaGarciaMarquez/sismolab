import numpy as np
import matplotlib.pyplot as plt

from obspy import read

def autoP_icking(z_data, z_time):
    window_size = 5
    window_step=window_size

    num_windows = int((len(z_data) - window_size) / window_step) + 1
    classification = []
    
    found_wave = False
    p_wave_time = None
    p_wave_amplitude = None

    for i in range(num_windows):
        start = i * window_step
        end = start + window_size

        if end >= len(z_data):
            break

        window = z_data[start:end]
        mean_w = np.mean(window)

        start_next = start + window_step
        end_next = start_next + window_size
        if end_next >= len(z_data):
            break

        window_next = z_data[start_next:end_next]
        mean_next = np.mean(window_next)

        n = 3

        if abs(mean_next) < n * abs(mean_w):
            classification.append(0)
        elif abs(mean_next) > n * abs(mean_w):
            classification.append(1)
            if not found_wave:
                found_wave = True
                p_wave_amplitude = np.max(window_next)
                max_index = np.argmax(window_next)
                p_wave_time = z_time[start_next + max_index]

    return p_wave_time, p_wave_amplitude

def main():
    st = read('data/ee140418.o27')
    tr = st[0]
    z_data = tr.data
    z_time = tr.times('utcdatetime')

    p_wave_time, p_wave_amplitude = autoP_icking(z_data, z_time)
    print('P-wave time:', p_wave_time)
    print('P-wave amplitude:', p_wave_amplitude)

main()