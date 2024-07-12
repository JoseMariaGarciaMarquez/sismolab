import numpy as np
import matplotlib.pyplot as plt
from obspy import read

st = read('data/ee140418.o27')
z_data = st[0].data
z_time = st[0].times('utcdatetime')
z_sample_rate = st[0].stats.sampling_rate


window_size = 5
window_step = window_size

num_windows = int((len(z_data) - window_size) / window_step) + 1

window_means = []
window_stds = []
clasification = []

# Found wave
found_wave = False
p_wave_time = None
p_wave_amplitude = None

for i in range(num_windows):
    start = i * window_step
    end = start + window_size

    start2 = start + window_step
    end2 = start2 + window_size

    if end >= len(z_data):
        break

    window = z_data[start:end]
    window_next = z_data[start2:end2]

    mean_w = np.mean(window)
    std_w = np.std(window)

    mean_next = np.mean(window_next)
    std_next = np.mean(window_next)

    n = 3

    if abs(mean_next) < n * abs(mean_w):
        valor = 0
        clasification.append(valor)



    if abs(mean_next) > n * abs(mean_w):
        valor = 1

        clasification.append(valor)
        diferencias = []

        if not found_wave:
            found_wave = True
            p_wave_amplitude = np.max(window_next)
            max_index = np.argmax(window_next)
            p_wave_time = z_time[start2 + max_index]

            fig, ax = plt.subplots()
            ax.plot(window_next)
            ax.set_title('Significant window with P wave')
            plt.show()

    window_means.append(mean_w)
    window_stds.append(std_w)

window_means = np.array(window_means)
window_stds = np.array(window_stds)
clasification_array = np.array(clasification)

print("El valor de la onda P {}\n Se encontro a las {} ".format(p_wave_amplitude, p_wave_time))
"""
fig, ax = plt.subplots(4, 1, figsize=(10, 6))
ax[0].plot(z_time, z_data, label='Data')
if p_wave_time and p_wave_amplitude:
    ax[0].plot(p_wave_time, p_wave_amplitude, 'ro', label='P wave')
ax[1].plot(window_means, label='Mean')
ax[2].plot(window_stds, label='Standard Deviation')
ax[3].plot(clasification_array, label='Classification')

ax[3].set_xlabel('Window')
ax[3].set_ylabel('Classification')
ax[3].legend()

plt.show()
"""
fig, ax = plt.subplots(dpi=200, sharex=True)
ax.plot(z_time, z_data, label='Data')
if p_wave_time and p_wave_amplitude:
    ax.plot(p_wave_time, p_wave_amplitude, 'ro', label='P wave\n Time: {}'.format(p_wave_time))
ax.set_xlabel('Time UTC')
ax.set_ylabel('Amplitude')
ax.legend()
ax.grid()
plt.show()