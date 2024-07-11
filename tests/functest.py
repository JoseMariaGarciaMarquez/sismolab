import numpy as np
import matplotlib.pyplot as plt
from obspy import read

st = read('data/ee140418.o27')
z_data = st[0].data
z_sample_rate = st[0].stats.sampling_rate

window_size = 5
window_step = 1

num_windows = int((len(z_data) - window_size) / window_step) + 1

window_means = []
window_stds = []

for i in range(num_windows):
    start = int(i * window_step * z_sample_rate)
    end = start + int(window_size * z_sample_rate)

    if end >= len(z_data):
        break

    window = z_data[start:end]
    window_means.append(np.mean(window))
    window_stds.append(np.std(window))

window_means = np.array(window_means)
window_stds = np.array(window_stds)

fig, ax = plt.subplots(3,1,figsize=(10, 6))
ax[0].plot(z_data, label='Data')
ax[1].plot(window_means, label='Mean')
ax[2].plot(window_stds, label='Standard Deviation')

plt.show()
