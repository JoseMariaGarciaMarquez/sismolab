import numpy as np
from obspy import read

import matplotlib.pyplot as plt

# Leer los datos
st = read('data/ee140418.o27')
z_data = st[0].data
z_time = st[0].times("utcdatetime")

diferencias = []

# Clasificación
clasificacion = []
for i in range(len(z_data) - 1):
    next = z_data[i+1]
    valor = z_data[i]

    norm = next * (1/valor)
    diferencias.append(norm)
diferencias.append(0)

diferencias = np.array(diferencias)
z_data_norm = (z_data - np.min(diferencias)) / (np.max(diferencias) - np.min(diferencias))
dif_norm = (diferencias - np.min(diferencias)) / (np.max(diferencias) - np.min(diferencias))

"""
print(len(z_data), len(diferencias))
fig, ax = plt.subplots(dpi=200)

scatter = ax.scatter(z_time, z_data, c=diferencias, cmap='viridis', s=0.4)
cbar = plt.colorbar(scatter)
cbar.set_label('Differences')

plt.xlabel('Time')
plt.show()
"""

fig, ax = plt.subplots(2,1, dpi=200, sharex = True)

ax[0].plot(z_time, z_data, label='Original data')
ax[1].plot(z_time, diferencias, label='Differences')

plt.show()
