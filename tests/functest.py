import pandas as pd
import folium
from scipy.optimize import minimize
from geopy.distance import geodesic
from obspy import read
import matplotlib.pyplot as plt

# Leer los datos del seismograma
st = read('data/ee140418.o27')
st.plot()

# Extraer los datos de los componentes
z_data = st[0].data
n_data = st[1].data
e_data = st[2].data

# Plotear los componentes del seismograma
fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
axs[0].plot(z_data, color='black')
axs[0].set_title('Componente Z')
axs[1].plot(n_data, color='black')
axs[1].set_title('Componente N')
axs[2].plot(e_data, color='black')
axs[2].set_title('Componente E')
plt.xlabel('Tiempo')
plt.show()

import numpy as np

def autopicking(z_data, n_data, e_data, threshold_p=0.5, threshold_s=0.3):
    """
    Identifica los arribos de las ondas P y S en un seismograma.
    
    :param z_data: Datos del componente Z.
    :param n_data: Datos del componente N.
    :param e_data: Datos del componente E.
    :param threshold_p: Umbral para la detección de la onda P.
    :param threshold_s: Umbral para la detección de la onda S.
    :return: Indices de los arribos de las ondas P y S.
    """
    # Normalizar los datos
    z_norm = z_data / np.max(np.abs(z_data))
    n_norm = n_data[:len(z_norm)] / np.max(np.abs(n_data[:len(z_norm)]))
    e_norm = e_data[:len(z_norm)] / np.max(np.abs(e_data[:len(z_norm)]))

    # Identificación del arribo de la onda P
    p_index = np.argmax(z_norm > threshold_p)
    
    # Identificación del arribo de la onda S
    # Usamos la suma de los componentes horizontales para detectar la onda S
    horizontal_sum = np.abs(n_norm) + np.abs(e_norm)
    s_index = p_index + np.argmax(horizontal_sum[p_index:] > threshold_s)

    return p_index, s_index

# Usar la función de autopicking en los datos del seismograma
p_index, s_index = autopicking(z_data, n_data, e_data)

# Mostrar los resultados
print(f"Arribo de la onda P: {p_index}")
print(f"Arribo de la onda S: {s_index}")

# Plotear los resultados
fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
axs[0].plot(z_data, color='black')
axs[0].axvline(p_index, color='red', linestyle='--', label='P Arrival')
axs[0].set_title('Componente Z')
axs[1].plot(n_data, color='black')
axs[1].axvline(s_index, color='blue', linestyle='--', label='S Arrival')
axs[1].set_title('Componente N')
axs[2].plot(e_data, color='black')
axs[2].axvline(s_index, color='blue', linestyle='--', label='S Arrival')
axs[2].set_title('Componente E')
plt.xlabel('Tiempo')
plt.show()