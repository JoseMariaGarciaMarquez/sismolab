""""
Sismofast: Software para el cálculo de la magnitud de energía de un sismo de manera rápida y eficiente.
Implementando la función de auto-picking para la detección de la onda P y el cálculo de la magnitud de energía.
Solo usarse en sismos mayores a 6.0 grados.
"""

import time
import numpy as np
from obspy import read


def auto_picking(z_data):
    """
    Esta función implementa el algoritmo de auto-picking para la detección de la onda P en una señal sísmica.
    Recibe como entrada los datos de la componente Z de la señal (z_data) y los tiempos correspondientes (z_time).
    Devuelve el tiempo de llegada de la onda P (p_wave_time) y su amplitud (p_wave_amplitude).
    """

    window_size = 5
    window_step = window_size

    num_windows = (len(z_data) - window_size) // window_step + 1

    found_wave = False
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
            pass
        elif abs(mean_next) > n * abs(mean_w):
            if not found_wave:
                found_wave = True
                p_wave_amplitude = np.max(window_next)

    return p_wave_amplitude

def momentum0(A, A0):
    M0 = (A / A0) * 1e16 #N-m
    return M0

def magnitude_w(st, A0):
    z_data, n_data, e_data = st[0].data, st[1].data, st[2].data
    amp_z, amp_n, amp_e = np.max(z_data), np.max(n_data), np.max(e_data)
    A = np.sqrt(amp_z**2 + amp_n**2 + amp_e**2)
    M0 = momentum0(A, A0)  # Momento sísmico en N m
    MW = (2/3) * (np.log10(M0) - 9.1)
    return MW

def magnitude_energy(st):
    tr = st[0]
    z_data = tr.data

    p_wave_amplitude = auto_picking(z_data)
    if p_wave_amplitude is not None:
        A0 = p_wave_amplitude

        Es = 10**(1.5 * magnitude_w(st, A0) + 12.66)
        ME = 2/3 * np.log10(Es) - 8.45

        return ME
    else:
        raise ValueError("No se pudo encontrar la amplitud de la onda P")

def magnitude_amplitude(st):
    z_data, n_data, e_data = st[0].data, st[1].data, st[2].data
    amp_z, amp_n, amp_e = np.max(z_data), np.max(n_data), np.max(e_data)
    A = np.sqrt(amp_z**2 + amp_n**2 + amp_e**2)
    p_wave_amplitude = auto_picking(z_data)
    if p_wave_amplitude is not None:
        A0 = p_wave_amplitude
        M0 = momentum0(A, A0)  # Momento sísmico en N m
        MA = (np.log10(M0*1e7) - 16.1) / 1.5
        return MA

def main():
    start_time = time.time() 
    st = read('data/ee140418.o27')
    magnitude = magnitude_amplitude(st)
    print(f"Magnitud de Amplitud: {round(magnitude,1)}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Tiempo de ejecución: {elapsed_time:.4f} segundos")

if __name__ == "__main__":
    main()
