"""
SISMOLAB:
Aplicación para el análisis de seismogramas y cálculo de magnitudes de momento y energía radiada.
Desarrollada por:
- José María García Márquez
- Email: josemariagarciamarquez2.72@gmail.com
- Github: https://www.github.com/JoseMariaGarciaMarquez
- LinkedIn: https://www.linkedin.com/in/josé-maría-garcía-márquez-556a75199/
- Webpage: https://www.josemaria.me
- PayPal: https://www.paypal.com/paypalme/Chemitas96
- Patreon: https://patreon.com/chemitas

Esta aplicación permite cargar datos de seismogramas, visualizarlos y realizar cálculos de distancia, magnitud de momento y magnitud de energía radiada. Los resultados se guardan en un archivo CSV.

Para utilizar la aplicación, siga los siguientes pasos:
1. Haga clic en el botón "Load Data" para cargar un archivo de seismograma.
2. Utilice el botón "Toggle Picking Mode" para activar o desactivar el modo de selección.
3. Haga clic en el seismograma para seleccionar un punto de interés y guardar la amplitud y el tiempo en el archivo CSV.
4. Utilice el botón "Calculate Distance" para calcular la distancia basada en los tiempos P y S guardados en el archivo CSV.
5. Utilice los botones "Calculate Moment Magnitude" y "Calculate Energy Magnitude" para calcular las magnitudes de momento y energía radiada, respectivamente, utilizando los datos del seismograma y el valor de amplitud de referencia (A0) guardado en el archivo CSV.

Nota: Asegúrese de tener instaladas las siguientes bibliotecas antes de ejecutar la aplicación:
- obspy
- numpy
- pandas
- tkinter
- matplotlib

Para más información y detalles sobre el funcionamiento de la aplicación, consulte la documentación en el repositorio de GitHub mencionado anteriormente.

"""

import csv
import math
import numpy as np
import pandas as pd
import tkinter as tk
from obspy import read
from tkinter import filedialog, messagebox, simpledialog, ttk
import matplotlib.pyplot as plt
import sys

# Coordenadas de las estaciones de referencia
CUIG = [19.329, -99.178]  # El Pozo (Ciudad Universitaria, Coyoacán, Ciudad de México)
CAIG = [17.048, -100.267]  # El Cayaco (Coyuca de Benítez, Guerrero)
CMIG = [17.091, -94.884]   # Matías Romero (Oaxaca)
HUIG = [15.769, -96.108]   # Huatulco (Oaxaca)
OXIG = [17.072, -96.733]   # Oaxaca (Oaxaca)
PLIG = [18.392, -99.502]   # Platanillo (Iguala, Guerrero)
LPIG = [24.101, -110.309]  # La Paz (Baja California Sur)
SRIG = [27.32, -112.241]   # Santa Rosalía (Baja California Sur)
PPIG = [19.067, -98.628]   # Popocatépetl


output_path = "/Users/Chemitas/Desktop/Desk/proyectos/seislab/picking.csv"

class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.configure(state="normal")
        self.widget.insert("end", str, (self.tag,))
        self.widget.configure(state="disabled")
        self.widget.see("end")

    def flush(self):
        pass

class SeismogramApp:
    def __init__(self, root):
        self.root = root
        self.root.title("SISMOLAB")
        
        self.st = None
        self.picking_mode = False

        # Crear los botones
        self.load_button = ttk.Button(self.root, text="Load Data", command=self.load_data)
        self.load_button.pack(pady=10)

        self.picking_button = ttk.Button(self.root, text="Toggle Picking Mode", command=self.toggle_picking)
        self.picking_button.pack(pady=10)

        self.distance_button = ttk.Button(self.root, text="Calculate Distance", command=self.calculate_distance)
        self.distance_button.pack(pady=10)

        self.magnitude_w_button = ttk.Button(self.root, text="Calculate Moment Magnitude", command=lambda: self.calculate_magnitude('w'))
        self.magnitude_w_button.pack(pady=10)

        self.magnitude_e_button = ttk.Button(self.root, text="Calculate Energy Magnitude", command=lambda: self.calculate_magnitude('e'))
        self.magnitude_e_button.pack(pady=10)

 
        self.text_terminal = tk.Text(self.root, height=10, state="disabled")
        self.text_terminal.pack(pady=10)


        sys.stdout = TextRedirector(self.text_terminal, "stdout")

    def load_data(self):
        ruta = filedialog.askopenfilename()
        try:
            self.st = read(ruta)
            self.seismogram()
        except Exception as e:
            messagebox.showerror("Error", "Error al cargar el archivo\n" + str(e))

    def toggle_picking(self):
        self.picking_mode = not self.picking_mode
        status = "activado" if self.picking_mode else "desactivado"
        messagebox.showinfo("Modo Picking", f"Modo Picking {status}")

    def seismogram(self):
        fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        components = ['Z', 'N', 'E']
        for i, tr in enumerate(self.st):
            axs[i].plot(tr.times("utcdatetime"), tr.data, c='black', linewidth=0.5)
            axs[i].set_ylabel("Amplitud")
            axs[i].set_title(f"Componente {components[i]}")
            axs[i].grid()
        axs[2].set_xlabel("Tiempo [s]")
        fig.suptitle(f"Estación: {self.st[0].stats.station}\n Start Time: {self.st[0].stats.starttime}\n End Time: {self.st[0].stats.endtime}")
        fig.canvas.mpl_connect('button_press_event', lambda event: self.click_event(event, fig))
        plt.show()

    def click_event(self, event, fig):
        if fig.canvas.toolbar.mode != '': 
            return 

        if not self.picking_mode:
            return 

        if event.inaxes is not None:
            ax = event.inaxes
            index = ax.get_subplotspec().rowspan.start
            tr = self.st[index]
            times = tr.times("utcdatetime")
            time_diff = np.abs(times - event.xdata)
            idx = np.argmin(time_diff)
            
            if 0 <= idx < len(tr.data):
                amplitude = tr.data[idx]
                time = times[idx]
                self.save_to_csv(amplitude, time)
        else:
            print("Click was outside the axes")

    def save_to_csv(self, amplitude, time):
        try:
            data = {'Amplitude': [amplitude], 'Time': [time]}
            df = pd.DataFrame(data)
            df.to_csv(output_path, mode='a', header=not pd.io.common.file_exists(output_path), index=False)
        except Exception as e:
            print(f"Error saving data to CSV: {e}")

    def distance(self, csv_path):
        """
        Calcula la distancia basada en los tiempos P y S a partir de un archivo CSV.
        :param csv_path: ruta del archivo CSV
        :return: distancia calculada
        """
        try:
            df = pd.read_csv(csv_path)
            
            if len(df) < 2:
                print("El archivo CSV debe contener al menos dos registros de tiempo.")
                return None

            df['Time'] = pd.to_datetime(df['Time'])

            df = df.sort_values(by='Time')

            p_time = df['Time'].iloc[0]
            s_time = df['Time'].iloc[1]

            dt = (s_time - p_time).total_seconds()
            print(f"DT = {dt}[s]")

            if dt <= 0:
                print("El tiempo S debe ser posterior al tiempo P.")
                return None

            vp = 5.5
            vs = 0.7 * vp
            
            dslow = 1/vs - 1/vp

            distance = abs(dt / dslow)
            print(f"Distancia: {distance} km")
            return distance
        except Exception as e:
            print(f"Error al calcular la distancia: {e}")
            return None

    def magnitude_w(self, A0):
        z_data = self.st[0].data
        n_data = self.st[1].data
        e_data = self.st[2].data

        amp_z = np.max(z_data)
        amp_n = np.max(n_data)
        amp_e = np.max(e_data)
        
        A = np.sqrt(amp_z**2 + amp_n**2 + amp_e**2)
        
        print(f"A = {A}")
        print(f"A0 = {A0}")
        M0 = self.momentum0(A, A0) # N m
        print(f"M0 = {M0} N m")
        MW = (2/3) * (math.log10(M0) - 9.1)
        
        return MW

    def momentum0(self, A, A0):
        M0_dyne_cm = (A/A0) * 10**23  # dyne-cm
        M0_Nm = M0_dyne_cm * 1e-7  # Convertir a N m
        return M0_Nm

    def magnitude_energy(self, A0):
        Es = self.energy_SandO(A0)
        ME = 2/3 * math.log10(Es) - 8.45
        print("Energy: ", Es)
        print("Magnitude Energy: ", ME)
        return ME

    def energy_SandO(self, A0):
        Es = 10**(1.5 * self.magnitude_w(A0) + 12.66)
        return Es

    def calculate_distance(self):
        try:
            dist = self.distance(output_path)
            if dist is not None:
                messagebox.showinfo("Distancia Calculada", f"Distancia: {dist} km")
        except Exception as e:
            messagebox.showerror("Error", f"Error al calcular la distancia\n{e}")

    def calculate_magnitude(self, type):
        try:
            df = pd.read_csv(output_path)
            if len(df) == 0:
                messagebox.showerror("Error", "El archivo CSV no contiene datos.")
                return

            A0 = df['Amplitude'][0]

            if self.st is None:
                messagebox.showerror("Error", "No se han cargado datos del seismograma.")
                return

            if type == 'w':
                magnitude = self.magnitude_w(A0)
                messagebox.showinfo("Magnitud de Momento", f"Magnitud de Momento: {magnitude}")
            elif type == 'e':
                magnitude = self.magnitude_energy(A0)
                messagebox.showinfo("Magnitud de Energía", f"Magnitud de Energía: {magnitude}")
            else:
                messagebox.showerror("Error", "Tipo de magnitud no válido.")
        except Exception as e:
            messagebox.showerror("Error", f"Error al calcular la magnitud\n{e}")

def main():
    root = tk.Tk()
    app = SeismogramApp(root)
    root.mainloop()

main()

