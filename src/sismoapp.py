import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from obspy import read
import math
from geopy.distance import geodesic
import folium
from scipy.optimize import minimize

output_path = "picks/picking.csv"
distance_output_path = "files/distances.csv"

# Coordenadas de las estaciones de referencia
stations_coords = {
    'CUIG': [19.329, -99.178],  # El Pozo (Ciudad Universitaria, Coyoacán, Ciudad de México)
    'CAIG': [17.048, -100.267],  # El Cayaco (Coyuca de Benítez, Guerrero)
    'CMIG': [17.091, -94.884],   # Matías Romero (Oaxaca)
    'HUIG': [15.769, -96.108],   # Huatulco (Oaxaca)
    'OXIG': [17.072, -96.733],   # Oaxaca (Oaxaca)
    'PLIG': [18.392, -99.502],   # Platanillo (Iguala, Guerrero)
    'LPIG': [24.101, -110.309],  # La Paz (Baja California Sur)
    'SRIG': [27.32, -112.241],   # Santa Rosalía (Baja California Sur)
    'PPIG': [19.067, -98.628]    # Popocatépetl
}

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
        self.gain = 1.0  # Inicializar el gain en 1.0

        self.control_frame = ttk.Frame(self.root)
        self.control_frame.grid(row=0, column=0, padx=10, pady=10)

        self.plot_frame = ttk.Frame(self.root)
        self.plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

        self.load_button = ttk.Button(self.control_frame, text="Load Data", command=self.load_data)
        self.load_button.pack(pady=10)

        self.picking_button = ttk.Button(self.control_frame, text="Toggle Picking Mode", command=self.toggle_picking)
        self.picking_button.pack(pady=10)

        self.distance_button = ttk.Button(self.control_frame, text="Calculate Distance", command=self.calculate_distance)
        self.distance_button.pack(pady=10)

        self.magnitude_w_button = ttk.Button(self.control_frame, text="Calculate Moment Magnitude", command=lambda: self.calculate_magnitude('w'))
        self.magnitude_w_button.pack(pady=10)

        self.magnitude_e_button = ttk.Button(self.control_frame, text="Calculate Energy Magnitude", command=lambda: self.calculate_magnitude('e'))
        self.magnitude_e_button.pack(pady=10)

        self.map_button = ttk.Button(self.control_frame, text="Generar Mapa", command=self.generate_map)
        self.map_button.pack(pady=10)

        self.text_terminal = tk.Text(self.control_frame, height=10, state="disabled")
        self.text_terminal.pack(pady=10)

        self.gain_frame = ttk.Frame(self.plot_frame)
        self.gain_frame.pack(side=tk.BOTTOM, pady=10)

        self.increase_gain_button = ttk.Button(self.gain_frame, text="+", command=self.increase_gain)
        self.increase_gain_button.pack(side=tk.LEFT, padx=5)

        self.decrease_gain_button = ttk.Button(self.gain_frame, text="-", command=self.decrease_gain)
        self.decrease_gain_button.pack(side=tk.LEFT, padx=5)

        sys.stdout = TextRedirector(self.text_terminal, "stdout")

    def load_data(self):
        ruta = filedialog.askopenfilename()
        try:
            self.st = read(ruta)
            self.plot_seismogram()
        except Exception as e:
            messagebox.showerror("Error", "Error al cargar el archivo\n" + str(e))

    def toggle_picking(self):
        self.picking_mode = not self.picking_mode
        status = "activado" if self.picking_mode else "desactivado"
        messagebox.showinfo("Modo Picking", f"Modo Picking {status}")

    def plot_seismogram(self):
        for widget in self.plot_frame.winfo_children():
            if widget != self.gain_frame:
                widget.destroy()

        self.fig, self.axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        components = ['Z', 'N', 'E']
        for i, tr in enumerate(self.st):
            self.axs[i].plot(tr.times("utcdatetime"), self.gain * tr.data, c='black', linewidth=0.5)
            self.axs[i].set_ylabel("Amplitud")
            self.axs[i].set_title(f"Componente {components[i]}")
            self.axs[i].grid()
        self.axs[2].set_xlabel("Tiempo [s]")
        self.fig.suptitle(f"Estación: {self.st[0].stats.station}\n Start Time: {self.st[0].stats.starttime}\n End Time: {self.st[0].stats.endtime}")
        self.fig.canvas.mpl_connect('button_press_event', self.click_event)

        canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def click_event(self, event):
        if not self.fig.canvas.toolbar or self.fig.canvas.toolbar.mode != '': 
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
                amplitude = self.gain * tr.data[idx]
                time = times[idx]
                self.save_to_csv(amplitude, time)
                self.mark_picking_point(ax, time, amplitude)
        else:
            print("Click was outside the axes")

    def mark_picking_point(self, ax, time, amplitude):
        ax.plot(time, amplitude, 'ro')  # Marca el punto de picking en la gráfica
        self.fig.canvas.draw()  # Actualiza la gráfica

    def save_to_csv(self, amplitude, time):
        try:
            data = {'Amplitude': [amplitude], 'Time': [time]}
            df = pd.DataFrame(data)
            df.to_csv(output_path, mode='a', header=not pd.io.common.file_exists(output_path), index=False)
        except Exception as e:
            print(f"Error saving data to CSV: {e}")

    def distance(self, csv_path):
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
            
            # Guardar distancia en CSV
            data = {'distance': [distance], 'station': [self.st[0].stats.station]}
            df_distance = pd.DataFrame(data)
            df_distance.to_csv(distance_output_path, mode='a', header=not pd.io.common.file_exists(distance_output_path), index=False)
            
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

    def generate_map(self):
        try:
            maping(distance_output_path)
            messagebox.showinfo("Mapa Generado", "El mapa se ha guardado en 'examples/maping.html'")
        except Exception as e:
            messagebox.showerror("Error", f"Error al generar el mapa\n{e}")

    def increase_gain(self):
        self.gain *= 1.2
        self.plot_seismogram()

    def decrease_gain(self):
        self.gain /= 1.2
        self.plot_seismogram()

# Función para añadir círculos y marcadores
def add_circle_and_marker(location, radius, color, popup, map_obj):
    folium.Marker(location=location, popup=popup).add_to(map_obj)
    folium.Circle(
        location=location,
        radius=radius * 1000,  # Convertir km a metros
        color=color,
        fill=True,
        fill_color=color
    ).add_to(map_obj)

# Función objetivo para minimizar
def objective_function(point, stations, radii):
    return sum(abs(geodesic(point, stations[i]).km - radii[i]) for i in range(len(stations)))

def maping(csv_file):
    df = pd.read_csv(csv_file)
    stations = df['station']
    distances = df['distance']

    # Crear un mapa centrado en México
    m = folium.Map(location=[23.6345, -102.5528], zoom_start=5)

    station_locations = []
    station_radii = []

    for station, distance in zip(stations, distances):
        if station in stations_coords:
            location = stations_coords[station]
            add_circle_and_marker(location, distance, 'crimson', station, m)
            station_locations.append(location)
            station_radii.append(distance)
        else:
            print(f'Estación {station} no encontrada')

    if len(station_locations) > 1:
        # Punto inicial para la optimización (promedio de las coordenadas)
        initial_point = [
            sum(coord[0] for coord in station_locations) / len(station_locations),
            sum(coord[1] for coord in station_locations) / len(station_locations)
        ]

        # Optimización para encontrar el punto de intersección
        result = minimize(objective_function, initial_point, args=(station_locations, station_radii), method='Nelder-Mead')
        epicenter = result.x

        # Agregar marcador del epicentro
        folium.Marker(location=epicenter, popup=f'Epicentro: {epicenter[0]:.4f}, {epicenter[1]:.4f}', icon=folium.Icon(icon='star', color='orange')).add_to(m)

    # Guardar el mapa con los círculos
    m.save('examples/maping.html')
    print("Map saved as examples/maping.html")

# Llamar a la función principal
def main():
    root = tk.Tk()
    root.grid_columnconfigure(1, weight=1)
    root.grid_rowconfigure(0, weight=1)
    app = SeismogramApp(root)
    root.mainloop()

main()
