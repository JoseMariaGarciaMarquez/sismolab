import numpy as np
import matplotlib.pyplot as plt
from obspy import read
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import sys
import pandas as pd
import math
from geopy.distance import geodesic
import folium
from scipy.optimize import minimize
from PIL import ImageTk, Image

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

# Función para auto-picking de la onda P
def autoP_icking(z_data, z_time):
    window_size = 5
    window_step = window_size

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

class TextRedirector(object):
    """Clase para redirigir la salida estándar a un widget de texto en la interfaz gráfica."""
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
    """This class represents the main application window."""
    def __init__(self, root):
        self.root = root
        self.root.title("SISMOLAB")
        
        self.st = None
        self.picking_mode = False
        self.gain = 1.0  
        self.lines = []  

        self.main_frame = ttk.Frame(self.root)
        self.main_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        self.control_frame = ttk.Frame(self.main_frame)
        self.control_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        self.gain_frame = ttk.Frame(self.main_frame)
        self.gain_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

        self.plot_frame = ttk.Frame(self.main_frame)
        self.plot_frame.grid(row=0, column=2, rowspan=2, padx=10, pady=10, sticky="nsew")

        self.root.grid_columnconfigure(2, weight=1)
        self.root.grid_rowconfigure(0, weight=1)

        style = ttk.Style()
        style.configure("TButton", padding=5, relief="flat", font=('Helvetica', 10), width=8)

        self.load_button = ttk.Button(self.control_frame, text="Load", command=self.load_data)
        self.load_button.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        self.picking_button = ttk.Button(self.control_frame, text="Picking", command=self.toggle_picking)
        self.picking_button.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

        self.distance_button = ttk.Button(self.control_frame, text="Distance", command=self.calculate_distance)
        self.distance_button.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        self.magnitude_w_button = ttk.Button(self.control_frame, text="Mw", command=lambda: self.calculate_magnitude('w'))
        self.magnitude_w_button.grid(row=1, column=1, padx=5, pady=5, sticky="nsew")

        self.magnitude_e_button = ttk.Button(self.control_frame, text="ME", command=lambda: self.calculate_magnitude('e'))
        self.magnitude_e_button.grid(row=2, column=0, padx=5, pady=5, sticky="nsew")

        self.map_button = ttk.Button(self.control_frame, text="Maps", command=self.generate_map)
        self.map_button.grid(row=2, column=1, padx=5, pady=5, sticky="nsew")

        self.text_terminal = tk.Text(self.control_frame, height=10, state="disabled")
        self.text_terminal.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky="nsew")

        self.increase_gain_button = ttk.Button(self.gain_frame, text="+", command=self.increase_gain)
        self.increase_gain_button.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

        self.decrease_gain_button = ttk.Button(self.gain_frame, text="-", command=self.decrease_gain)
        self.decrease_gain_button.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")

        self.reset_gain_button = ttk.Button(self.gain_frame, text="Reset", command=self.reset_gain)
        self.reset_gain_button.grid(row=2, column=0, padx=5, pady=5, sticky="nsew")

        # Load and display the image
        self.load_image("images/icon.jpeg")

        sys.stdout = TextRedirector(self.text_terminal, "stdout")

    def load_image(self, path):
        try:
            img = Image.open(path)
            img = img.resize((150, 150), Image.ANTIALIAS)
            self.photo = ImageTk.PhotoImage(img)
            self.image_label = ttk.Label(self.gain_frame, image=self.photo)
            self.image_label.grid(row=3, column=0, padx=5, pady=5, sticky="nsew")
        except Exception as e:
            messagebox.showerror("Error", f"Error al cargar la imagen\n{e}")

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
        if hasattr(self, 'fig'):
            self.fig.clear()
        else:
            self.fig, self.axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
            self.fig.canvas.mpl_connect('button_press_event', self.click_event)

        self.lines.clear()
        components = ['Z', 'N', 'E']
        for i, tr in enumerate(self.st):
            line, = self.axs[i].plot(tr.times("utcdatetime"), self.gain * tr.data, c='black', linewidth=0.5)
            self.lines.append(line)
            self.axs[i].set_ylabel("Amplitud")
            self.axs[i].set_title(f"Componente {components[i]}")
            self.axs[i].grid()

            self.axs[i].set_ylim(np.min(tr.data) * 1.2, np.max(tr.data) * 1.2)

            # Detección automática de la onda P
            p_wave_time, p_wave_amplitude = autoP_icking(tr.data, tr.times("utcdatetime"))
            if p_wave_time and p_wave_amplitude:
                self.axs[i].plot(p_wave_time, self.gain * p_wave_amplitude, 'ro', label='P Wave Arrival')
                self.axs[i].legend()

                # Guardar la onda P detectada
                self.save_to_csv(p_wave_amplitude, p_wave_time)

        self.axs[2].set_xlabel("Tiempo [s]")
        self.fig.suptitle(f"Estación: {self.st[0].stats.station}\n Start Time: {self.st[0].stats.starttime}\n End Time: {self.st[0].stats.endtime}")
        plt.show()

    def click_event(self, event):
        if self.fig.canvas.toolbar.mode != '': 
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
                print(f"Picked point at time: {time}, amplitude: {amplitude}")
                self.save_to_csv(amplitude, time)
                self.mark_picking_point(ax, time, amplitude)
        else:
            print("Click was outside the axes")

    def mark_picking_point(self, ax, time, amplitude):
        ax.plot(time, amplitude, 'ro')  # Marca el punto de picking en la gráfica
        self.fig.canvas.draw()  # Actualiza la gráfica
        print(f"Marked point at time: {time}, amplitude: {amplitude}")

    def save_to_csv(self, amplitude, time):
        try:
            data = {'Amplitude': [amplitude], 'Time': [time]}
            df = pd.DataFrame(data)
            df.to_csv(output_path, mode='a', header=not pd.io.common.file_exists(output_path), index=False)
            print(f"Saved point to CSV: time: {time}, amplitude: {amplitude}")
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
        for i, tr in enumerate(self.st):
            self.lines[i].set_ydata(self.gain * tr.data)
        self.fig.canvas.draw()

    def decrease_gain(self):
        self.gain /= 1.2
        for i, tr in enumerate(self.st):
            self.lines[i].set_ydata(self.gain * tr.data)
        self.fig.canvas.draw()

    def reset_gain(self):
        self.gain = 1.0
        for i, tr in enumerate(self.st):
            self.lines[i].set_ydata(self.gain * tr.data)
        self.fig.canvas.draw()

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
        if (location := stations_coords.get(station)) is not None:
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
    app = SeismogramApp(root)
    root.mainloop()

main()
