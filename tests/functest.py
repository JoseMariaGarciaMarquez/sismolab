import pandas as pd
import folium
from scipy.optimize import minimize
from geopy.distance import geodesic

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
    m.save('examples/map.html')
    print("Map saved as examples/map.html")

# Llamar a la función principal
def main():
    maping('examples/mapfile.csv')

main()
