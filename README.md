### SISMOLAB

**SISMOLAB** is an interactive application for analyzing seismological data using Tkinter.

#### Features

- **Load Data** 📂: Import seismogram files from your computer.
- **Interactive Visualization** 📊: Display Z, N, and E components with zoom and pan capabilities.
- **Record Amplitudes and Times** 🖱️: Use the picking mode to register clicks on graphs and save amplitudes and times to a CSV file.
- **Distance Calculation** 📏: Calculate the distance between the earthquake origin and the station using P and S wave times.
- **Magnitude Calculation** 📐: Compute Moment Magnitude (Mw) and Energy Magnitude (Me).
- **Integrated Terminal** 💻: Display messages and results within the application interface.

SISMOLAB simplifies the analysis and visualization of seismological data, aiding in the study of earthquakes and interpretation of seismographic data.

![_eef4fcc7-f5b5-4f95-9d9f-d8b7f90878df](https://github.com/JoseMariaGarciaMarquez/sismolab/assets/30852961/fc80f1ba-020a-48e5-95de-c184f4f27e02)

#### Usage

##### Main window:
The main window has buttons:
- **Load Data**: Opens a file dialog to select a seismogram file (MiniSEED format) from your computer and loads it for analysis.
- **Picking**: Toggles the picking mode, allowing you to click on the graph to record amplitude and time values.
- **Distance**: Calculates the distance between the earthquake origin and the station based on the P and S wave arrival times.
- **Mw**: Computes the Moment Magnitude (Mw) based on the recorded amplitudes.
- **ME**: Computes the Energy Magnitude (Me) based on the recorded amplitudes.
- **Maps**: Generates a map with the locations of the seismic stations and the estimated epicenter.
- **+ Gain**: Increases the gain (amplification) of the seismogram plot.
- **- Gain**: Decreases the gain (amplification) of the seismogram plot.
- **Reset Gain**: Resets the gain to the default value.

![main_window](https://github.com/JoseMariaGarciaMarquez/sismolab/assets/30852961/b38f9565-3245-4918-8754-443c261b36b9)

#### Load Data

Click the **Load Data** button to open a file dialog. Select a seismogram file in MiniSEED format. The application will load the file and display the Z, N, and E components of the seismogram in separate plots.

![load](https://github.com/JoseMariaGarciaMarquez/sismolab/assets/30852961/71cd4b37-4788-40f1-b1e7-2714cfd57a74)

#### Picking Mode

Toggle the picking mode by clicking the **Picking** button. In this mode, you can click on the seismogram plots to record the amplitude and time of specific points. The recorded data is saved to a CSV file.

![pick_p](https://github.com/JoseMariaGarciaMarquez/sismolab/assets/30852961/609d0638-6c92-4887-b952-0c2ff17cefeb)
![pick_ps](https://github.com/JoseMariaGarciaMarquez/sismolab/assets/30852961/7ca1d808-7447-412d-b719-28de94309fd0)

#### Distance Calculation

Click the **Distance** button to calculate the distance between the earthquake origin and the station. This is done by using the recorded P and S wave arrival times. The distance is calculated using the formula:

\[ \text{Distance} = \frac{\Delta t}{\left(\frac{1}{V_s} - \frac{1}{V_p}\right)} \]

where \(\Delta t\) is the difference between the P and S wave arrival times, \(V_s\) is the S wave velocity, and \(V_p\) is the P wave velocity.

#### Magnitude Calculation

- **Mw (Moment Magnitude)**: Click the **Mw** button to calculate the Moment Magnitude. The formula used is:

\[ M_w = \frac{2}{3} \left( \log_{10}(M_0) - 9.1 \right) \]

where \(M_0\) is the seismic moment calculated from the recorded amplitudes.

- **ME (Energy Magnitude)**: Click the **ME** button to calculate the Energy Magnitude. The formula used is:

\[ M_E = \frac{2}{3} \log_{10}(E_s) - 8.45 \]

where \(E_s\) is the seismic energy calculated from the recorded amplitudes.

![me](https://github.com/JoseMariaGarciaMarquez/sismolab/assets/30852961/6a5306be-5d68-4375-8d34-6fdbcb3b24ad)

#### Map Generation

Click the **Maps** button to generate a map showing the locations of the seismic stations and the estimated epicenter. The map is created using the Folium library and saved as an HTML file.

![mapgenerator](https://github.com/JoseMariaGarciaMarquez/sismolab/assets/30852961/180d840b-2799-40af-b18b-6a0b83dd1aa0)

![mapshot](https://github.com/JoseMariaGarciaMarquez/sismolab/assets/30852961/0a35b7fb-a093-4a8f-9579-6c2dcb35577f)


