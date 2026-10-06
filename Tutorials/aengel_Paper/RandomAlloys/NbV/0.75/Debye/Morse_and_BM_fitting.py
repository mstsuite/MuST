from chempy.util import periodic as prd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# USER INPUT NEEDED HERE! --------------------------------------------------------------------------------------------

# Specify the input and output file paths
input_file_path = "/jet/home/aengel/Paper/RandomAlloys/NbV/0.75/Debye/NbV_input.txt"
output_plot_path = "/jet/home/aengel/Paper/RandomAlloys/NbV/0.75/Debye/NbV_bcc.png"
output_file_path = "/jet/home/aengel/Paper/RandomAlloys/NbV/0.75/Debye/NbV_bcc.txt"

# # Or, the user is prompted to input the file paths
# input_file_path = input("Enter the path of the input file (e.g., /path/to/data.txt): ")
# output_plot_path = input("Enter the path to save the output plot (e.g., /path/to/output_plot.png): ")
# output_file_path = input("Enter the path to save the output file (e.g., /path/to/output_values.txt): ")

# Data reading -------------------------------------------------------------------------------------------------------

# Function to read system-related data from the input file
def read_system_data(file_path):
    system_data = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                # Skip comments
                continue
            if line.strip() == '':
                # Skip empty lines
                continue
            if ':' not in line:
                break
            key, value = line.strip().split(':', 1)
            system_data[key.strip()] = value.strip()
    return system_data

# Read system-related data from the input file
system_data = read_system_data(input_file_path)

# Extract system-related data from the dictionary
system = system_data.get('System', '')
variation = system_data.get('Variation', '')
elements = system_data.get('Elements', '').split()
numbers = list(map(float, system_data.get('Numbers', '').split()))

# Calculate the total number of atoms in the unit cell
number_of_atoms = sum(numbers)

# Print the system
print("System:", system, variation)

# Look up the atomic masses for the elements
atomic_masses = [prd.relative_atomic_masses[prd.symbols.index(element)] for element in elements]

# Calculate the total atomic_mass
atomic_mass = sum(weight * atomic_mass for weight, atomic_mass in zip(numbers, atomic_masses)) / number_of_atoms
print("Atomic Mass:", atomic_mass)

# Function to read data from the input file with two columns (energy and volume per atom)
def read_data(file_path):
    E_data, V_data = [], []
    with open(file_path, 'r') as file:
        for line in file:
            stripped_line = line.strip()
            if stripped_line and (stripped_line[0].isdigit() or stripped_line[0] == '-'):
                try:
                    energy, volume = map(float, stripped_line.split())
                    E_data.append(energy)
                    V_data.append(volume)
                except ValueError:
                    print(f"Error: Invalid format in line: {line}")
    return E_data, V_data

# Read the energy and volume per atom data from the input file
E_data, V_data = read_data(input_file_path)

# Print the fitting data
print("E_data =", E_data)
print("V_data =", V_data)

# Find the index of the minimum energy value
min_energy_index = np.argmin(E_data)

# Morse fitting ------------------------------------------------------------------------------------------------------

# Define the Morse function
def morse(V, A, D, lmbda, V0):
    term1 = A - 2 * D * np.exp(- lmbda * ((3 / (4 * np.pi)) ** (1 / 3)) * (V ** (1 / 3) - V0 ** (1 / 3)))
    term2 = D * np.exp(-2 * lmbda * ((3 / (4 * np.pi)) ** (1 / 3)) * (V ** (1 / 3) - V0 ** (1 / 3)))
    E = term1 + term2
    return E

# Initial guesses for the Morse parameters [A, D, lambda, r0]
p0_morse = [E_data[min_energy_index], 0, 1, V_data[min_energy_index]]

# Perform the Morse least square fit
params_morse, _ = curve_fit(morse, V_data, E_data, p0_morse, maxfev=50000)

# Extract the fitted Morse parameters
E0_morse_fit, D_fit, lambda_fit, V0_morse_fit = params_morse

# Print the fitted Morse parameters
print("Fitted Morse Parameters:")
print("E0_morse =", E0_morse_fit)
print("D =", D_fit)
print("lambda =", lambda_fit)
print("V0_morse =", V0_morse_fit)

# Birch-Murnaghan fitting --------------------------------------------------------------------------------------------

# Define the Birch-Murnaghan function
def bm(V, E0, V0, B0, B0_prime):
    term1 = (V / V0) ** (2 / 3) - 1
    term2 = 6 - 4 * (V / V0) ** (2 / 3)
    E = E0 + (9 * V0 * B0 / 16) * ((term1 ** 3) * B0_prime + term1 ** 2 * term2)
    return E

# Initial guesses for the Birch-Murnaghan fitting parameters []
p0_bm = [E_data[min_energy_index], V_data[min_energy_index], 0.01, 1]

# Perform the Birch-Murnaghan least square fit
params_bm, _ = curve_fit(bm, V_data, E_data, p0_bm, maxfev=50000)

# Extract the fitted Birch-Murnaghan parameters
E0_bm_fit, V0_bm_fit, B0_fit, B0_prime_fit = params_bm

# Print the fitted Birch-Murnaghan parameters
print("Fitted Birch-Murnaghan Parameters:")
print("E0_bm =", E0_bm_fit)
print("V0_bm =", V0_bm_fit)
print("B0 =", B0_fit)
print("B0_prime =", B0_prime_fit)

# Pressures Calculations ---------------------------------------------------------------------------------------------

# Conversion factor (from Ry/a.u.^3 to GPa)
conversion_factor = 13.6 / 0.5291 ** 3 * 160.21766208
# conversion_factor = 0.5 * 2.9421912 * 10 ** 13 / 10 ** 9  # another way of conversion, with the same result

# Calculate P as the derivative of E_morse with respect to V
def calculate_P_morse(D, lmbda, V, V0):
    term1 = D * lmbda * ((9 * np.pi * V ** 2) / 2) ** (- 1 / 3)
    term2 = np.exp(- lmbda * (3 / (4 * np.pi)) ** (1 / 3) * (V ** (1 / 3) - V0 ** (1 / 3)))
    term3 = 1 - np.exp(- lmbda * (3 / (4 * np.pi)) ** (1 / 3) * (V ** (1 / 3) - V0 ** (1 / 3)))
    P = - term1 * term2 * term3
    return P

# # Calculate P_morse
# P_morse = [
#     calculate_P_morse(D_fit, lambda_fit, v, V0_morse_fit)
#     for v in V_data
# ]

#  third-order Birch–Murnaghan isothermal equation of state
def calculate_P_bm(B0, B0_prime, V, V0):
    term1 = (3 * B0 / 2) * ((V0 / V) ** (7 / 3) - (V0 / V) ** (5 / 3))
    term2 = 1 + (3 / 4) * (B0_prime - 4) * ((V0 / V) ** (2 / 3) - 1)
    P = term1 * term2
    return P

# # Calculate P_bm
# P_bm = [
#     calculate_P_bm(B0_fit, B0_prime_fit, v, V0_bm_fit)
#     for v in V_data
# ]

# Bulk Moduli Calculations -------------------------------------------------------------------------------------------

# Define the function to calculate the Bulk Modulus at r0
def calculate_B0_morse(D, lambd, V0):
    B = (D * lambd ** 2) / ((162 * np.pi ** 2 * V0) ** (1 / 3))
    return B

# Calculate B0 at V0_morse_fit and convert to GPa and kbar
B0_morse = calculate_B0_morse(D_fit, lambda_fit, V0_morse_fit)
B0_morse_GPa = B0_morse * conversion_factor  # Convert to GPa
B0_morse_kbar = B0_morse_GPa * 10  # Convert to kbar
print(f"B0 (Morse fit) = {B0_morse_kbar:.2f} kbar")
print(f"B0 (Morse fit) = {B0_morse_GPa:.2f} GPa")

# Convert B0_bm_fit to GPa and kbar
B0_bm = B0_fit
B0_bm_GPa = B0_fit * conversion_factor  # Convert to GPa
B0_bm_kbar = B0_bm_GPa * 10  # Convert to kbar
print(f"B0 (BM fit) = {B0_bm_kbar:.2f} kbar")
print(f"B0 (BM fit) = {B0_bm_GPa:.2f} GPa")

# Plotting -----------------------------------------------------------------------------------------------------------

# Create a finer grid for plotting
V_plot = np.linspace(min(V_data), max(V_data), 100)

# Plot the Morse function and Birch–Murnaghan function
E_morse = morse(V_plot, E0_morse_fit, D_fit, lambda_fit, V0_morse_fit)  # Compute E_morse values for the plot
E_bm = bm(V_plot, E0_bm_fit, V0_bm_fit, B0_fit, B0_prime_fit)  # Compute E_bm values for the plot

# Create the plot
# fig, ax1 = plt.subplots()
fig, ax1 = plt.subplots(figsize=(20, 12))

# Plot the energy data
ax1.plot(V_data, E_data, 'bo', markersize=8, label='$E$ Data')  # Plot the E data points
ax1.plot(V_plot, E_morse, 'r-', linewidth=1.5, label='$E$ Fit (Morse)')  # Plot the Morse function
ax1.plot(V_plot, E_bm, 'm-', linewidth=1.5, label='$E$ Fit (BM)')  # Plot the Birch-Murnaghan function
ax1.set_xlabel('$V$ [(a.u.)$^3$]', fontsize=20)
ax1.set_ylabel('$E$ [Ry]', color='blue', fontsize=20)
ax1.tick_params(axis='x', labelsize=18)
ax1.tick_params(axis='y', colors='blue', labelsize=18)

# Get the current y-axis tick locations
y_ticks = ax1.get_yticks()
formatted_labels = [f'{label:.2f}' for label in y_ticks]
# Set the tick locations and labels
ax1.set_yticks(y_ticks)
ax1.set_yticklabels(formatted_labels, fontsize=18)

# Create a secondary y-axis for pressure in GPa
ax2 = ax1.twinx()

# Plot the pressure data
P_morse = calculate_P_morse(D_fit, lambda_fit, V_plot, V0_morse_fit)  # Calculate P_morse
P_morse_GPa = [p * conversion_factor for p in P_morse]  # Convert P_morse to GPa
P_bm = calculate_P_bm(B0_fit, B0_prime_fit, V_plot, V0_bm_fit)  # Calculate P_bm
P_bm_GPa = [p * conversion_factor for p in P_bm]  # Convert P_bm to GPa
# ax2.plot(V_data, P_data_GPa, 'g--', linewidth=1.5, label='P (3PV)')  # Plot P data in GPa
ax2.plot(V_plot, P_morse_GPa, 'c--', linewidth=1.5, label='$P$ (Morse)')  # Plot P_morse in GPa
ax2.plot(V_plot, P_bm_GPa, 'y--', linewidth=1.5, label='$P$ (BM)')  # Plot P_bm in GPa
ax2.set_ylabel('$P$ [GPa]', color='green', fontsize=20)
ax2.tick_params(axis='y', colors='green', labelsize=18)

# Add horizontal line at P = 0
ax2.axhline(y=0, color='black', linestyle=':', linewidth=1.5)

# Mark the point (V0_morse_fit, E0_morse_fit) with a red triangle
ax1.plot(V0_morse_fit, 
         morse(V0_morse_fit, E0_morse_fit, D_fit, lambda_fit, V0_morse_fit), 
         'r^', 
         markersize=8, 
         label='$(V_0, E_0)$ (Morse)')

# Add a vertical line at V0_morse_fit
ax1.axvline(x=V0_morse_fit, color='gray', linestyle=':', linewidth=1.5)

# Add text annotation for V0_morse_fit
ax1.annotate(f'{V0_morse_fit:.1f}', 
             (V0_morse_fit, ax1.get_ylim()[0]), 
             xytext=(-15, 10), 
             textcoords='offset points', 
             ha='center', 
             fontsize=18, 
             color='red')

# Mark the point (V0_bm_fit, E0_bm_fit) with a magenta triangle
ax1.plot(V0_bm_fit, 
         bm(V0_bm_fit, E0_bm_fit, V0_bm_fit, B0_fit, B0_prime_fit), 
         'm^', 
         markersize=8, 
         label='$(V_0, E_0)$ (BM)') 

# Add a vertical line at V0_bm_fit
ax1.axvline(x=V0_bm_fit, color='gray', linestyle=':', linewidth=1.5)

# Add text annotation for V0_bm_fit
ax1.annotate(f'{V0_bm_fit:.1f}', 
             (V0_bm_fit, ax1.get_ylim()[0]), 
             xytext=(15, 10), 
             textcoords='offset points', 
             ha='center', 
             fontsize=18, 
             color='magenta')

# Adjust the layout to prevent overlapping labels and increase spacing
plt.tight_layout(pad=2)

# Set the title of the plot
plt.title(f"Energy-Volume and Pressure-Volume Curves for {system} {variation}", fontsize=20)

# Add legend
ax1.legend(loc="upper left", bbox_to_anchor=(0.1, 1), fontsize=20)
ax2.legend(loc="upper right", bbox_to_anchor=(0.9, 1), fontsize=20)

# Adjust the layout to make space for the title
plt.subplots_adjust(top=0.95)

# Save the plot
plt.savefig(output_plot_path)

# Debye Temperatures Calculations ------------------------------------------------------------------------------------

# Define the function to calculate the Debye Temperature
def calculate_theta(V0, B, M):
    theta = 41.63 * ((3 * V0 / (4 * np.pi)) ** (1 / 6)) * (B ** (1 / 2)) * (M ** (- 1 / 2))
    return theta

# Calculate the Debye temperature evaluated at V0_morse_fit
theta0_morse = calculate_theta(V0_morse_fit, B0_morse_kbar, atomic_mass)
print(f"Debye temperature (Morse fit) = {theta0_morse:.2f} K")

# Calculate the Debye temperature evaluated at V0_bm_fit
theta0_bm = calculate_theta(V0_bm_fit, B0_bm_kbar, atomic_mass)
print(f"Debye temperature (BM fit) = {theta0_bm:.2f} K")

# Text file output ---------------------------------------------------------------------------------------------------

# Save the values to an output file
with open(output_file_path, "w") as output_file:
    output_file.write("E_data = " + str(E_data) + "\n")
    output_file.write("V_data = " + str(V_data) + "\n")
    output_file.write("Fitted Morse Parameters:\n")
    output_file.write("E0_morse = " + str(E0_morse_fit) + "\n")
    output_file.write("D = " + str(D_fit) + "\n")
    output_file.write("lambda = " + str(lambda_fit) + "\n")
    output_file.write("V0_morse = " + str(V0_morse_fit) + "\n")
    output_file.write("Fitted Birch-Murnaghan Parameters:\n")
    output_file.write("E0_bm = " + str(E0_bm_fit) + "\n")
    output_file.write("V0_bm = " + str(V0_bm_fit) + "\n")
    output_file.write("B0 = " + str(B0_fit) + "\n")
    output_file.write("B0_prime = " + str(B0_prime_fit) + "\n")
    output_file.write("B0 (Morse fit) = {:.2f} kbar\n".format(B0_morse_kbar))
    output_file.write("B0 (Morse fit) = {:.2f} GPa\n".format(B0_morse_GPa))
    output_file.write("B0 (BM fit) = {:.2f} kbar\n".format(B0_bm_kbar))
    output_file.write("B0 (BM fit) = {:.2f} GPa\n".format(B0_bm_GPa))
    output_file.write("Debye temperature (Morse fit) = {:.2f} K\n".format(theta0_morse))
    output_file.write("Debye temperature (BM fit) = {:.2f} K\n".format(theta0_bm))
