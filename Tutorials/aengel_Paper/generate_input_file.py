import numpy as np

# USER INPUT NEEDED HERE! --------------------------------------------------------------------------------------------

# Specify the path to the file to generate
input_file_path = "/home/yuqinglin/MuST/MST/Morse/Bridges2/Li/Li_input.txt"

# Specify the system, elements inside, and their corresponding numbers
system = "Li"
variation = "" # Different parameters used in calculations of the same system
elements = ["Li"]  # Elements in the system
numbers = [1]  # Specify the number of atoms corresponding to the elements

# Define the range of lattice constants and step size
lattice_constants = np.arange(5.00, 9.00, 0.10)

# Provide the path to the parent directory containing the {system} directory
base_path = "/home/yuqinglin/MuST/MST/Morse/Bridges2/"

# Data reading -------------------------------------------------------------------------------------------------------

# Calculate the number of atoms in the unit cell
number_of_atoms = sum(numbers) 
print("System:", system, variation)

# Initialize empty lists to store energy, volume, and pressure-volume data
E_data = []
V_cell = [] # The volume of the unit cell

# Initialize boolean for convergence check
success = True

# Loop over lattice constants and file names
for lattice_constant in lattice_constants:
    # File path for the k and o output file
    if variation:
        k_file = f"{base_path}{system}/{variation}/{lattice_constant:.2f}/k_n0000000_{system}"
        o_file = f"{base_path}{system}/{variation}/{lattice_constant:.2f}/o_n0000000_{system}"
    else:
        k_file = f"{base_path}{system}/{lattice_constant:.2f}/k_n0000000_{system}"
        o_file = f"{base_path}{system}/{lattice_constant:.2f}/o_n0000000_{system}"

    # Read E data from k file
    with open(k_file, 'r') as f:
        lines = f.readlines()

        # Find the line that starts with "** " (indicating the start of iteration data)
        for line in lines:
            if line.startswith('** '):
                break

        # Check if the lines list is not empty
        if lines:
            # Extract the E and 3PV value from the line after the last iteration
            last_iteration_line = lines[-1]
            if last_iteration_line.startswith('** '):
                E_value = float(last_iteration_line.split()[2])  # Only one iteration, E value is in the third column
            else:
                E_value = float(last_iteration_line.split()[1])  # Otherwise E value is in the second column
        else:
            # Handle the case when the lines list is empty
            print("No data found in file:", k_file)

        # Read the energy offset value
        for line in lines:
            if "Energy Offset" in line:
                E_offset = float(line.split(':')[1].strip())
                break

        # Calculate the final E value by adding the Energy Offset
        final_E_value = E_offset + E_value
        E_data.append(final_E_value)

    with open(o_file, 'r') as f:
        lines = f.readlines()

        # Check if the calculations converged
        convergence = any("SCF Convergence is reached !!!" in line for line in lines)
        if not convergence:
            print("Calculation did not converge in:", o_file)
            success = False

        vectors = []
        for i, line in enumerate(lines):
            if line.startswith("# Bravais Lattice Vector"):
                # Read the three lines following the vector header
                v1 = list(map(float, lines[i+1].strip().split()[1:]))
                v2 = list(map(float, lines[i+2].strip().split()[1:]))
                v3 = list(map(float, lines[i+3].strip().split()[1:]))
                vectors = [v1, v2, v3]
                break

        # Calculate the volume
        if len(vectors) == 3:
            v1, v2, v3 = vectors
            volume = np.abs(np.dot(np.cross(v1, v2), v3))
        else:
            print("Could not find all three Bravais lattice vectors in:", o_file)

        V_cell.append(volume)

if success:
    print("All calculations converged!")

V_data = [v / number_of_atoms for v in V_cell]  # Volume per atom

# File generation ----------------------------------------------------------------------------------------------------

# Write system, variation, elements, and numbers data to the input file
with open(input_file_path, "w") as input_file:
# Write system information
    input_file.write(
        f"# Name of the system\n"
        f"System: {system}\n\n"
        f"# Different parameters used in calculations of the same system\n"
        f"Variation: {variation}\n\n"
        f"# Elements in the system (separated by spaces)\n"
        f"Elements: {' '.join(elements)}\n\n"
        f"# Number of atoms in the unit cell corresponding to the elements (separated by spaces)\n"
        f"Numbers: {' '.join(map(str, numbers))}\n\n"
    )
    # Write the data with two columns: energy and volume per atom
    input_file.write(f"# Energy and volume per atom data (one data point per line)\n")
    for energy, volume_per_atom in zip(E_data, V_data):
        input_file.write(f"{energy} {volume_per_atom}\n")

print("Data saved to:", input_file_path)