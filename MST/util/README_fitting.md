# Energy-Volume and Pressure-Volume Analysis with Morse and Birch-Murnaghan fitting

This script is designed to analyze energy-volume and pressure-volume curves obtained from calculations and perform Morse and Birch-Murnaghan fittings. It reads data from an input file containing energy and volume per atom data and fits the data to derive parameters such as equilibrium volume, bulk modulus, and Debye temperature. 

The user only needs to input the paths in the section marked with "USER INPUT NEEDED HERE!" The rest of the script handles data reading, fitting, plotting, and analysis. The results are displayed on the console and saved to the output files specified.

## Prerequisites

Make sure you have the following Python libraries installed:
- `numpy`
- `scipy`
- `matplotlib`
- `chempy`

## Usage

You only need to make changes to the section labeled "USER INPUT NEEDED HERE!" in order to use the code.

### User Input Needed

Specify the paths to the input file, output plot, and output values. For example:
```python
input_file_path = "/path/to/AgCu3_data.txt"
output_plot_path = "/path/to/AgCu3_fitting_plot.png"
output_file_path = "/path/to/AgCu3_fitting_values.txt"
```
Alternatively, uncomment the lines to provide them interactively during runtime.

<span style="color:green">**Note:**</span> The input file should follow a specific format, which is specified in the section **Input File Format** at the end of this README.

### Data Reading

The script reads the system-related data and energy-volume data from the input file.

### Energy-Volume and Pressure-Volume Analysis

The script reads the energy and volume per atom data from the input file and performs the fitings:

Fit the data to the Morse potential to extract optimal parameters.
The Morse function is defined as:
$$
E_{\text{morse}}(V) = A - 2 D \exp\left(-\lambda \left(\frac{3}{4 \pi}\right)^{\frac{1}{3}} \left(V^{\frac{1}{3}} - V_0^{\frac{1}{3}}\right)\right) + D \exp\left(-2 \lambda \left(\frac{3}{4 \pi}\right)^{\frac{1}{3}} \left(V^{\frac{1}{3}} - V_0^{\frac{1}{3}}\right)\right)
$$
- $A$: Additive constant
- $D$: Dissociation energy
- $\lambda$: Constant determining the curvature near the minimum energy
- $V$: Volume of the system
- $V_0$: Equilibrium volume

Fit the data to the Birch-Murnaghan equation of state to extract parameters such as equilibrium energy, equilibrium volume, bulk modulus, and its pressure derivative.

The Birch-Murnaghan function is defined as:
$$
E_{\text{bm}}(V) = E_0 + \frac{9 V_0 B_0}{16} \left(\left(\frac{V}{V_0}\right)^{\frac{2}{3}} - 1\right)^3 B_0^{\prime} + \left(\left(\frac{V}{V_0}\right)^{\frac{2}{3}} - 1\right)^2 \left(6 - 4 \left(\frac{V}{V_0}\right)^{\frac{2}{3}}\right)
$$
- $E_0$: Equilibrium energy
- $V$: Volume of the system
- $V_0$: Equilibrium volume
- $B_0$: Bulk modulus
- $B_0^{\prime}$: Pressure derivative of the bulk modulus

### Pressure Calculations

The pressure `P_morse` is calculated as the negative derivative of the Morse potential energy function `E_morse` with respect to volume `V`. The function is defined as follows:

$$
P_{\text{morse}}(V) = - D \lambda \left(\frac{9 \pi V^2}{2}\right)^{-\frac{1}{3}} \exp\left(-\lambda \left(\frac{3}{4\pi}\right)^{\frac{1}{3}} \left(V^{\frac{1}{3}} - V_0^{\frac{1}{3}}\right)\right) \left[1 - \exp\left(-\lambda \left(\frac{3}{4\pi}\right)^{\frac{1}{3}} \left(V^{\frac{1}{3}} - V_0^{\frac{1}{3}}\right)\right)\right]
$$

The parameters `V0`, `D` and `lambda` are the fitted Morse parameters obtained earlier. 

The pressure `P_bm` is calculated using the third-order Birch–Murnaghan isothermal equation of state and the fitted parameters. The function is defined as follows:

$$
P_{\text{bm}}(V) = \frac{3B_{0}}{2}\left[\left(\frac{V_{0}}{V}\right)^{\frac{7}{3}} - \left(\frac{V_{0}}{V}\right)^{\frac{5}{3}}\right] \left\{1 + \frac{3}{4}\left(B_{0}^{\prime} - 4\right)\left[\left(\frac{V_{0}}{V}\right)^{\frac{2}{3}} - 1\right]\right\}
$$

Here, `V0`, `B0` and `B0_prime` are the fitted Birch-Murnaghan parameters. 

The pressure values are computed for each volume value `V` in the dataset and converted from $\text{Ry}/\text{(a.u.)}^3$ to $\text{GPa}$.

### Bulk Moduli Calculations

The bulk modulus `B0` for the Morse potential is calculated at the equilibrium volume `V0_morse_fit`. The formula used to compute `B0` is based on the Morse parameters `D` and `lambda`, and the equilibrium volume `V0_morse_fit`. The formula is as follows:

$$
B_0 = \frac{D \lambda^2}{\left(162 \pi^2 V_0\right)^{\frac{1}{3}}}
$$

The bulk modulus `B0_bm` for the Birch-Murnaghan equation of state is directly obtained from the fitted parameter `B0_fit`.

The calculated bulk moduli are then converted to $\text{kbar}$ and $\text{GPa}$.

### Plotting

The script generates a plot showing the energy-volume and pressure-volume curves for the system. The plot contains:
- Energy-Volume data points.
- Fitted Morse function curve.
- Fitted Birch-Murnaghan function curve.
- Pressure values (GPa) calculated using the fitted curves.
- The points corresponding to the equilibrium volume marked for each fit.

The plot is saved to an image file specified by `output_plot_path`.

### Debye Temperatures Calculations

The Debye temperature $\theta_D$ is calculated with:
$$
\theta_D = 41.63 \times \left(\frac{3V_0}{4\pi}\right)^{\frac{1}{6}} B^{\frac{1}{2}} M^{-\frac{1}{2}}
$$
where $V_0$ is the volume per atom at the minimum energy configuration, $B$ is the bulk modulus at the minimum energy configuration, and $M$ is the average atomic mass per atom.

Using the calculated Morse and Birch-Murnaghan parameters, we can determine the Debye temperature for each fit.

### Output

The script saves the energy, volume, fitted parameters, and calculation results to an output file specified by `output_file_path`.

## Input File Format

In the input file, the parameters and data points should be specified as follows:

- Lines starting with "#" are treated as comments and ignored.
- Regarding the system data, each non-empty line should contain a parameter key and its value separated by a colon. For example, `System: AgCu3`
- The following parameter keys are expected in the input file:
    - `System`: The name of the system being analyzed (e.g., `AgCu3`).
    - `Variation`: (Optional) different parameters used in calculations of the same system (e.g., `MFT`).
    - `Elements`: Elements in the system, separated by spaces. (e.g., `Ag Cu`).
    - `Numbers`: The number of atoms in the unit cell corresponding to each element, listed in the same order and separated by spaces (e.g., `1 3`).
- After the parameter section, the data section should contain two columns. Each line contains two values separated by a space: the energy value and the volume per atom value. These data points are used for the Morse and Birch-Murnaghan fitting.

Below is an example of the input file:

```plaintext
# Input file example for the Morse and Birch-Murnaghan fitting
# Comments start with '#'

# Name of the system
System: AgCu3
# Different parameters used in calculations of the same system
Variation: MFT
# Elements in the system (separated by spaces)
Elements: Ag Cu
# Number of atoms in the unit cell corresponding to the elements (separated by spaces)
Numbers: 1 3

# Energy and volume per atom data (one data point per line)
-5060.20357 54.0
-5060.251018 56.745249999999984
-5060.290581 59.58200000000001
-5060.323369 62.51174999999999
-5060.350346 65.53600000000002
-5060.372286 68.65625
-5060.389891 71.874
-5060.403795 75.19075000000001
-5060.414537 78.60799999999999
-5060.422594 82.12725000000002
-5060.428408 85.75
-5060.432263 89.47774999999999
-5060.434445 93.31200000000001
-5060.435222 97.25425
-5060.434862 101.30600000000001
-5060.433505 105.46875
-5060.431338 109.74399999999999
-5060.428473 114.13325000000002
-5060.424968 118.63799999999999
-5060.421313 123.25975000000001
-5060.416902 128.0
-5060.412411 132.86024999999998
-5060.40767 137.84199999999998
-5060.402788 142.94675000000004
-5060.397797 148.17600000000002
-5060.39277 153.53125
-5060.387686 159.01399999999998
-5060.382582 164.62574999999995
-5060.377485 170.36800000000005
-5060.372422 176.24225
```
