# Phase-Stability-Tool-for-High-Entropy-Alloys
This tool has been created to understand the single phase stability of high entropy alloys, tool can also be used to create high entropy alloy model for molecular dynamics simulation in LAMMPS. I will try to develop this code with new features in future. The code language is created with Python 3 language.


⚛️ HEA Thermodynamic and Structural Calculation Tool

This Python-based tool is developed to calculate the thermodynamic and structural properties of High-Entropy Alloys (HEAs).
It takes alloy elements and their compositions as input from the user, then calculates various parameters based on the Miedema model and literature criteria.
The results are exported into an Excel file and a LAMMPS input file.

📂 Required Files

element_properties.xlsx
Contains the fundamental physical and chemical properties of elements:
Atomic number
Atomic weight
Density
Atomic radius
Melting point
Crystal lattice parameter
Vickers hardness
Young’s modulus
Thermal neutron absorption cross-section
Valence electron concentration (VEC)
Pauling electronegativity

miedema_matrix.xlsx
Contains mixing enthalpy (ΔH) values for binary systems.

⚙️ Installation

Install the required Python libraries:

pip install pandas numpy matplotlib openpyxl

▶️ Usage

When the program runs, the user is prompted step by step:

Element Selection
Enter the elements in the format:

24Cr,27Co,40Zr


(This matches the atomic number + symbol combination in the file.)

Composition Input
Enter the percentage ratios of the selected elements (total = 100).

Temperature Input (K)
Input the temperature value for Gibbs free energy calculations.

📊 Calculated Parameters

ΔS_mix → Mixing Entropy (J/mol·K)
ΔH_mix → Mixing Enthalpy (kJ/mol)
ΔG_mix → Gibbs Free Energy (J/mol)
Tm → Average Melting Temperature (K)
Ω (Omega) → Maximum Entropy Ratio
δ → Atomic Size Difference (%)
VEC → Valence Electron Concentration
χ̄ and Δχ → Mean electronegativity and difference

📑 Stability Criteria

ΔH_mix → between −10 and +5 kJ/mol
ΔS_mix ≥ 1.5 · R
Ω ≥ 1.1
δ ≤ 6.6%
Δχ% ≤ 8 → increases the probability of forming a single-phase solid solution
The program evaluates whether the alloy is likely to form a single-phase solid solution based on these criteria.

📤 Outputs

alloy_summary.xlsx
All calculated parameters are saved in an Excel table.
LAMMPS input file (in_Element1_Element2_...lmp)
A simple LAMMPS input file is generated for the selected elements.
Includes:
Simulation cell definition
Alloy atom definitions
Masses and initial parameters
Output: structure.lmpdat

🧑‍🔬 Example Run

Input:

Keys: 24Cr, 27Co, 40Zr
Percentage for Cr: 33.3
Percentage for Co: 33.3
Percentage for Zr: 33.3
Temperature (K): 1200


Output (summary):

ΔS_mix = 9.13 J/mol·K
ΔH_mix = −2.85 kJ/mol
ΔG_mix = −12,192 J/mol
Phase prediction = BCC + FCC
δ = 5.7 % (valid)
Ω = 2.14 (valid)

📌 Notes

Interactive user input is required.

os.startfile() is suitable for Windows; alternative methods may be needed on Linux/macOS.

The LAMMPS section is a simplified template; for real simulations, additional pair_style and potential files must be included.
