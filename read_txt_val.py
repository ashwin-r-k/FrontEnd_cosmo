import sys
import os

# Define variables to store values
seed = Nbin = hh = omega_m = omega_l = spectral_index = None
omega_baryon = sigma_8 = box = fraction_fill = LL = None
output_flag = pk_flag = a_initial = delta_a = None
Len_redshift = None
redshift_values = []

# Check if the filename is provided as an argument
#if len(sys.argv) < 2:
#    print("Please provide the filename as an argument.")
#    sys.exit(1)

# Get the filename from command line argument
print(os.getcwd())
location_home=os.getcwd()

file_path = '../../Save/Run/input.nbody_comp'

# Read data from the file
with open(file_path, 'r') as file:
    data = file.read().splitlines()

# Assign values to variables
seed, Nbin = map(int, data[0].split())
hh, omega_m, omega_l, spectral_index = map(float, data[1].split())
omega_baryon, sigma_8 = map(float, data[2].split())
box, _, _, fraction_fill, LL = map(float, data[3].split())
output_flag, pk_flag = map(int, data[4].split())
a_initial, delta_a = map(float, data[5].split())
Len_redshift = int(data[6])
redshift_values = data[7]

# Convert redshift values from string to floats
redshift_values = list(map(float, redshift_values.split()))

# Print the values to verify
print("seed:", seed)
print("Nbin:", Nbin)
print("hh:", hh)
print("omega_m:", omega_m)
print("omega_l:", omega_l)
print("spectral_index:", spectral_index)
print("omega_baryon:", omega_baryon)
print("sigma_8:", sigma_8)
print("box:", box)
print("fraction_fill:", fraction_fill)
print("LL:", LL)
print("output_flag:", output_flag)
print("pk_flag:", pk_flag)
print("a_initial:", a_initial)
print("delta_a:", delta_a)
print("Len_redshift:", Len_redshift)
print("redshift_values:", redshift_values)
