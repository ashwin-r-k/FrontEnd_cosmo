import sys
import os
import subprocess

sys.path.insert(1,"../../Save/Run/")
#import imp
#module = imp.load_source("read_txt_val", '../../Save/Run/input.nbody_comp')
from values import *

print(os.getcwd())

location_home=os.getcwd() 

redshift_values = redshift_values[::-1]



redshift_values_str = ' '.join(str(x) for x in redshift_values)

omegabh2=0.022

lines_r_z = f"""{omega_m} {omegabh2} {hh}
{Len_redshift}
{redshift_values_str}"""

print(lines_r_z)
with open('../../Save/Run/inputs/input.r_z', 'w') as file:  
        # Using the writelines function  
        file.writelines(lines_r_z)
 
run_process = subprocess.run("ls", stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        # Get the output and error messages from the program 
output = run_process.stdout.decode() 
print(output)
# using os.system to run the command and redirect the output to a file
run_process = subprocess.run("./r_z", stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        # Get the output and error messages from the program 
output = run_process.stdout.decode() 

print(output)

"""
with open('output', 'w') as file:  
        # Using the writelines function  
        file.writelines(output)
"""

lines = f"""{int(seed)} {int(Nbin)}
{hh} {omega_m} {omega_l} {spectral_index}
{omega_baryon} {sigma_8}
{box} {box} {box} {int(fraction_fill)} {LL}
{output_flag} {pk_flag}
{a_initial} {delta_a}
"""+output

print(lines)

with open('../../Save/Run/inputs/input.sampling', 'w') as file:  
        # Using the writelines function  
        file.writelines(lines)  
