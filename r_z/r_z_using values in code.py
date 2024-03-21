
import subprocess
import os
print(os.getcwd())


def Zarray(start,stop,step):
    list_Z=[]
    i=start
    while i <= stop :
        list_Z.append(round(i,5))
        i=i+step
    return list_Z
 
seed = -100012
Nbin = 10
hh = 0.6704
omega_m = 0.3183
omega_l = 0.6817
spectral_index =0.9619
omega_baryon = 0.04902
sigma_8 = 0.8347
box=pow(2,8)
fraction_fill = 1
LL= 0.07
a_initial = 0.008
delta_a= 0.004
output_flag=0
pk_flag=1

redshift_values = Zarray(8,10,0.1)


Len_redshift = len(redshift_values)
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
