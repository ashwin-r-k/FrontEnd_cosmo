
import subprocess
import os
from values import *

print(os.getcwd())

location_home=os.getcwd()
location_nbody = "N-body"
location_fof = "FoF-Halo-finder"
location_regio_yuga = "ReionYuga"

location_r_z = "r_z"
location_sampling = "sampling"
location_lightcone="lightcone"
location_rion_lc="reionz_lc"

location_py_env_act="fe-env/bin"


    
def compile_prog(location,compile_command,location_home):
# Compile the C++ program 
        os.chdir(os.path.join(location_home,location))
        print(os.getcwd())
        home = os.path.expanduser("~")
        compile_command = "make"+" "+compile_command+" "+f"FFTW={home}/softwares/fftw"
        print(compile_command)
        os.system(compile_command)
        '''
        #compile_process = subprocess.Popen([location, "make",compile_command])
        #compile_process = subprocess.run(compile_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         
        # Check if the compilation was successful 
        if compile_process.returncode == 0:
            print("Compilation successful.") 
        else: 
            # Print the compilation error messages 
            print("Compilation failed.") 
            print(compile_process.stderr.decode())
        '''
        os.chdir(location_home)

def run_prog(location,run_command,location_home):
            # Run the compiled program 
        os.chdir(os.path.join(location_home,location))

        run_process = subprocess.run(run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
         
        # Get the output and error messages from the program 
        output = run_process.stdout.decode() 
        error = run_process.stderr.decode() 
         
        # Print the output and error messages 
        print("Output:") 
        print(output) 
         
        print("Error (if any):") 
        print(error) 
        os.chdir(location_home)


prog = "nbody_comp"
compile_prog(location_nbody,prog,location_home)

prog = "fof_main"
compile_prog(location_fof,prog,location_home)

prog = "ionz_main"
compile_prog(location_regio_yuga,prog,location_home)

prog = "r_z"
compile_prog(location_r_z,prog,location_home)

prog = "random"
compile_prog(location_sampling,prog,location_home)

prog = "lightcone"
compile_prog(location_lightcone,prog,location_home)

prog = "ionz_main"
compile_prog(location_rion_lc,prog,location_home)

redshift_values_str = ' '.join(str(x) for x in redshift_values)

lines = f"""{int(seed)} {int(Nbin)}
{hh} {omega_m} {omega_l} {spectral_index}
{omega_baryon} {sigma_8}
{box} {box} {box} {int(fraction_fill)} {LL}
{output_flag} {pk_flag}
{a_initial} {delta_a}
{Len_redshift}
{redshift_values_str}"""

print(lines)

os.chdir(os.path.join(location_home,location_nbody))
with open('input.nbody_comp', 'w') as file:  
        # Using the writelines function  
        file.writelines(lines)  
        
