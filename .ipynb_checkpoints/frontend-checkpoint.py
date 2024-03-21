import subprocess
import os
print(os.getcwd())

location_home=os.getcwd()
location_nbody = "N-body"
location_fof = "FoF-Halo-finder"
location_regio_yuga = "ReionYuga"
location_py_env_act="fe-env/bin"

def compile_prog(location,compile_command,location_home):
# Compile the C++ program 
        os.chdir(os.path.join(location_home,location))
        print(os.getcwd())

        if (os.getlogin()=="ashwin"):
        	compile_command = "make"+" "+compile_command+" "+"FFTW=/usr/local"
        else:
        	compile_command = "make"+" "+compile_command+" "+"FFTW=/gpfs-scratch/m220590ph/softwares/fftw"
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

seed = 3141592
Nbin = 10
hh = 0.6704
omega_m = 0.3183
omega_l = 0.6817
spectral_index =0.9619
omega_baryon = 0.04902
sigma_8 = 0.8347
box=128
fraction_fill = 2
LL= 0.07
jcheck 
a_initial = 0.008

delta_a= 0.004
output_flag=0
pk_flag=0

redshift_values = [i for i in range(50,10,-1)]

Len_redshift = len(redshift_values)
redshift_values_str = ' '.join(str(x) for x in redshift_values)

lines = f"""{int(seed)} {int(Nbin)}
{hh} {omega_m} {omega_l} {spectral_index}
{omega_baryon} {sigma_8}
{box} {box} {box} {int(fraction_fill)} {LL}
{output_flag} {pk_flag}
{a_initial} {delta_a}
{Len_redshift}
{redshift_values_str}"""

os.chdir(os.path.join(location_home,location_nbody))
with open('input.nbody_comp', 'w') as file:  
        # Using the writelines function  
        file.writelines(lines)  
with open('input.nbody_comp', 'r') as file:  
        # Using the writelines function  
        lines = file.readlines()
        print(lines)


