
import subprocess
import os
from values import *
import shlex, subprocess

print(os.getcwd())

location_home=os.getcwd()

location_r_z = "r_z"
location_sampling = "sampling"
location_lightcone="lightcone"
location_rion_lc="reionz_lc"


    
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


lines = f"""{int(seed)} {int(Nbin)}
{hh} {omega_m} {omega_l} {spectral_index}
{omega_baryon} {sigma_8}
{box} {box} {box} {int(fraction_fill)} {LL}
{output_flag} {pk_flag}
{a_initial} {delta_a}
{Len_redshift}
{redshift_values}"""

print(lines)


location="lightcone"
nion = Zarray(20,30.1,5)



def run_command( command ):
    subprocess.call(shlex.split(command))

os.chdir(location_home)

print("nion: ",nion)
for nion_i in nion:
    #run_command = "bash test.sh {nion}"
    command=f"bash lightcone_nion_var.sh {nion_i}"
    print(command)
    run_command(command)
    #then test code to handel all the running of files and then moving files and saving it in appropriate destinations.
print("DONEEEEEe")
