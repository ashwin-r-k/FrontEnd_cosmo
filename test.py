import subprocess
import os
print(os.getcwd())

location_home=os.getcwd()
location_nbody = "N-body"
location_fof = "FoF-Halo-finder"
location_regio_yuga = "ReionYuga"
    
def compile_prog(location,compile_command,location_home):
# Compile the C++ program 
        os.chdir(os.path.join(location_home,location))
        print(os.getcwd())

        compile_command = "make"+" "+compile_command
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

def run_prog(run_command,location_home):
            # Run the compiled program 

        run_process = subprocess.run(run_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
         
        # Get the output and error messages from the program 
        output = run_process.stdout.decode() 
        error = run_process.stderr.decode() 
        print(os.getcwd())

        # Print the output and error messages 
        print("Output:") 
        print(output) 
         
        print("Error (if any):") 
        print(error) 
        os.chdir(location_home)

prog="/home/ashwin/HPC/FrontEnd_cosmo/fe-env/bin/activate"
run_prog(prog,location_home)