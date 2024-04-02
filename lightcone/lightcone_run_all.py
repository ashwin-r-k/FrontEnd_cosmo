import sys
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline


#sys.path.insert(1,"../../Save/Run/")
#from values import *

print(os.getcwd())
location_home=os.getcwd() 


def edit_line_in_file(file_path, line_number, new_content):
    try:
        # Open the file in read mode
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Check if line number is within the range of lines in the file
        if 0 < line_number <= len(lines):
            # Edit the desired line
            lines[line_number - 1] = new_content + '\n'

            # Open the file in write mode to overwrite the original content
            with open(file_path, 'w') as file:
                file.writelines(lines)
            print("Line edited successfully.")
        else:
            print("Line number is out of range.")

    except FileNotFoundError:
        print("File not found. Please check the file path.")

def read_int_from_line(file_path, line_number):
    try:
        # Open the file in read mode
        with open(file_path, 'r') as file:
            # Read all lines from the file
            lines = file.readlines()

            # Check if line number is within the range of lines in the file
            if 0 < line_number <= len(lines):
                # Get the line content at the specified line number
                line_content = lines[line_number - 1].strip()
                
                # Attempt to convert the line content to an integer
                try:
                    int_value = int(line_content)
                    return int_value
                except ValueError:
                    print("The content of the line is not an integer.")
                    return None
            else:
                print("Line number is out of range.")
                return None

    except FileNotFoundError:
        print("File not found. Please check the file path.")
        return None

line_number_to_read = 5  # Line number to read the integer f


# Usage example
file_path = '../../Save/Run/inputs/input.lightcone'  # Path to the text file
line_number_to_edit = 4  # Line number to edit
edit_line_in_file(file_path, line_number_to_edit, str(1))
end=read_int_from_line(file_path, line_number_to_read)

#Run Once and then for remaing process later 
#First get the value of needed.
subprocess.run("./lightcone", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
file_needed="../../Save/Run/inputs/input.needed"


with open(file_needed, 'r') as file:
        lines = file.readlines()
        needed=int(lines[0])
print("needed ",needed)
print("end ",end)


curent_file = 1+needed  #Running for all the values of x_HI
while curent_file+needed<end:
     print("curent_file : " ,curent_file )
     edit_line_in_file(file_path, line_number_to_edit, str(curent_file))
     subprocess.run("./lightcone", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
     curent_file=curent_file+needed
print("Compleated !")



# Defining parametrized arctangent function
def f(x, k, w, x0, y0):
    return k * np.arctan(w*(x-x0)) + y0

def search_and_plot(directory, search_string):
    # Get a list of all files in the directory
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    
    # Filter files by the provided search string
    files = [f for f in files if search_string in f][::-1]
    
    if not files:
        print("No files found in the directory matching the search string.")
        return
    
    z=[]
    xH=[]
    for file_name in files:
        file_path = os.path.join(directory, file_name)
        # Process each file
        with open(file_path, 'r') as file:
            # Assuming the file contains two columns of data separated by whitespace
            data = [list(map(float, line.strip().split())) for line in file]
            # Extract x and y values
            x_values = [row[0] for row in data]
            y_values = [row[1] for row in data]

            z=z+x_values
            xH=xH+y_values

            # Plot the data
            plt.plot(x_values, y_values, label=file_name)
    # Add legend and labels
  
    plt.xlabel('X-axis label')
    plt.ylabel('Y-axis label')
    plt.title('Data from Files')
    plt.grid(True)

    # Defining data points
    xdata = np.array(z)
    ydata = np.array(xH )

    sorted_indices = np.argsort(xdata)
    # Sort the first array
    xdata = xdata[sorted_indices]
    # Rearrange the elements of the second array based on the sorting order of the first array
    ydata = ydata[sorted_indices]

    try:
        # Estimating initial parameter values
        p0 = [1, 1]
        p0.append(0.5*(max(xdata)+min(xdata)))
        p0.append(0.5*(max(ydata)+min(ydata)))
        # Fitting data
        output = curve_fit(f, xdata, ydata, p0=p0, full_output=True)
        # Extracting relevant information from output
        copt = output[0]
        res = output[2]['fvec']
        numeval = output[2]['nfev']
        msg = output[3]
        # Plotting data and fitted function
        xi = np.linspace(min(xdata), max(xdata), 100)
        yi = f(xi, *copt)
        plt.plot(xi,yi,label="ArcTan")
    except:
        print("Cannot Optimise using curvefit.")
        spl = UnivariateSpline(xdata,ydata,k=2)
        xs = np.linspace(min(xdata), max(xdata), 5*len(xdata))
        plt.plot(xs, spl(xs), 'g', lw=3,label="Smooth Fit")


    plt.legend()
    plt.title('Reionization History')
    output_file = "z_xHI.png"  # Specify the desired file name and format
    plt.savefig(os.path.join(directory, output_file))
    return xdata,ydata

directory_path = "../../Save/Run/lightcone/"  # Directory path provided as command line argument
search_string = "z_xHI"  # Search string provided as command line argument

z,xH=search_and_plot(directory_path, search_string)

# Convert lists to NumPy arrays
array1 = np.array(z)
array2 = np.array(xH)

# Stack arrays horizontally
stacked_array = np.column_stack((array1, array2))

# Output file path
output_file = "z_xHI.txt"

output_file = os.path.join(directory_path, output_file)
# Save stacked array to a text file
np.savetxt(output_file, stacked_array, delimiter='\t', fmt='%f')

print("Lists saved into output.txt file with two columns.")

