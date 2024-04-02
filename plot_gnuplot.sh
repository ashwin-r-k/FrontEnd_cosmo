#!/bin/bash

# Directory containing the files
directory="/home/ashwin/HPC/Save/Run/lightcone/"

# File name pattern to search for
file_name_pattern="z_xHI_*"

# Output combined plot file
output_plot="combined_plot.png"

# Find files with the specified name pattern in the directory
files=$(find "$directory" -name "$file_name_pattern")

# Gnuplot script file
gnuplot_script="plot_script.gp"

# Create a Gnuplot script to plot all files
echo "set term png" > "$gnuplot_script"
echo "set output '$output_plot'" >> "$gnuplot_script"
echo "plot \\" >> "$gnuplot_script"

# Append plot commands for each file
first_file=true
for file in $files; do
    if [ "$first_file" = true ]; then
        echo "'$file' with lines title '$file'" >> "$gnuplot_script"
        first_file=false
    else
        echo ", \\" >> "$gnuplot_script"
        echo "'$file' with lines title '$file'" >> "$gnuplot_script"
    fi
done

# Run Gnuplot with the generated script
gnuplot "$gnuplot_script"

# Clean up the temporary Gnuplot script
rm "$gnuplot_script"

echo "Combined plot created: $output_plot"
