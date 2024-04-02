#!/bin/bash

# Directory containing the files
directory="/home/ashwin/HPC/Save/Run/lightcone/"

# File name pattern to search for
file_name_pattern="z_xHI_*"

# Output file
output_file="concatenated_files.txt"

# Find files with the specified name pattern in the directory
files=$(find "$directory" -name "$file_name_pattern")

# Concatenate all matching files into the output file
cat $files > "$output_file"

echo "Files concatenated into $output_file"
home=`pwd`
gnuplot -e "set terminal png size 800,600; set output 'plot.png'; plot '$loc/$output_file' "
