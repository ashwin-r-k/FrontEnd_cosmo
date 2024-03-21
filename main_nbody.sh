#!/bin/bash

# Define the CSV file path
csv_file="task_times.csv"

# Check if the CSV file exists, otherwise create it with a header
if [ ! -e "$csv_file" ]; then
    echo "Size,loading,N-body,FOF,ReionYuga,PLOT" > "$csv_file"
fi

# Function to calculate time difference
calculate_duration() {
    start_time=$1
    end_time=$2
    duration=$((end_time - start_time))
    echo "$duration"
}

task_name="loading"
# Get the current timestamp in seconds
start_timestamp=$(date +"%s")


# Get the current username
current_user=$(whoami)

# Check if the username is "ashwin"
if [ "$current_user" == "ashwin" ]; then

echo "Running command for user ashwin"
__conda_setup="$('/opt/anaconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/anaconda/etc/profile.d/conda.sh" ]; then
        . "/opt/anaconda/etc/profile.d/conda.sh"
    else
        export PATH="/opt/anaconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

LOC=/home/ashwin/HPC/FrontEnd_cosmo/
    # Add your command for the "ashwin" user here
    # Example: command_for_ashwin
else
echo "Running command for HPC"

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/gpfs-home/m220590ph/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/gpfs-home/m220590ph/miniconda/etc/profile.d/conda.sh" ]; then
        . "/gpfs-home/m220590ph/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/gpfs-home/m220590ph/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate frontend || exit 1


LOC=/gpfs-scratch/m220590ph/FrontEnd_cosmo/
    # Add your command for other users here
    # Example: command_for_others
fi


#source /gpfs-scratch/m220590ph/FrontEnd_cosmo/fe-env/bin/activate

#LOC=/gpfs-scratch/m220590ph/FrontEnd_cosmo/

cd "$LOC" || exit 1

bash clean  || exit 1

python3 frontend.py || exit 1


end_timestamp=$(date +"%s")

# Calculate the duration in seconds
loading_duration=$(calculate_duration "$start_timestamp" "$end_timestamp")




task_name="N-body"
# Get the current timestamp in seconds
start_timestamp=$(date +"%s")

cd N-body || exit 1
echo "$(pwd)"
./nbody_comp || exit 1 &
NBODY_PID=$!  # Save the process ID of the background process

wait $NBODY_PID
cd "$LOC" || exit 1
bash movefilesN-F  || exit 1

end_timestamp=$(date +"%s")

# Calculate the duration in seconds
n_duration=$(calculate_duration "$start_timestamp" "$end_timestamp")


echo "256,$loading_duration,$n_duration" >> "$csv_file"



