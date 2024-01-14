#!/bin/bash

#!/bin/bash

# Get the current username
current_user=$(whoami)

# Check if the username is "ashwin"
if [ "$current_user" == "ashwin" ]; then
echo "Running command for user ashwin"
source ./fe-env/bin/activate
LOC=$(pwd)
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

python3 frontend.py || exit 1

cd N-body || exit 1
echo "$(pwd)"
./nbody_comp || exit 1 &
NBODY_PID=$!  # Save the process ID of the background process

wait $NBODY_PID

cd "$LOC" || exit 1
bash movefilesN-F  || exit 1

cd FoF-Halo-finder || exit 1
echo "$(pwd)"
./fof_main || exit 1 &
FOF_PID=$!  # Save the process ID of the background process
wait "$FOF_PID"

cd "$LOC" || exit 1
bash movefilesF-R  || exit 1

# Wait for the background process to finish


cd ReionYuga || exit 1
echo "$(pwd)"
./ionz_main || exit 1 &

RION_PID=$!  # Save the process ID of the background process
wait "$RION_PID"
cd "$LOC" || exit 1
echo "Exicution Completed"

bash backupOutput
