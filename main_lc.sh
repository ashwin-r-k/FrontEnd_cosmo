#!/bin/bash

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

mkdir -p ../Save/Run/inputs || exit 1
echo "Inputes dir created"

cd "r_z" || exit 1
echo "$(pwd)"
python3 r_z.py || exit 1

#./r_z || exit 1 &
#RZ_PID=$!  # Save the process ID of the background process

#wait $RZ_PID
#mv -f input.sampling "../sampling/input.sampling" || exit 1

echo "R_Z lightcone Done"

cd "$LOC" || exit 1
cd sampling || exit 1
echo "$(pwd)"
./random || exit 1 &
S_PID=$!  # Save the process ID of the background process
wait $S_PID
cd "$LOC" || exit 1
echo "Sampling Done"

#mv -f sampling/input.sampling reionz_lc/input.sampling || exit 1

cd reionz_lc || exit 1
./ionz_main || exit 1 &
Ionz_lc_PID=$!
wait $Ionz_lc_PID
echo "ionz lc Done"


cd "$LOC" || exit 1
cd lightcone || exit 1
./lightcone || exit 1 &
LC_PID=$!
wait $LC_PID
echo "lightcone Done"
