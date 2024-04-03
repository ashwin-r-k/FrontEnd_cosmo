#!/bin/bash

nion=$1

echo "nion value in bash $nion "

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

#clean the files
cd ../Save/Run/
rm -rv lightcone
rm -rv ionz_out_lc

cd "$LOC" || exit 1

cd reionz_lc || exit 1
./ionz_main_var_nion $nion || exit 1 &
Ionz_lc_PID=$!
wait $Ionz_lc_PID
echo "ionz lc Done"


cd "$LOC" || exit 1
cd lightcone || exit 1
python3 lightcone_run_all.py || exit 1
echo "Waiting for python lightcone run all"
wait
echo "lightcone Done"

# Moving the files
cd "$LOC" || exit 1
cd ../Save/Run/
#ls -la
echo $nion > ./lightcone/nion

# Creating directories
mkdir -p "./nion/$nion/"

#ls -la

# Moving folders
mv lightcone "./nion/$nion/"
mv ionz_out_lc "./nion/$nion/"

echo "Move compleated for nion $nion "

#change

