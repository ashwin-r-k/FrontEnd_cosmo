HELP For ALL Small things

Sign in GITHUB:
1. use gh auth login
2. login using https (nitc dosent block ssh)
3. then set the git config --global user.email 35029567+ashwin-r-k@users.noreply.github.com public email
4. git commit --amend --reset-author
5. pit gush / pull anthing 

Gnuplot
single line powerspectra plot
p_file is name of file
gnuplot -e "set terminal png size 800,600; set log; set output 'powerspectra_log.png'; plot '$p_file' w l "

