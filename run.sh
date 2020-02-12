# Script to run multiGWAS in a new terminal
# $1 must the filename of the config file
cfg=$1
xterm -e "multiGWAS-r.R $1;p" &
