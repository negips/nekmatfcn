system("grep Qk logfile > wts")

set logscale y

set format y "%4.2E"
set tics font "Verdana,12"
#set terminal enhanced font "Verdana,10"

# set ylabel '{/Symbol q}' font "Verdana,25"
# set xlabel 'Time' font "Verdana,16"
set term x11 0 enhanced
plot 'wts' using 3:4 with points axes x1y1 ps 1 pt 13 

pause 10

reread
