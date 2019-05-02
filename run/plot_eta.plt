#system("grep Alpha out.512  > alpha.log")
#system("grep Terms out.512 > terms.log")
#system("grep torqz out.512 > torqz.log")

set y2tics
set mxtics
set mytics
set grid x
set grid y
set format y "%4.2E"
set format y2 "%4.2E"
set tics font "Verdana,12"
#set terminal enhanced font "Verdana,10"

#stats 'alpha.log' using 6 name "vel" nooutput
#stats 'alpha.log' using 5 name "pos" nooutput

#set yrange [-3*pos_stddev:3*pos_stddev]
#set yrange [-5.7:5.7]  # degrees
#set y2range [-0.05:0.05]
ttl1 = "pos"
ttl2 = "vel"

set ylabel '{/Symbol q}' font "Verdana,25"
set xlabel 'Time' font "Verdana,16"
set term x11 0 enhanced 
plot 'nlfsi_io.out' every 10::10 using 2:3 with lines axes x1y1 lw 2 title ttl1, \
 'nlfsi_io.out' every 10::10 using 2:4 with lines axes x1y2 lw 2 title ttl2

#set term x11 1 enhanced
##set yrange [-0.2:0.2]
#set y2range []
#unset y2tics
#set xlabel 'ISTEP'
#set ylabel 'Forces'
#plot 'fsi_terms.out' every 50::100 using 2:3 with lines axes x1y1 lw 2 title 'Fs', \
# 'fsi_terms.out' every 50::100 using 2:5 with lines axes x1y1 lw 2 title 'Ks'
#
#set term x11 2
#plot 'torqz.log' every 10::10 using 2:4 with lines axes x1y1 title 'torqz'

# 'etav.log' using (0.004*$2):3 axes x2y2 title 'vy'

pause 10

reread
