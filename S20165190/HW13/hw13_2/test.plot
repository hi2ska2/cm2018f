set term post
set output 'test.eps'
set xlab "Time(hour)"
set ylab "Distance(meter)"
set title "Velocity of My Car"
set nokey
plot "data.txt" with linespoints


