set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig2.tex"

set xlabel 'Voltage (V)'
set ylabel '$J$ (A/m$^2$)'

set xtics 0.2
set mxtics 2
set yrange [0:1.42e+12]

plot "current_short.dat" u ($1):($2) w p pt 6 notitle,\

set output
system('latex fig2.tex && dvips fig2.dvi && ps2pdf fig2.ps')
