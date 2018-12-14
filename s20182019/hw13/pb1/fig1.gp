set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig1.tex"

set xlabel 'Voltage (V)'
set ylabel '$J$ (A/m$^2$)'

set xtics 0.2
set mxtics 2

plot "current_long.dat" u ($1):($2) w p pt 6 notitle,\

set output
system('latex fig1.tex && dvips fig1.dvi && ps2pdf fig1.ps')
