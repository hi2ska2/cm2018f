set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig3.tex"

set xlabel '$t$ (s)'
set ylabel '$I$ ($10^{-9}$ A)'

set xtics 0.02
set mxtics 2
set ytics 1
set mytics 2

plot "sol_100.dat" u ($1):($4*1e+9) w l notitle,\
"sol_100_exact.dat" u ($1):($2*1e+9) w l notitle

set output
system('latex fig3.tex && dvips fig3.dvi && ps2pdf fig3.ps')
