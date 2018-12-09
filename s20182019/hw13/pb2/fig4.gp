set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig4.tex"

set xlabel '$t$ ($\mu$s)'
set ylabel '$I$ ($10^{-7}$ A)'

set xtics 10
set mxtics 2
set ytics 4
set mytics 2

plot "sol_500000.dat" u ($1*1e+6):($4*1e+7) w l notitle,\
"sol_500000_exact.dat" u ($1*1e+6):($2*1e+7) w l notitle

set output
system('latex fig4.tex && dvips fig4.dvi && ps2pdf fig4.ps')
