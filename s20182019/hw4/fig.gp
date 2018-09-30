set term epslatex 8 standalone size 6.0,3.0 font ",12"
set output "fig0.tex"

set xlabel '$x$ (nm)'
set ylabel '$\Delta\phi$'

set rmargin 13

set xtics 1

set key at screen 1.0,0.87


plot for [i=0:9] 'data'.i.'.dat' u 1:($3) title '0.'.i.' V',\
'data10.dat' u 1:($3) title '1.0 V'

set output
system('latex fig0.tex && dvips fig0.dvi && ps2pdf fig0.ps')
