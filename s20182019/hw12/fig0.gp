set term epslatex 8 standalone size 6.0,3.0 font ",12"
set output "fig1.tex"

set multiplot

set size 0.5,1.0
set origin 0.0,0.0

set xlabel '$x$ (nm)'
set ylabel '$n$'

set rmargin 2

set xtics 1

set key samplen 2
set key at screen 0.51,0.87
set logscale y 10

set title 'Long structure'
plot 'EDlong_62_init.dat' u 1:($2+0.001) title 'Poisson',\
'ED_long_62_after.dat' u 1:($2+0.001) title 'Coupled'

set size 0.5,1.0
set origin 0.5,0.0
set title 'Short structure'
set key at screen 1.0,0.87
plot 'EDshort_26_init.dat' u 1:($2+0.001) title 'Poisson',\
'ED_short_26_init.dat' u 1:($2+0.001) title 'coupled'

unset multiplot
set output
system('latex fig1.tex && dvips fig1.dvi && ps2pdf fig1.ps')
