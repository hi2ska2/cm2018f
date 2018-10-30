set term epslatex 8 standalone size 6.0,3.0 font ",12"
set output "fig0.tex"

set multiplot

set size 0.5,1.0
set origin 0.0,0.0

set xlabel '$z$ (nm)'
set ylabel '$n$ $(10^{20}$cm$^3)$'

set rmargin 2

set xtics 1

set key samplen 2
set key at screen 0.51,0.87
set logscale y 10

set title 'Semi-classical Poisson solver'
plot for [i=0:20] 'ed_cl_'.i.'.dat' u 1:($2*1e-26+0.01) notitle

set size 0.5,1.0
set origin 0.5,0.0
set title 'Schrodinger-Poisson solver'
set key at screen 1.0,0.87
plot for [i=0:20] 'ed_cor_'.i.'.dat' u 1:($2*1e-26+0.01) notitle

unset multiplot
set output
system('latex fig0.tex && dvips fig0.dvi && ps2pdf fig0.ps')
