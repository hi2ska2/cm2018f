set term epslatex 8 standalone size 8.0,3.0 font ",12"
set output "fig0.tex"

set multiplot

set size 0.5,1.0
set origin 0.0,0.0

set xlabel '$z$ (nm)'
set ylabel '$n$ $(10^{18}/$cm$^3)$'

set rmargin 10

set xtics 1

set key samplen 2
set key at screen 0.51,0.87


plot for [i=0:10] 'density_'.i.'.dat' u 1:($2*1e-18) title sprintf("%.2f V",(i-10)*0.01) 

set size 0.5,1.0
set origin 0.5,0.0
set key at screen 1.0,0.87
plot for [i=10:20] 'density_'.i.'.dat' u 1:($2*1e-18) title sprintf("%.2f V",(i-10)*0.01)

unset multiplot
set output
system('latex fig0.tex && dvips fig0.dvi && ps2pdf fig0.ps')
