set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig4.tex"

set xlabel '$x$ (nm)'
set ylabel '$n$ (/m$^3$)'

set xtics 30
set mxtics 2

plot for [i=0:10] "ed".i."_x.dat" u ($1):($2) notitle

set output
system('latex fig4.tex && dvips fig4.dvi && ps2pdf fig4.ps')
