set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig1.tex"

set xlabel '$y$ (nm)'
set ylabel '$x$ (nm)'
set zlabel '\small $\phi$ (V)' offset 1.0,0.0

set xtics 2
set mxtics 2
set ytics 40
set mytics 2

set ztics 0.2

set view 65,100

set key at screen 1.05,0.98


splot for [i=0:3] "phi".i.".dat" u 1:2:($3) title '\small $V_g=$0.'.i.' V'
#"phi10.dat" u 1:2:($3) title '\small $V_g$=1 V'

set output
system('latex fig1.tex && dvips fig1.dvi && ps2pdf fig1.ps')
