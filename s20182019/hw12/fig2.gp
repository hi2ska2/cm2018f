set term epslatex 8 standalone size 6.0,3.0 font ",12"
set output "fig2.tex"

set multiplot

set size 0.5,1.0
set origin 0.0,0.0

set xlabel '$x$ (nm)'
set ylabel '$n$' offset 1.5,0.0

set lmargin 6
set rmargin 6

set xtics 1

set key samplen 2
set key at screen 0.51,0.87
set logscale y 10

l1 = 1
l2 = 1

set xtics 200
set xrange [0:600]

set ytics('$10^{24}$'1e+24,'$10^{23}$'1e+23,'$10^{22}$'1e+22,'$10^{21}$'1e+21)

set title 'Long structure'
plot 'EDlong_601_init.dat' u ($1*l1):($2+0.001) title '\small P',\
'ED_long_601_after.dat' u ($1*l1):($2+0.001) title '\small C'

set xrange [0:120]
set ytics('$10^{24}$'1e+24,'$10^{23}$'1e+23,'$10^{25}$'1e+25,'$10^{26}$'1e+26,'$10^{27}$'1e+27)
set xtics 30
set size 0.5,1.0
set origin 0.5,0.0
set title 'Short structure'
set key at screen 1.01,0.87
plot 'EDshort_121_init.dat' u ($1*l2):($2+0.001) title '\small P',\
'ED_short_121_after.dat' u ($1*l2):($2+0.001) title '\small C'

unset multiplot
set output
system('latex fig2.tex && dvips fig2.dvi && ps2pdf fig2.ps')
