set term epslatex 8 standalone size 6.0,3.0 font ",12"
set output "fig0.tex"

set multiplot layout 1,2 columnsfirst\
            margins screen .1, .95, .2, .95 spacing screen 0.08

set xlabel '$\log_{10}N^+$'
set ylabel '$\phi$'

set yrange [0:20]
set xrange [15:25]

set key at graph 0.6,0.9

plot 'sol.dat' u (log10($1)):2 pt 4 lc rgb "black" title '$\phi_{num}$',\
'sol.dat' u (log10($1)):3 pt 6 lc rgb "red" title '$\phi_{exact}$

set key at graph 0.9,0.9
set yrange [-20:0]
unset ylabel
set xlabel '$\log_{10}(-N^+)$'
plot 'sol.dat' u (log10(-$1)):2 pt 4 lc rgb "black" title '$\phi_{num}$',\
'sol.dat' u (log10(-$1)):3 pt 6 lc rgb "red" title '$\phi_{exact}$

unset multiplot
set output
system('latex fig0.tex && dvips fig0.dvi && ps2pdf fig0.ps')
