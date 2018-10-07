set term epslatex 8 standalone size 6.0,3.0 font ",12"
set output "fig0.tex"

set multiplot layout 1,2 columnsfirst\
            margins screen .15, .95, .2, .95 spacing screen 0.08

set xlabel '$\log_{10}N^+$'
set ylabel '$\log_{10}|\phi_{exact}-\phi_{num}|$'

set xrange [15:25]

unset key

plot 'sol.dat' u (log10($1)):(log10(abs($2-$3))) pt 4 lc rgb "black" title '$\phi_{num}$'

set key at graph 0.9,0.9
unset ylabel
set xlabel '$\log_{10}(-N^+)$'
plot 'sol.dat' u (log10(-$1)):(log10(abs(($2-$3)))) pt 4 lc rgb "black" title '$\phi_{num}$'

unset multiplot
set output
system('latex fig0.tex && dvips fig0.dvi && ps2pdf fig0.ps')
