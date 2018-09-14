set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig0.tex"

set xlabel '$x$'
set ylabel '$\phi$'

set ytics 0.2
set xtics 5
unset key
plot "vec.dat" pt 6 notitle

set output
system('latex fig0.tex && dvips fig0.dvi && ps2pdf fig0.ps')
