set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig0.tex"

set xlabel '$V_G$'
set ylabel '$n_{2D}$ (/cm$^2$)'
set ytics('$10^{9}$'1e9,'$10^{10}$'1e10,'$10^{11}$'1e11,'$10^{12}$'1e12,'$10^{13}$'1e13,'$10^{14}$'1e14)


set xtics 0.2
set mxtics 2
unset key
set logscale y 10
plot "integrated_density.dat" u 1:($2/10000) pt 6 notitle

set output
system('latex fig0.tex && dvips fig0.dvi && ps2pdf fig0.ps')
