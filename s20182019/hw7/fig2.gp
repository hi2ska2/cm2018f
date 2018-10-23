set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig1.tex"

set xlabel '$E_F$ (eV)'
set ylabel '$n_{2D}$ (/cm$^2$)'
set ytics('$10^{9}$'1e9,'$10^{10}$'1e10,'$10^{11}$'1e11,'$10^{12}$'1e12,'$10^{13}$'1e13,'$10^{14}$'1e14)


set xtics 0.05
set mxtics 2
set xrange [-0.12:0.12]
unset key
set logscale y 10
set pointsize 1.5
plot "integrated_density.dat" u 1:($2) pt 6 notitle

set output
system('latex fig1.tex && dvips fig1.dvi && ps2pdf fig1.ps')
