set term epslatex 8 standalone size 4.0,3.0 font ",12"
set output "fig1.tex"

set xlabel '$V_G$ (V)'
set ylabel '$n_{2D}$ (/cm$^2$)'
set ytics('$10^{9}$'1e9,'$10^{10}$'1e10,'$10^{11}$'1e11,'$10^{12}$'1e12,'$10^{13}$'1e13,'$10^{14}$'1e14)


set key left bottom
set xtics 0.2
set mxtics 2
set xrange [0:1]
set logscale y 10
set pointsize 1.5
plot "integ_density_cl.dat" u 1:($2*1e-6) w lp pt 6 lc rgb "light-red" title 'Poisson',\
"integ_density_qm.dat" u 1:($2*1e-6) w lp pt 4 lc rgb "blue" title 'Schrodinger-Poisson'

set output
system('latex fig1.tex && dvips fig1.dvi && ps2pdf fig1.ps')
