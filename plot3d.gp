set logscale z

set xlabel 'x'
set ylabel 'y'

set title 'L=100'

splot './data/L100.dat'

set terminal pdf
set output "plot3d.pdf"

splot './data/L100.dat'
