set terminal pdf
set output "plot.pdf"

set logscale y
set xlabel 't'
set title 'L=100'

plot for [i=0:5] './data/L100p'.i.'.dat' title '<~{/Symbol e}{.6\~} _'.i.'~{/Symbol e}{.6\~} _'.i.'>_c'
