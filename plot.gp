set terminal pdf
set output "plot.pdf"

set logscale y
set xlabel 't'

plot for [i=0:5] './data/L100p'.i.'.dat' title '<~{/Symbol e}{.6\~} _'.i.'~{/Symbol e}{.6\~} _'.i.'>_c'
