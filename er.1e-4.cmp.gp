set log x
unset key
set xlabel "ea"
set log y
set grid
set ylabel "abs diff"
plot "er.1e-4.cmp" u 2:9
set title "abs (ptilde - ptilde2) @ er=1e-4,p=0.1,...,0.9"
set term png
set output "er.1e-4.cmp.png"
repl
