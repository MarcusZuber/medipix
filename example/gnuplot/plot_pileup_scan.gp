set terminal png
set output "pileup_scan.png"
set grid
set xlabel 'input photon flux / (1/(mm^2 s))'
set ylabel 'detected photon flux'


set logscale xy

plot sprintf("%s/pileup_scan.txt", input_dir) using ($1):($3) with lines title 'SPM', \
    sprintf("%s/pileup_scan.txt", input_dir) using ($1):($1) with lines title 'Ideal'
