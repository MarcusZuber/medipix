set terminal png
set output "threshold_scan.png"
set grid
set xlabel 'energy threshold / keV'
set ylabel 'differential counts'
set yrange [0:*]
d(y) = ($0 == 0) ? (y1 = y, 1/0) : (y2 = y1, y1 = y, y2-y1)

plot input_file_spm using ($1):(d($2)) with lines title 'SPM', input_file_csm using ($1):(d($2)) with lines title 'CSM'