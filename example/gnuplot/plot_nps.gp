set terminal pdf
set output "nps.pdf"
set grid
#set xlabel 'time / µs'
set xlabel 'Frequency / Nyquist frequency'
#set yrange [0:2000]
#set xrange[0:120]
set logscale y
# 1E4, 5E4, 1E5, 5E5, 1E6, 5E6, 1E7, 5E7, 1E8, 5E8, 1E9

stats sprintf("%s/flux_%i_fourier.raw", input_dir, 1E4) binary format="%float" using 1:2 name "A"

plot sprintf("%s/flux_%i_fourier.raw", input_dir, 1E4) binary format="%float" using ($0/A_max_y):1 w l title "1E4", \
     sprintf("%s/flux_%i_fourier.raw", input_dir, 1E5) binary format="%float" using ($0/A_max_y):1 w l title "1E5", \
     sprintf("%s/flux_%i_fourier.raw", input_dir, 1E6) binary format="%float" using ($0/A_max_y):1 w l title "1E6", \
     sprintf("%s/flux_%i_fourier.raw", input_dir, 1E7) binary format="%float" using ($0/A_max_y):1 w l title "1E7", \
     sprintf("%s/flux_%i_fourier.raw", input_dir, 1E8) binary format="%float" using ($0/A_max_y):1 w l title "1E8", \
     sprintf("%s/flux_%i_fourier.raw", input_dir, 1E9) binary format="%float" using ($0/A_max_y):1 w l title "1E9"