set terminal png
set output "pixel_signal.png"
set grid
set xlabel 'time / µs'
set ylabel 'effective charge / keV'
#set yrange [0:*]
#set logscale x
plot sprintf("%s/pixel_signal_0_0.raw", input_dir) binary format="%float" using ($0/100.):1 w l title "pixel 0, 0", \
        sprintf("%s/pixel_signal_0_1.raw", input_dir) binary format="%float" using ($0/100.):1 w l title "pixel 0, 1", \
        sprintf("%s/pixel_signal_1_0.raw", input_dir) binary format="%float" using ($0/100.):1 w l title "pixel 1, 0", \
        sprintf("%s/pixel_signal_1_1.raw", input_dir) binary format="%float" using ($0/100.):1 w l title "pixel 1, 1", \
        8 w l title "threshold"
