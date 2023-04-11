set terminal pdf
set output "frequency_pattern.pdf"
set xlabel "Pixel"
set ylabel "Pixel"
set autoscale fix
set link x
set link y
set size square

plot input_file binary array=(256,256) format='%uint32' with image notitle