set terminal png
set output "frequency_pattern.png"
set xlabel "Pixel"
set ylabel "Pixel"
set autoscale fix
set link x
set link y
set size square

plot input_file binary array=(256,256) format='%uint32' with image notitle