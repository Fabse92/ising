reset
set title "2D-Ising Model" 
set autoscale fix
set tics scale 0
set terminal png size 600,600
set output '2D_Ising_Matrix.png'
set palette defined (-1 "black", 1 "white")

unset key
unset colorbox

plot "data" matrix with image
