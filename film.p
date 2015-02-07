reset

set title "2D-Ising Model" 
set autoscale fix
set tics scale 0
set terminal png size 600,600
set palette defined (-1 "black", 1 "white")
unset key
unset colorbox


do for [i=1:25] {
set output 'png/plot'.i.'.png'
plot 'film/data'.i matrix with image
}
