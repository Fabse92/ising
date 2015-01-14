set title "Magnetisierung pro Spin" 
set xlabel "Temperatur [Kelvin ? ]"
set ylabel "M/N"
set autoscale;
set yrange [-1:1]

set terminal png size 600,600
set output 'MagperSpin.png'


plot "MagperSpin" using 1:2 title 'M/N im Temperaturverlauf' with linespoints lt rgb "red"
