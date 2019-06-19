a = 0.6*(1-sqrt(1-0.670))
d=0.46
f(x) = (0.6-a*d*d/((d+2*k*x)*(d+2*k*x)))/0.6
fit f(x) 'area_averaged_axial_mean_velocity_TI_15.txt' using 1:6 via k
# plotting
set terminal postscript eps font 24
set out 'k_fit_ti_15_tsr5.eps'
set autoscale
unset log
unset label
unset pm3d
set xtic auto
set ytic auto
unset grid
# set title 'Normalized velocity recover in the wake'
set xlabel  "Normalized axial distance, X/D"
set xrange [*:12]
# r0 initial pulse
set yrange [0.5:1]
set ylabel "Normalized mean axial velocity, ~U{0.8-} / U{/Symbol \245}"
set style line 1 lt 1 lc rgb "black" lw 4 pt 1 ps 2
set style line 2 lt 2 lc rgb "black" lw 4  pt 3 ps 2
set style line 3 lt 3 lc rgb "black" lw 4  pt 5 ps 2
set style line 4 lt 4 lc rgb "black" lw 4  pt 7 ps 2
set style line 5 lt 5 lc rgb "black" lw 4
set style line 6 lt 6 lc rgb "brown" lw 4
k_value = sprintf("k = %.3f", k)
set label 1 at 1, 0.85  k_value 
set key at graph 0.9, 0.3
set key spacing 1
plot 'area_averaged_axial_mean_velocity_TI_15.txt' using 1:6 ls 1 with points title 'RANS', f(x) lw 3 title "best fitted"
