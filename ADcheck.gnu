set term pdf enhanced font 'Helvetica,14'
set output "spec_AD.pdf" 
set size 1.0,1.0
set origin 0,0
set multiplot
set size 1.0, 0.5
set origin 0.0, 0.52
set format y '10^{%T}'
set xrange [1:2048]
set yrange [1e-48:1e-40]
set ytics (1e-46,1e-44,1e-42,1e-40)
unset xtics
set logscale y
set ylabel 'S(f) (Hz^{-1})'
plot  "PSD.dat" using 1:3 notitle with lines lc rgb "grey", "PSD.dat" using 1:2 notitle with lines lc rgb "black"
set origin 0.04, 0.0
set size 0.96, 0.57
set xtics (0,200,400,600,800,1000,1200,1400,1600,1800,2000)
set xlabel 'f (Hz)'
unset logscale y
unset ylabel
unset yrange
unset ytics
set ytics (0,1,2,3,4)
set format y '%.0f'
set yrange [0:2]
set ylabel 'AD'
plot "ADtest.dat" using 1:2 notitle with points pt 7 ps 0.5  lc rgb "blue", 0.752 notitle lc rgb "red"
unset multiplot
