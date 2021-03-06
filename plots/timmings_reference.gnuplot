set term term_type
set output output_file
set key below
set ylabel "Time (ns)"
set xlabel "N"

set logscale y
set logscale x 2
set format x "2^{%L}"
set format y "10^{%L}"

set grid mxtics mytics
set yrange [10**3:10**11]

plot 'data/timmings_reference.dat' using 1:2 w l title "Allocations",\
  '' u 1:3 w l title "Create array",\
  '' u 1:4 w l title "Sort array",\
  '' u 1:5 w l title "Fill matrix",\
  '' u 1:6 w l title "Fill tally",\
  '' u 1:7 w l title "Create histogram",\
  '' u 1:8 w l title "Total";
set output
