set term pdfcairo
set outp 'nulls2.pdf'

set pm3d map
set cbr [0:40]
set title 'number of iterations of Newton-Raphsdon method'
set xlabel 'Re(z)'
set ylabel 'Im(z)'
set size square
sp 'nullstellen.dat' us 2:3:1


unset outp
# pause -1

