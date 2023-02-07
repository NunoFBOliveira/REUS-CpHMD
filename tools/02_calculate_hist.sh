#grom=/gromacs/gromacs-2018.6/bin/gmx
grom=/gromacs/gromacs-5.1.5_pH_I/bin/gmx
Dir=`pwd` #"/home/noliveira/02_ADP_ATP_Transport/05_AAC+Sub"
mkdir -p analysis/histog/data

system="penta"

min=-1.2
max=1.5
bin=0.05

USs=`sed -n '/{/,/}/p' REUS_input.py | awk '$1!="USs" && $1!="}" {print $1}' | awk -F "'" '{printf("%4s", $2)}'`

for u in $USs
do
    echo "Doing US " $u
    awk '{print NR/10, $2}' $Dir/US${u}/*pullx.xvg > ${Dir}/analysis/histog/data/US${u}.dat

    awk '$1>=500 {print $2}' ${Dir}/analysis/histog/data/US${u}.dat | histog -r $min,$max -d $bin > ${Dir}/analysis/histog/data/histog_${u}.dat

    tot=`awk '{s+=$2}END{print s}' ${Dir}/analysis/histog/data/histog_${u}.dat `

    awk -v tot=$tot -v b=$bin '{print $1, ($2/tot)*b}' ${Dir}/analysis/histog/data/histog_${u}.dat > ${Dir}/analysis/histog/data/histog_${u}_norm.dat
	
    rm -f ${Dir}/analysis/histog/data/\#* ${Dir}/analysis/histog/data/*\~
done

gpfile="${Dir}/analysis/histog/plot_hist.gp"
cat <<EOF >> $gpfile
set term postscript enhanced color solid "Helvetica" 20
set encoding iso_8859_1
set border 31 lt -1 lw 1
set output "${Dir}/analysis/histog/plots_hist.ps"

set size 1.1,1.1
set lmargin 7
set tmargin 1
set bmargin 1

set style line  1 lt rgb "#FF0000" lw 3 pt 7 ps 1.0
set style line  2 lt rgb "#0000FF" lw 3 pt 7 ps 1.0
set style line  3 lt rgb "#000000" lw 3 pt 5 ps 1.0
set style line  4 lt rgb "#FF0000" lw 3 pt 7 ps 1.0
set style line  5 lt rgb "#0000FF" lw 3 pt 7 ps 1.0 dt 2
set style line  6 lt rgb "black" lw 3 pt 7 ps 1.0
set style line  7 lt rgb "black" lw 3 pt 7 ps 1.0 dt 2

 set linetype  1 lc rgb "dark-violet" lw 5
 set linetype  2 lc rgb "#009e73" lw 5
 set linetype  3 lc rgb "#56b4e9" lw 5
 set linetype  4 lc rgb "#e69f00" lw 5
 set linetype  5 lc rgb "#468c11" lw 5
 set linetype  6 lc rgb "#0072b2" lw 5
 set linetype  7 lc rgb "#e51e10" lw 5
 set linetype  8 lc rgb "black"   lw 5
 set linetype  9 lc rgb "gray50"  lw 5
 set linetype cycle  9

set grid back
set tics front

set xlabel "Distance (nm)" font "Helvetica, 25"
set xrange [-1.5:1.5] 
set xtics -2,0.5 font "Helvetica,22" nomirror
set mxtics 5 

set ylabel "Probability density" font "Helvetica, 25"
unset yrange
set yrange [0:0.02]
set ytics 0,0.002 font "Helvetica,22" nomirror
set mytics 4

unset key

#set xzeroaxis

filelist=system("ls ${Dir}/analysis/histog/data/histog_*_norm.dat")
plot for [filename in filelist] filename using 1:2 w filledcurve x1 notitle

unset arrow

EOF

gnuplot $gpfile
! ps2pdf ${Dir}/analysis/histog/plots_hist.ps ${Dir}/analysis/histog/hist.pdf

rm -f  *.ps  $gpfile 


