#grom=/gromacs/gromacs-2018.6/bin/gmx
grom=/gromacs/gromacs-5.1.5_pH_I/bin/gmx
Dir=`pwd`
mkdir -p $Dir/analysis/PMF/data

system="penta" # system name

zzmin=-1.2 #smallest position
zzmax=1.5  #largest position

min=5000   # first point to be read
max=250000  # last point to be read

zprof0=1.4  #value to use as zero reference
bootstrap=10 #

/bin/rm -f aux*

USs=`sed -n '/{/,/}/p' REUS_input.py | awk '$1!="USs" && $1!="}" {print $1}' | awk -F "'" '{printf("%4s", $2)}'`

rm -f $Dir/analysis/PMF/data/tpr-files.dat
for u in $USs
do
    if [ -f $Dir/US${u}/${system}_001.tpr.gz ]; then
	gunzip $Dir/US${u}/${system}_001.tpr.gz
    fi
    
    ls -1 $Dir/US${u}/${system}_001.tpr >> $Dir/analysis/PMF/data/tpr-files.dat
    
    cat  $Dir/US${u}/${system}_???_pullf.xvg | \
	awk -v min=$min -v max=$max '$1 > min && $1 <= max' > $Dir/analysis/PMF/data/aux_${u}_f.xvg
    ls -1 $Dir/analysis/PMF/data/aux_${u}_f.xvg >> $Dir/analysis/PMF/data/aux_f-files.dat
done
	
#cd $Dir
# Global PMFs
${grom} wham -it $Dir/analysis/PMF/data/tpr-files.dat \
	-if $Dir/analysis/PMF/data/aux_f-files.dat \
	-o $Dir/analysis/PMF/data/aux_profilef_pH${pH}.xvg \
	-hist $Dir/analysis/PMF/data/aux_histf_pH${pH}.xvg \
	-bsres $Dir/analysis/PMF/pmf_${system}.xvg \
	-bsprof $Dir/analysis/PMF/data/pmf_bsprofilesf_pH${pH}.xvg \
	-unit kCal -xvg none\
	-tol 1e-6 \
	-bins 100 \
	-zprof0 $zprof0 \
	-b $min \
	-e $max \
	-min $zzmin \
	-max $zzmax \
	-nBootstrap $bootstrap
	
rm -f $Dir/analysis/PMF/data/aux_* $Dir/analysis/PMF/data/\#* $Dir/analysis/PMF/data/*\~

####### Make the blue line with only the M-state #####

gpfile="$Dir/analysis/PMF/plot.gp"
cat <<EOF >> $gpfile
set term postscript enhanced color solid "Helvetica" 20
set encoding iso_8859_1
set border 31 lt -1 lw 1
set output "$Dir/analysis/PMF/plots_profiles.ps"

set style line  1 lt rgb "#FF0000" lw 3 pt 7 ps 1.0
set style line  2 lt rgb "#0000FF" lw 3 pt 7 ps 1.0
set style line  3 lt rgb "#000000" lw 3 pt 5 ps 1.0
set style line  4 lt rgb "#FF0000" lw 3 pt 7 ps 1.0 dt 2
set style line  5 lt rgb "#0000FF" lw 3 pt 7 ps 1.0 dt 2
set style line  6 lt rgb "black" lw 3 pt 7 ps 1.0
set style line  7 lt rgb "black" lw 3 pt 7 ps 1.0 dt 2

set grid back

set xlabel "Distance to monolayer P (nm)" font "Helvetica, 25"
set xrange [-2:+4] 
set xtics -2,1.0
set mxtics 5
set ylabel "Free energy  (kcal mol^{-1})" font "Helvetica, 25"
unset yrange
set yrange [-5:1]
set ytics -5,1
set mytics 2
set key left

set xzeroaxis

plot '$Dir/analysis/PMF/pmf_${system}.xvg'  u 1:2   title 'pH07.00'      w l        ls 6 ,\
     '$Dir/analysis/PMF/pmf_${system}.xvg'  u 1:2:3 every 5 notitle w errorbar ls 6 ,\

EOF

gnuplot $gpfile
! ps2pdf $Dir/analysis/PMF/plots_profiles.ps $Dir/analysis/PMF/pmf.pdf

rm -f $Dir/analysis/PMF/data/pull_*.xvg $Dir/analysis/PMF/data/*.ps $Dir/analysis/PMF/data/*0_seg-*.xvg $Dir/analysis/PMF/data/*bsprofilesf*  $Dir/analysis/PMF/data/tpr-files.dat $gpfile  $Dir/analysis/PMF/data/\#* $Dir/analysis/PMF/data/*\~



